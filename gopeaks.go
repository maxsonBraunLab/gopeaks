package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"regexp"
	"strings"
	"sync"

	gn "github.com/pbenner/gonetics"
	"github.com/sirupsen/logrus"
	"gonum.org/v1/gonum/stat/distuv"
)

func main() {

	bam := flag.String("bam", "", "Bam file with ")
	cs := flag.String("cs", "", "Supply chromosome sizes for the alignment genome if not found in the bam header")
	within := flag.Int("mdist", 150, "Merge distance for nearby peaks")
	outfile := flag.String("of", "peaks.bed", "Output file to write peaks to")
	minreads := flag.Int("mr", 15, "Min reads per coverage bin to be considered")
	pval := flag.Float64("pval", 0.05, "Pvalue threshold for keeping a peak bin")
	step := flag.Int("step", 100, "Bin size for coverage bins")
	slide := flag.Int("slide", 50, "Slide size for coverage bins")
	minwidth := flag.Int("minwidth", 150, "Minimum width to be considered a peak")

	flag.Parse()

	// require bamfile
	if *bam == "" {
		flag.Usage()
		os.Exit(1)
	}

	// read bamfile to GRanges
	r := gn.GRanges{}
	if err := r.ImportBamPairedEnd(*bam); err != nil {
		logrus.Errorf("Error %s", err.Error())
		os.Exit(1)
	}

	g := gn.Genome{}
	var err error
	if *cs != "" {
		err := g.Import(*cs)
		if err != nil {
			logrus.Errorln("Failed to import chromsizes file")
			os.Exit(1)
		}
	}

	if *cs == "" {
		fmt.Println("Reading chromsizes from bam header...")
		g, err = gn.BamImportGenome(*bam)
		if err != nil {
			fmt.Println("Genome could not be determined from bam file")
		}
	}

	gf := KnownChroms(&g)
	fr := r.FilterGenome(gf)

	writeBigWig(fr, *bam)

	// calculate coverage
	binRanges := binGenome(g, *step, *slide)
	binCounts := countOverlaps(binRanges, fr)
	nreads := fr.Length()

	// callpeaks
	peaks := callpeaks(binCounts, float64(nreads), *within, *minwidth, *minreads, *pval)

	fmt.Printf("Number of peaks found: %d\n", peaks.Length())
	err = peaks.ExportBed6(*outfile, false)
	if err != nil {
		logrus.Errorln(err)
	}
}

func callpeaks(coverage gn.GRanges, total float64, within, width, minreads int, pval float64) gn.GRanges {

	covCounts := coverage.GetMetaInt("overlap_counts")
	// calculate coverage in non-zero bins
	sum := 0
	n := 0
	for i := 0; i < len(covCounts); i++ {
		if covCounts[i] != 0 {
			sum += covCounts[i]
			n += 1
		}
	}

	// probability of read in non-zero bin
	p := (float64(sum) / float64(n)) / float64(total)
	max := MaxIntSlice(covCounts)

	// calculate probability map
	max1p := max + 1
	probMap := map[int]float64{}
	for i := 0; float64(i) < max1p; i++ {
		probMap[i] = BinomTest(i, total, p)
	}

	// filter coverage bins
	var keepSlice []int
	for i := 0; i < len(covCounts); i++ {
		cnt := covCounts[i]
		if keep := filterBins(cnt, probMap, minreads, pval); keep {
			keepSlice = append(keepSlice, i)
		}
	}

	// merge overlapping and nearby peaks
	binsKeep := coverage.Subset(keepSlice)
	binsKeepMerge := binsKeep.Merge()
	peaks := mergeWithin(binsKeepMerge, within)
	peaksFilt := filterPeakWidth(peaks, width)
	return peaksFilt
}

// filterPeakWidth returns a granges object with ranges having width
// greater than the provided width
func filterPeakWidth(peaks gn.GRanges, width int) gn.GRanges {
	var keepIdx []int
	for i := 0; i < len(peaks.Seqnames); i++ {
		if (peaks.Ranges[i].To - peaks.Ranges[i].From) > width {
			keepIdx = append(keepIdx, i)
		}
	}
	return peaks.Subset(keepIdx)
}

// return true if bin is significant
func filterBins(c int, cMap map[int]float64, minReads int, threshold float64) bool {
	p := 1.0
	if c > minReads {
		p = cMap[c]
	}
	if p < threshold {
		return true
	}
	return false
}

// BinomTest returns the p-value testing the null hypothesis that the
// probability of a positive Bernoulli trial is p
func BinomTest(count int, total float64, p float64) float64 {
	dist := distuv.Binomial{N: float64(total) - float64(count), P: p}
	return dist.Prob(float64(count))
}

// MaxIntSlice returns the Max of an []Int
// cast as a float64
func MaxIntSlice(slice []int) float64 {
	max := 0
	for _, i := range slice {
		if max < i {
			max = i
		}
	}
	return float64(max)
}

// merges ranges in obj that are "within" base pairs apart
func mergeWithin(obj gn.GRanges, within int) gn.GRanges {

	out := []gn.Range{}
	outSeqs := []string{}

	in := obj.Ranges
	inSeqs := obj.Seqnames

	for i := 0; i < len(in); i++ {

		outLen := len(out)
		if i == 0 {
			out = append(out, in[i])
			outSeqs = append(outSeqs, inSeqs[i])
			continue
		}

		if outSeqs[len(outSeqs)-1] == inSeqs[i] {
			if (out[outLen-1].To + within) >= in[i].From {
				out[outLen-1].To = in[i].To
			} else {
				// append
				out = append(out, in[i])
				outSeqs = append(outSeqs, inSeqs[i])
			}
		} else {
			out = append(out, in[i])
			outSeqs = append(outSeqs, inSeqs[i])
		}
	}

	of := []int{}
	ot := []int{}
	os := []byte{}
	for _, r := range out {
		of = append(of, r.From)
		ot = append(ot, r.To)
		os = append(os, '*')
	}
	return gn.NewGRanges(outSeqs, of, ot, os)
}

// countOverlaps counts the overlapping in r2 and reports them as
// a new metadata column "overlap_counts" on r2
func countOverlaps(r1 gn.GRanges, r2 gn.GRanges) gn.GRanges {
	s, _ := gn.FindOverlaps(r1, r2)
	idxMap := map[int]int{}
	for i := 0; i < len(s); i++ {
		idxMap[s[i]] += 1
	}
	var olaps []int
	for i := 0; i < r1.Length(); i++ {
		var cnt int
		cnt, ok := idxMap[i]
		if !ok {
			cnt = 0
		}
		olaps = append(olaps, cnt)
	}
	r1.AddMeta("overlap_counts", olaps)
	return r1
}

func binChrom(genome gn.Genome, chr string, step, slide int) gn.GRanges {
	var seqnames []string
	var ranges []gn.Range
	var strand []byte
	start := 0
	len, _ := genome.SeqLength(chr)
	count := 0
	for start <= len-step {
		end := start + step
		ranges = append(ranges, gn.Range{From: start, To: end})
		seqnames = append(seqnames, chr)
		start += slide
		count += 1
	}

	strand = make([]byte, count)
	for i := 0; i < count; i++ {
		strand[i] = '*'
	}

	ret := gn.GRanges{
		Seqnames: seqnames,
		Ranges:   ranges,
		Strand:   strand,
		Meta:     gn.Meta{},
	}

	return ret
}

// bin Result stores the chromosome bin result
// and it's chromosome sort order
type BinnedRangesOrder struct {
	Order  int
	Ranges gn.GRanges
}

// read results into output channel
func binChromToChan(g gn.Genome, chr string, out chan BinnedRangesOrder, step, slide int) {
	var res BinnedRangesOrder
	for i, s := range g.Seqnames {
		if s == chr {
			res.Order = i
			res.Ranges = binChrom(g, chr, step, slide)
			out <- res
		}
	}
}

func handleChromBins(input chan BinnedRangesOrder, output chan gn.GRanges, wg *sync.WaitGroup) {

	// parse input into slice
	var gRes []BinnedRangesOrder
	for r := range input {
		gRes = append(gRes, r)
		wg.Done()
	}

	var ret gn.GRanges
	// append to output preserving chr order
	for i := 0; i < len(gRes); i++ {
		for _, g := range gRes {
			if g.Order == i {
				ret = ret.Append(g.Ranges)
			}
		}
	}
	output <- ret
}

// bin genome into overlapping ranges with step and slide
// bin genome creates coverages for each chromosome in separate go routines
func binGenome(genome gn.Genome, step int, slide int) gn.GRanges {
	input := make(chan BinnedRangesOrder)
	output := make(chan gn.GRanges)
	var wg sync.WaitGroup
	go handleChromBins(input, output, &wg)
	defer close(output)
	for _, chr := range genome.Seqnames {
		wg.Add(1)
		go binChromToChan(genome, chr, input, step, slide)
	}

	wg.Wait()
	close(input)
	return <-output
}

// filters unknown chromosome names from a strings slice
// using a regex of unwanted string matches
func filterUnkownChroms(start []string) []string {
	var ret []string
	filt := `Un|_|EBV|N|M`
	for _, s := range start {
		r := regexp.MustCompile(filt)
		if !r.MatchString(s) {
			ret = append(ret, s)
		}
	}
	return ret
}

// returns a genome of filtered chromosomes
func KnownChroms(genome *gn.Genome) gn.Genome {

	// make map of known seqs
	knownMap := map[string]bool{}
	knownSeqs := filterUnkownChroms(genome.Seqnames)
	for _, s := range knownSeqs {
		knownMap[s] = true
	}

	// return new genome with only known chroms
	seqnames := []string{}
	lengths := []int{}
	for i := 0; i < genome.Length(); i++ {
		if b, _ := knownMap[genome.Seqnames[i]]; b {
			seqnames = append(seqnames, genome.Seqnames[i])
			lengths = append(lengths, genome.Lengths[i])
		}
	}
	return gn.NewGenome(seqnames, lengths)
}

func makeTrackName(infile string) string {
	return strings.Replace(infile, ".bam", ".bw", 1)
}

func writeBigWig(bamRanges gn.GRanges, bamfilename string) {
	outfilename := makeTrackName(bamfilename)
	var filenamesTreatment []string
	var filenamesControl []string
	var fraglenTreatment []int
	var fraglenControl []int

	optionsList := []interface{}{
		gn.OptionEstimateFraglen{Value: true},
		gn.OptionNormalizeTrack{Value: "cpm"},
	}
	filenamesTreatment = append(filenamesTreatment, bamfilename)
	result, _, _, _ := gn.BamCoverage(outfilename, filenamesTreatment, filenamesControl, fraglenTreatment, fraglenControl, optionsList...)
	fmt.Printf("Writing track `%s'... ", outfilename)
	parameters := gn.DefaultBigWigParameters()
	if err := (gn.GenericTrack{Track: result}).ExportBigWig(outfilename, parameters); err != nil {
		log.Fatal(err)
	} else {
		fmt.Println("done")
	}

}
