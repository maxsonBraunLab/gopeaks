package main

import (
	"flag"
	"fmt"
	"os"
	"regexp"
	"sync"

	gn "github.com/pbenner/gonetics"
	"github.com/sirupsen/logrus"
	"gonum.org/v1/gonum/stat/distuv"
)

// number of base-pairs to clsuter nearby peaks
const clusterWithin = 150

func main() {

	bam := flag.String("bam", "", "Bam file with ")
	cs := flag.String("cs", "", "Supply chromosome sizes for the alignment genome if not found in the bam header")

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
		}
	}

	if *cs == "" {
		logrus.Infoln("Reading chromsizes from bam header...")
		g, err = gn.BamImportGenome(*bam)
		if err != nil {
			logrus.Infoln("Genome could not be determined from bam file")
		}
	}

	gf := KnownChroms(&g)
	fr := r.FilterGenome(gf)

	total := float64(fr.Length())

	// calculating coverage
	binRanges := binGenome(gf)
	binCounts := countOverlaps(binRanges, r)

	overlaps := binCounts.GetMetaInt("overlap_counts")

	// calculate coverage in non-zero bins
	sum := 0
	n := 0
	for i := 0; i < len(overlaps); i++ {
		if overlaps[i] != 0 {
			sum += overlaps[i]
			n += 1
		}
	}

	// probability of read in non-zero bin
	p := (float64(sum) / float64(n)) / float64(total)

	max := MaxIntSlice(overlaps)

	fmt.Printf("Total Reads: %f\n", total)
	fmt.Printf("Total Coverage: %d\n", sum)
	fmt.Printf("Probability of read in non-zero bin: %e\n", p)
	fmt.Printf("Max counts in bin: %f\n", max)

	max1p := max + 1
	probMap := map[int]float64{}
	for i := 0; float64(i) < max1p; i++ {
		bt := BinomTest(i, total, p)
		probMap[i] = bt
	}

	var keepSlice []int
	fmt.Printf("The number of overlaps: %d\n", len(overlaps))
	for i := 0; i < len(overlaps); i++ {
		cnt := overlaps[i]
		if keep := filterBins(cnt, probMap, 15, 0.05); keep {
			keepSlice = append(keepSlice, i)
		}
	}

	binsKeep := binCounts.Subset(keepSlice)
	binsKeepMerge := binsKeep.Merge()
	peaks := mergeWithin(binsKeepMerge, 150)
	fmt.Printf("Numbe of peaks founds: %d\n", peaks.Length())
	err = peaks.ExportBed6("peaks.bed", false)
	if err != nil {
		logrus.Errorln(err)
	}
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
			if (out[outLen-1].To + 150) >= in[i].From {
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

// appendAll apppends a slice of ranges sequentially
func appendAll(ranges []gn.GRanges) gn.GRanges {
	ret := ranges[0]
	for _, r := range ranges[1:] {
		ret = ret.Append(r)
	}
	return ret
}

// shift slice down, popping of the last element
// and using zero as the first
func shiftSlice(sl []int) []int {
	var shsl []int
	size := len(sl)
	shsl = append(shsl, 0)
	shsl = append(shsl, sl[:size-1]...)
	return shsl
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

func binChrom(genome gn.Genome, chr string) gn.GRanges {
	var seqnames []string
	var ranges []gn.Range
	var strand []byte
	step := 100
	slide := 50
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
func binChromToChan(g gn.Genome, chr string, out chan BinnedRangesOrder) {
	var res BinnedRangesOrder
	for i, s := range g.Seqnames {
		if s == chr {
			res.Order = i
			res.Ranges = binChrom(g, chr)
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
// then merges them
// TODO: parametarize step and slide
func binGenome(genome gn.Genome) gn.GRanges {
	input := make(chan BinnedRangesOrder)
	output := make(chan gn.GRanges)
	var wg sync.WaitGroup
	go handleChromBins(input, output, &wg)
	defer close(output)
	for _, chr := range genome.Seqnames {
		wg.Add(1)
		go binChromToChan(genome, chr, input)
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
