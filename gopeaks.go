package main

import (
	"flag"
	"fmt"
	"os"
	"regexp"

	gn "github.com/pbenner/gonetics"
	"github.com/sirupsen/logrus"
	"gonum.org/v1/gonum/stat/distuv"
)

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

	// calculate coverage in non-zer bins
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
	fmt.Printf("Max counts in bin: %f", max)

	max1p := max + 1
	probMap := map[int]float64{}
	for i := 0; float64(i) < max1p; i++ {
		bt := BinomTest(i, total, p)
		probMap[i] = bt
	}

	var keepSlice []int
	for i := 0; i < len(overlaps); i++ {
		cnt := overlaps[i]
		if keep := filterBins(cnt, probMap, 15, 0.05); keep {
			keepSlice = append(keepSlice, i)
		}
	}

	binsKeep := fr.Subset(keepSlice)
	fmt.Println(binsKeep.Length())
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
		if max < slice[i] {
			max = slice[i]
		}
	}
	return float64(max)
}

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

func binGenome(genome gn.Genome) gn.GRanges {
	var seqnames []string
	var ranges []gn.Range
	var strand []byte
	step := 100
	slide := 50

	start := 0
	name := "chr1"
	for start <= genome.Lengths[0]-step {
		end := start + step
		ranges = append(ranges, gn.Range{From: start, To: end})
		seqnames = append(seqnames, name)
		start += slide
	}
	strand = make([]byte, len(seqnames))
	for i := 0; i < len(seqnames); i++ {
		strand[i] = '*'
	}
	ret := gn.GRanges{Seqnames: seqnames,
		Ranges: ranges,
		Strand: strand,
		Meta:   gn.Meta{}}

	return ret
}

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
