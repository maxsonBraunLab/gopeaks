package main

import (
	"flag"
	"fmt"
	"os"
	"regexp"

	gn "github.com/pbenner/gonetics"
	"github.com/sirupsen/logrus"
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
		logrus.Infoln("Attempting to read chromsizes from bam header...")
		g, err = gn.BamImportGenome(*bam)
		if err != nil {
			logrus.Infoln("Genome could not be determined from bam file")
		}
	}

	gf := KnownChroms(&g)
	r.FilterGenome(gf)

	logrus.Infoln("Total Coverage: %s", r.Length())
	logrus.Infoln("Total Coverage: %s", r.Length())

	fmt.Println(r.Row(1))

	// binRanges := binGenome(gf)
	// binCounts := countOverlaps(binRanges, r)

	// allCts := binCounts.MetaData

}

func countOverlaps(r1 gn.GRanges, r2 gn.GRanges) gn.GRanges {
	s, _ := gn.FindOverlaps(r1, r2)
	idxMap := map[int]int{}
	for i := 0; i < len(s); i++ {
		idxMap[s[i]] += 1
	}
	olaps := []int{}
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
