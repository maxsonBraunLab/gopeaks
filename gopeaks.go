package main

import (
	"encoding/json"
	"fmt"
	"os"
	"regexp"
	"sort"
	"sync"
	"time"
	"math"
	"strings"

	"github.com/akamensky/argparse"
	gn "github.com/pbenner/gonetics"
	"github.com/sirupsen/logrus"
	"gonum.org/v1/gonum/stat/distuv"
	"github.com/go-gota/gota/dataframe"
	"github.com/go-gota/gota/series"
)

const gopeaks_version = "1.0.0"

type Metrics struct {
	Version string `json:"gopeaks_version"`
	Date    string `json:"date"`
	Elapsed string `json:"elapsed"`
	Prefix  string `json:"prefix"`
	Command string `json:"command"`
	Peaks   int    `json:"peak_counts"`
}

func (m *Metrics) Log(op string) {
	resp, err := json.MarshalIndent(m, "", "\t")
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}

	f, err := os.Create(op + "_gopeaks.json")
	defer f.Close()
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}

	f.WriteString(string(resp))
	f.WriteString("\n")
}

func main() {

	// start time is what elapsed metric
	// is calculated from
	startTime := time.Now()

	parser := argparse.NewParser("GoPeaks",`GoPeaks is a peak caller designed for CUT&TAG/CUT&RUN sequencing data. GoPeaks by default works best with narrow peaks such as H3K4me3 and transcription factors. GoPeaks can be used with the "--broad" flag to call broad peaks like H3K27Ac/H3K4me1. We encourage users to explore the parameters of GoPeaks to analyze their data.`)
	bam := parser.String("b", "bam", &argparse.Options{Help: "Input BAM file (must be paired-end reads)"})
	control := parser.String("c", "control", &argparse.Options{Help: "Input BAM file with control signal to be normalized (e.g. IgG, Input)"})
	cs := parser.String("s", "chromsize", &argparse.Options{Help: "Chromosome sizes for the genome if not found in the bam header"})
	within := parser.Int("m", "mdist", &argparse.Options{Help: "Merge peaks within <mdist> base pairs", Default: 1000})
	minreads := parser.Int("r", "minreads", &argparse.Options{Help: "Test genome bins with at least <minreads> read pairs.", Default: 15})
	pval := parser.Float("p", "pval", &argparse.Options{Help: "Define significance threshold <pval> with multiple hypothesis correction via Benjamini-Hochberg", Default: 0.05})
	step := parser.Int("t", "step", &argparse.Options{Help: "Bin size for coverage bins", Default: 100})
	slide := parser.Int("l", "slide", &argparse.Options{Help: "Slide size for coverage bins", Default: 50})
	minwidth := parser.Int("w", "minwidth", &argparse.Options{Help: "Minimum width (bp) of a peak", Default: 150})
	outprefix := parser.String("o", "prefix", &argparse.Options{Help: "Output prefix to write peaks and metrics file", Default: "sample"})
	version := parser.Flag("v", "version", &argparse.Options{Help: "Print the current GoPeaks version"})
	broad := parser.Flag("", "broad", &argparse.Options{Help: "Run GoPeaks on broad marks (--step 5000 & --slide 1000)"})
	verbose := parser.Flag("", "verbose", &argparse.Options{Help: "Run GoPeaks in verbose mode."})
	// note: "Required" interface clashes with --version flag.
	err := parser.Parse(os.Args)

	// parse flags --------------------------------------------------------------------------------

	// check version
	if *version == true {
		fmt.Println("GoPeaks version:", gopeaks_version)
		os.Exit(0)
	}

	// check argparse errors
	if err != nil {
		fmt.Print(parser.Usage(err))
		os.Exit(1)
	}

	// require args
	if *bam == "" {
		fmt.Println(parser.Help(nil))
		os.Exit(1)
	}

	// read bamfile to GRanges
	r := gn.GRanges{}
	if err := r.ImportBamPairedEnd(*bam, gn.BamReaderOptions{ReadName: false, ReadCigar: false, ReadSequence: false}); err != nil {
		logrus.Errorf("Error %s", err.Error())
		os.Exit(1)
	}

	g := gn.Genome{}
	if *cs != "" {
		err := g.Import(*cs)
		if err != nil {
			logrus.Errorln("Failed to import chromsizes file")
			os.Exit(1)
		}
	}

	if *cs == "" {
		if *verbose {
			fmt.Println("Reading chromsizes from bam header...")
		}
		g, err = gn.BamImportGenome(*bam)
		if err != nil {
			fmt.Println("Genome could not be determined from bam file")
			os.Exit(1)
		}
	}

	if *broad == true {
		x := 5000
		step = &x
		y := 1000
		slide = &y
	}

	// import data --------------------------------------------------------------------------------

	gf := KnownChroms(&g)
	fr := r.FilterGenome(gf)
	nreads := fr.Length()
	
	// Clear full BAM data from memory - we'll process by chromosome
	r = gn.GRanges{}

	// Load control data if provided
	var ctrlRanges *gn.GRanges
	var ctrlReads int
	if *control != "" {
		c := gn.GRanges{}
		if err := c.ImportBamPairedEnd(*control, gn.BamReaderOptions{ReadName: false, ReadCigar: false, ReadSequence: false}); err != nil {
			logrus.Errorf("Error %s", err.Error())
			os.Exit(1)
		}
		cr := c.FilterGenome(gf)
		ctrlReads = cr.Length()
		ctrlRanges = &cr
		// Clear full control BAM data - we'll process by chromosome
		c = gn.GRanges{}
	}

	// Two-pass approach: 
	// Pass 1: Calculate coverage per chromosome (process separately to save memory)
	// Pass 2: Calculate global stats and call peaks per chromosome
	
	if *verbose {
		fmt.Printf("Pass 1: Calculating coverage for %d chromosomes...\n", len(gf.Seqnames))
	}

	// Store chromosome-specific bin counts (only metadata, not full GRanges)
	type chrBinData struct {
		chr       string
		binRanges gn.GRanges
		counts    []float64
	}
	var chrBinDataList []chrBinData

	// Pass 1: Calculate coverage per chromosome
	for chrIdx, chr := range gf.Seqnames {
		if *verbose {
			fmt.Printf("  Processing chromosome %s (%d/%d)...\n", chr, chrIdx+1, len(gf.Seqnames))
		}

		// Filter reads for this chromosome only
		chrReads := filterGRangesByChrom(fr, chr)
		
		// Skip if no reads for this chromosome
		if chrReads.Length() == 0 {
			chrReads = gn.GRanges{}
			continue
		}
		
		// Create bins for this chromosome only
		chrBinRanges := binChrom(gf, chr, *step, *slide)
		
		// Skip if no bins for this chromosome
		if chrBinRanges.Length() == 0 {
			chrReads = gn.GRanges{}
			chrBinRanges = gn.GRanges{}
			continue
		}
		
		// Count overlaps for this chromosome
		chrBinCounts := countOverlaps(chrBinRanges, chrReads)
		
		// Process control if provided
		if ctrlRanges != nil {
			chrCtrlReads := filterGRangesByChrom(*ctrlRanges, chr)
			if chrCtrlReads.Length() > 0 {
				chrCtrlCounts := countOverlaps(chrBinRanges, chrCtrlReads)
				chrBinCounts = normalizeToControl(chrBinCounts, chrCtrlCounts, nreads, ctrlReads)
				chrCtrlReads = gn.GRanges{} // Clear memory
				chrCtrlCounts = gn.GRanges{} // Clear memory
			}
		}
		
		// Extract counts and store (keep bin ranges for later peak calling)
		counts := chrBinCounts.GetMeta("overlap_counts").([]float64)
		chrBinDataList = append(chrBinDataList, chrBinData{
			chr:       chr,
			binRanges: chrBinRanges,
			counts:    counts,
		})
		
		// Clear chromosome-specific data (keep counts in chrBinDataList)
		chrReads = gn.GRanges{}
		chrBinCounts = gn.GRanges{}
	}
	
	// Clear full read data - we have all counts now
	fr = gn.GRanges{}
	if ctrlRanges != nil {
		*ctrlRanges = gn.GRanges{}
	}

	// Calculate global statistics from all chromosome counts
	if *verbose {
		fmt.Println("Pass 2: Calculating global statistics and calling peaks...")
	}
	
	var allCounts []float64
	for _, data := range chrBinDataList {
		allCounts = append(allCounts, data.counts...)
	}
	
	// Calculate global binomial parameters
	nzSignals, nzBins, globalNTests := binomialParameters(allCounts, *minreads)
	p := (float64(nzSignals) / float64(nzBins)) / float64(nreads)
	
	if *verbose {
		fmt.Printf("Global stats: nzSignals=%.0f, nzBins=%d, nTests=%d, p=%.6e\n", nzSignals, nzBins, globalNTests, p)
	}
	
	// Pass 2: Calculate p-values for all bins globally, then do FDR correction
	// Collect all eligible bins with their p-values
	type globalBinData struct {
		chrIdx   int
		binIdx   int
		count    float64
		pval     float64
		binRange gn.Range
		chr      string
	}
	
	var globalBins []globalBinData
	
	for chrIdx, data := range chrBinDataList {
		for binIdx := 0; binIdx < len(data.counts); binIdx++ {
			cnt := data.counts[binIdx]
			if cnt > float64(*minreads) {
				prob := BinomTest(cnt, float64(nreads), p)
				globalBins = append(globalBins, globalBinData{
					chrIdx:   chrIdx,
					binIdx:   binIdx,
					count:    cnt,
					pval:     prob,
					binRange: data.binRanges.Ranges[binIdx],
					chr:      data.chr,
				})
			}
		}
	}
	
	if *verbose {
		fmt.Printf("Calculated p-values for %d bins across all chromosomes\n", len(globalBins))
	}
	
	// Extract data for FDR correction
	bins := make([]int, len(globalBins))
	counts := make([]float64, len(globalBins))
	pvals := make([]float64, len(globalBins))
	for i, gb := range globalBins {
		bins[i] = i // Store index into globalBins
		counts[i] = gb.count
		pvals[i] = gb.pval
	}
	
	// Global FDR correction
	keepIndices := filterBinsbyFDR(bins, counts, pvals, *pval, globalNTests, *outprefix)
	
	if *verbose {
		fmt.Printf("After FDR correction: %d significant bins\n", len(keepIndices))
	}
	
	// Build GRanges from significant bins (grouped by chromosome)
	peaksByChr := make(map[string][]gn.Range)
	for _, idx := range keepIndices {
		gb := globalBins[idx]
		peaksByChr[gb.chr] = append(peaksByChr[gb.chr], gb.binRange)
	}
	
	// Clear global bins data
	globalBins = nil
	bins = nil
	counts = nil
	pvals = nil
	
	// Convert to GRanges and merge peaks per chromosome
	var allPeaks gn.GRanges
	allPeaks = gn.GRanges{
		Seqnames: []string{},
		Ranges:   []gn.Range{},
		Strand:   []byte{},
		Meta:     gn.Meta{},
	}
	
	for chrIdx, data := range chrBinDataList {
		if *verbose && len(peaksByChr[data.chr]) > 0 {
			fmt.Printf("  Merging peaks for chromosome %s (%d/%d)...\n", data.chr, chrIdx+1, len(chrBinDataList))
		}
		
		ranges := peaksByChr[data.chr]
		if len(ranges) == 0 {
			continue
		}
		
		// Create GRanges for this chromosome's peaks
		seqnames := make([]string, len(ranges))
		starts := make([]int, len(ranges))
		ends := make([]int, len(ranges))
		strand := make([]byte, len(ranges))
		for i, r := range ranges {
			seqnames[i] = data.chr
			starts[i] = r.From
			ends[i] = r.To
			strand[i] = '*'
		}
		
		chrPeaks := gn.NewGRanges(seqnames, starts, ends, strand)
		
		// Merge overlapping and nearby peaks
		chrPeaksMerge := chrPeaks.Merge()
		chrPeaksMerged := mergeWithin(chrPeaksMerge, *within)
		chrPeaksFilt := filterPeakWidth(chrPeaksMerged, *minwidth)
		
		// Append to master list
		if chrPeaksFilt.Length() > 0 {
			allPeaks = allPeaks.Append(chrPeaksFilt)
		}
		
		// Clear chromosome-specific data
		chrPeaks = gn.GRanges{}
		chrPeaksMerge = gn.GRanges{}
		chrPeaksMerged = gn.GRanges{}
		chrPeaksFilt = gn.GRanges{}
	}

	peaks := allPeaks

	outfile := *outprefix + "_peaks.bed"
	err = peaks.ExportBed3(outfile, false)
	if err != nil {
		logrus.Errorln(err)
	}

	// write output metrics -----------------------------------------------------------------------
	metrics := &Metrics{
		Version: gopeaks_version,
		Date:    time.Now().Format("2006-01-02 3:4:5 PM"),
		Elapsed: time.Since(startTime).String(),
		Prefix:  *outprefix,
		Command:  strings.Join(os.Args, " "),
		Peaks:   peaks.Length(),
	}

	// log metrics to file
	metrics.Log(*outprefix)
}

func scaleTreatToControl(counts []float64, s1 []float64, s2 []float64) []float64 {
	scale := make([]float64, len(s1))
	d1map := map[int]float64{}
	for i, s := range s1 {
		d1map[i] = s
	}
	var frac float64
	for i, o := range s2 {
		if d1map[i] > 0 {
			frac = o / d1map[i]
			if frac > 1 {
				frac = 1
			}
		} else {
			frac = 1
		}
		scale[i] = math.Floor(counts[i] * (1 - frac))
	}
	return scale
}

func cpm(in []float64, nreads float64) []float64 {
	var cpm []float64
	for _, o := range in {
		num := o * (1e6 / nreads)
		cpm = append(cpm, num)
	}
	return cpm
}

func normalizeToControl(treat gn.GRanges, ctrl gn.GRanges, treads, creads int) gn.GRanges {
	tcounts := treat.GetMeta("overlap_counts").([]float64)
	ccounts := ctrl.GetMeta("overlap_counts").([]float64)

	// calculate the cpm for each track
	tcountsNorm := cpm(tcounts, float64(treads))
	ccountsNorm := cpm(ccounts, float64(creads))

	// scale the treatment
	// scaled_counts = treat * 1-(control/treat)
	// NOTE: intervals of 0 signal includes actual 0 bins PLUS where IgG > treatment thx to scaleTreatToControl
	scale := scaleTreatToControl(tcounts, tcountsNorm, ccountsNorm)
	treat.AddMeta("overlap_counts", scale)
	return treat
}

func binomialParameters(counts []float64, minreads int) (float64, int, int) {

	// nzSignals = total signal in non-zero bins
	// nzBins = number of non-zero bins
	// nTests = number of tests (binCounts > minreads)

	nzSignals := 0.0
	nzBins := 0
	nTests := 0

	for i := 0; i < len(counts); i++ {
		binCounts := float64(counts[i])
		// a bin can satisfy non-zero signal AND > minreads. This is okay.
		if binCounts != 0.0 {
			nzSignals += counts[i]
			nzBins += 1
		}
		if binCounts > float64(minreads) {
			nTests += 1
		}
	}

	return nzSignals, nzBins, nTests
}

func callpeaks(coverage gn.GRanges, total float64, within, width, minreads int, pval float64, outprefix string, verbose bool) gn.GRanges {

	// coverage = GRanges of overlap counts in a bin
	// total = total number of paired-end reads

	ccts := coverage.GetMeta("overlap_counts").([]float64)
	nzSignals, nzBins, nTests := binomialParameters(ccts, minreads)

	// calculate probability of read in non-zero bin
	p := (float64(nzSignals) / float64(nzBins)) / float64(total)

	if verbose {
		fmt.Println("nTests:", nTests)
		fmt.Println("nzSignals:", nzSignals)
		fmt.Println("nzBins:", nzBins)
		fmt.Println("n:", total)
		fmt.Println("p:", p)
		fmt.Println("mu:", float64(nzBins) * float64(p))
		fmt.Println("var:", float64(nzBins) * float64(p) * (1-float64(p)))
	}

	var keepSlice []int
	var bins []int
	var counts []float64
	var pvals []float64

	nTests = 0
	for i := 0; i < len(ccts); i++ {
		cnt := ccts[i]
		if cnt > float64(minreads) {
			prob := BinomTest(cnt, total, p)
			nTests += 1
			bins = append(bins, i)
			counts = append(counts, cnt)
			pvals = append(pvals, prob)
		}
	}

	// `pvals` is list of p-values per eligible bin. `pval` is threshold for significance.
	keepSlice = filterBinsbyFDR(bins, counts, pvals, pval, nTests, outprefix)

	// merge overlapping and nearby peaks -----------------------------------------------
	binsKeep := coverage.Subset(keepSlice)
	binsKeepMerge := binsKeep.Merge()
	peaks := mergeWithin(binsKeepMerge, within)
	peaksFilt := filterPeakWidth(peaks, width)
	return peaksFilt
}

// callpeaksWithGlobalStats is similar to callpeaks but uses pre-calculated global probability and global nTests
func callpeaksWithGlobalStats(coverage gn.GRanges, total float64, globalP float64, globalNTests int, within, width, minreads int, pval float64, outprefix string, verbose bool) gn.GRanges {

	ccts := coverage.GetMeta("overlap_counts").([]float64)

	var keepSlice []int
	var bins []int
	var counts []float64
	var pvals []float64

	for i := 0; i < len(ccts); i++ {
		cnt := ccts[i]
		if cnt > float64(minreads) {
			prob := BinomTest(cnt, total, globalP)
			bins = append(bins, i)
			counts = append(counts, cnt)
			pvals = append(pvals, prob)
		}
	}

	// `pvals` is list of p-values per eligible bin. `pval` is threshold for significance.
	// Use global nTests for FDR correction
	keepSlice = filterBinsbyFDR(bins, counts, pvals, pval, globalNTests, outprefix)

	// merge overlapping and nearby peaks -----------------------------------------------
	binsKeep := coverage.Subset(keepSlice)
	binsKeepMerge := binsKeep.Merge()
	peaks := mergeWithin(binsKeepMerge, within)
	peaksFilt := filterPeakWidth(peaks, width)
	return peaksFilt
}

// BinPvalPair stores bin index and p-value for sorting
type BinPvalPair struct {
	Bin  int
	Pval float64
}

func filterBinsbyFDR(Bins []int, Counts []float64, Pvals []float64, Threshold float64, Tests int, outprefix string) []int {

	keepBins := []int{}

	// Use a more memory-efficient approach: sort indices instead of full DataFrame
	// Create pairs for sorting
	pairs := make([]BinPvalPair, len(Bins))
	for i := range Bins {
		pairs[i] = BinPvalPair{Bin: Bins[i], Pval: Pvals[i]}
	}

	// Sort by p-value (ascending) using Go's sort package
	sort.Slice(pairs, func(i, j int) bool {
		return pairs[i].Pval < pairs[j].Pval
	})

	// Calculate ranks (accounting for ties)
	pvalMap := make(map[float64]int)
	rank := 0
	for _, pair := range pairs {
		if _, ok := pvalMap[pair.Pval]; !ok {
			rank += 1
			pvalMap[pair.Pval] = rank
		}
	}

	// Calculate adjusted p-values and filter
	for _, pair := range pairs {
		r := float64(pvalMap[pair.Pval])
		padj := float64(Tests) * pair.Pval / r
		if padj >= 1 {
			padj = 1
		}
		if padj < Threshold {
			keepBins = append(keepBins, pair.Bin)
		}
	}

	return keepBins
}

func assignRanks(fdrDF dataframe.DataFrame) dataframe.DataFrame {

	// implement smart ranking scheme to account for same pvals

	// assume the pval col is sorted numerically
	rank := 0
	rankSeries := series.New([]int{}, series.Int, "rankSeries")
	pvalMap := make(map[float64]int)

	// create pval:rank map
	for i := 0; i < fdrDF.Nrow(); i++ {
		pval := fdrDF.Elem(i, 2).Float()
		_, ok := pvalMap[pval] // output = value, bool
		if !ok {
			rank += 1
			pvalMap[pval] = rank
		}
	}

	// assign rank to rankSeries
	for i := 0; i < fdrDF.Nrow(); i++ {
		pval := fdrDF.Elem(i, 2).Float()
		rankSeries.Append(pvalMap[pval])
	}

	// add rank column to DF
	fdrDF = fdrDF.Mutate(series.New(rankSeries, series.Int, "rank"))

	return fdrDF

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

// BinomTest returns the p-value testing the null hypothesis that the
// probability of a positive Bernoulli trial of probability p is p
func BinomTest(count float64, total float64, p float64) float64 {
	// dev notes: may need to use one-tailed binomial test. we're not interested in bins < expected.
	dist := distuv.Binomial{N: float64(total) - count, P: p}
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
	idxMap := map[int]float64{}
	for i := 0; i < len(s); i++ {
		idxMap[s[i]] += 1
	}
	var olaps []float64
	for i := 0; i < r1.Length(); i++ {
		var cnt float64
		cnt, ok := idxMap[i]
		if !ok {
			cnt = 0.0
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

// filterGRangesByChrom filters GRanges to only include entries for a specific chromosome
func filterGRangesByChrom(gr gn.GRanges, chr string) gn.GRanges {
	var keepIdx []int
	for i := 0; i < len(gr.Seqnames); i++ {
		if gr.Seqnames[i] == chr {
			keepIdx = append(keepIdx, i)
		}
	}
	if len(keepIdx) == 0 {
		return gn.GRanges{
			Seqnames: []string{},
			Ranges:   []gn.Range{},
			Strand:   []byte{},
			Meta:     gn.Meta{},
		}
	}
	return gr.Subset(keepIdx)
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
