package main

import (
	"testing"

	"github.com/matryer/is"
	gn "github.com/pbenner/gonetics"
)

func TestFilterUnknownChroms(t *testing.T) {
	is := is.New(t)
	chroms := []string{"ChrUn", "chrUn_test", "chrUn_GL000218v1", "chrM", "chr1_blah", "EBV", "chr3", "chr1", "chr4"}
	got := filterUnkownChroms(chroms)
	want := []string{"chr3", "chr1", "chr4"}
	is.Equal(got, want) // should equal
}

func TestBinGenome(t *testing.T) {
	is := is.New(t)
	seqnames := []string{"chr1", "chr2", "chr3", "chr4"}
	lengths := []int{1000, 2000, 5000, 300}
	g := gn.NewGenome(seqnames, lengths)
	bR := binGenome(g, 100, 50)

	is.Equal(bR.Seqnames[0], "chr1")                  // first should be chr1
	is.Equal(bR.Seqnames[len(bR.Seqnames)-1], "chr4") // last should be chr4

	// test that only the four chroms
	set := make(map[string]bool)
	for _, k := range bR.Seqnames {
		set[k] = true
	}
	is.Equal(len(set), 4) // should be only four chroms
}

func TestCountOverlaps(t *testing.T) {
	is := is.New(t)
	r1 := gn.NewGRanges(
		[]string{"chr1", "chr1", "chr1", "chr1"},
		[]int{1, 80, 90, 9},
		[]int{20, 100, 110, 40},
		[]byte{'*', '*', '*', '*'})

	r2 := gn.NewGRanges(
		[]string{"chr1", "chr1", "chr1"},
		[]int{4, 8, 90},
		[]int{8, 10, 95},
		[]byte{'*', '*', '*'})
	rCts := countOverlaps(r1, r2)
	dat := rCts.MetaData[0].([]int)
	is.Equal(len(dat), 4)
}

func TestMaxIntSlice(t *testing.T) {
	is := is.New(t)
	tslice := []int{1, 2, 3, 4, 5, 5, 7, 8, 0, 10, 300}
	want := float64(300)
	got := MaxIntSlice(tslice)
	is.Equal(got, want) // should equal 300
}

func TestMergeWithin(t *testing.T) {
	is := is.New(t)
	r1 := gn.NewGRanges(
		[]string{"chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2"},
		[]int{1, 80, 90, 9, 400, 550, 3, 30, 80},
		[]int{20, 100, 110, 40, 500, 600, 10, 40, 90},
		[]byte{'*', '*', '*', '*', '*', '*', '*', '*', '*'})
	want := gn.NewGRanges(
		[]string{"chr1", "chr1", "chr2"},
		[]int{1, 40, 3},
		[]int{400, 600, 90},
		[]byte{'*', '*', '*'})
	got := mergeWithin(r1, 150)
	is.True(got.Ranges[0].From == want.Ranges[0].From)
	is.True(got.Ranges[2].From == want.Ranges[2].From)
	is.True(got.Seqnames[2] == want.Seqnames[2])
}

func TestSubtractIntSlice(t *testing.T) {
	is := is.New(t)
	s1 := []int{1, 2, 3, 4, 5, 5, 6, 6, 7, 8, 8, 8, 9, 20}
	s2 := []int{1, 1.0, 2, 5, 4, 3, 2, 4, 2, 3, 3, 2, 2, 2}
	new := subtractIntSlices(s1, s2)
	is.Equal(len(new), len(s1))
}

func TestSumIntSlice(t *testing.T) {
	is := is.New(t)
	s1 := []int{1, 2, 3, 4}
	sum := sumIntSlice(s1)
	is.Equal(sum, 10)
}
