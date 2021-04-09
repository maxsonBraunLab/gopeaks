package main

import (
	"fmt"
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
	seqnames := []string{"chr1", "chr2"}
	lengths := []int{1000, 3000000}
	g := gn.NewGenome(seqnames, lengths)
	bR := binGenome(g)
	is.Equal(len(bR.Seqnames), 19) // should be 19
	is.Equal(bR.Ranges[0].From, 0) // should start at 0
	is.Equal(bR.Ranges[0].From, 0) // should start at 0
}

func TestCountOverlaps(t *testing.T) {
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
	fmt.Println(rCts.MetaData...)
}

func TestMaxIntSlice(t *testing.T) {
	is := is.New(t)
	tslice := []int{1, 2, 3, 4, 5, 5, 7, 8, 0, 10, 300}
	want := float64(300)
	got := MaxIntSlice(tslice)
	is.Equal(got, want) // should equal 300
}
