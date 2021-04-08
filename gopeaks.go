package main

import (
	"flag"
	"fmt"

	"github.com/pbenner/gonetics"
)

func main() {
	bam := flag.String("bamfile", "", "Bam file with ")
	flag.Parse()

	r := gonetics.GRanges{}
	if err := r.ImportBamPairedEnd(*bam); err != nil {
		fmt.Errorf("Error %s", err.Error())
	}

	fmt.Println(r)
	fmt.Println(r.Length())
}
