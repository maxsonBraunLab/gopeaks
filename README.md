[![Go](https://github.com/maxsonBraunLab/gopeaks/actions/workflows/go.yml/badge.svg?branch=main)](https://github.com/maxsonBraunLab/gopeaks/actions/workflows/go.yml)

# README

Simple peak caller for CUT&TAG data

# Configure

Download the latest release binary specific for your operating system:

```
wget -O gopeaks https://github.com/maxsonBraunLab/gopeaks/releases/download/v0.1.5/gopeaks-linux-amd64
chmod +x gopeaks
```

# Example Usage

```
./gopeaks -h 
Usage of ./gopeaks:
  -bam string
        Bam file with 
  -control string
        Bam file with contriol signal to be subtracted
  -cs string
        Supply chromosome sizes for the alignment genome if not found in the bam header
  -mdist int
        Merge distance for nearby peaks (default 150)
  -minwidth int
        Minimum width to be considered a peak (default 250)
  -mr int
        Min reads per coverage bin to be considered (default 15)
  -of string
        Output file to write peaks to (default "peaks.bed")
  -pval float
        Pvalue threshold for keeping a peak bin (default 0.05)
  -slide int
        Slide size for coverage bins (default 50)
  -step int
        Bin size for coverage bins (default 100)
```

## call peaks on a bam file using an IgG control

```
./gopeaks -bam ../chr1.bam -control chr1_igg.bam -cs data/hg38.known.chrom.sizes -of chr1.bed
```
