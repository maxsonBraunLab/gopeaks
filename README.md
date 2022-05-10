[![Go](https://github.com/maxsonBraunLab/gopeaks/actions/workflows/go.yml/badge.svg?branch=main)](https://github.com/maxsonBraunLab/gopeaks/actions/workflows/go.yml) ![conda](https://anaconda.org/jakevc/gopeaks/badges/installer/conda.svg)

# README

Simple peak caller for CUT&TAG data

# Configure

Download the latest release using conda: 

```
conda install -c jakevc gopeaks
```

Or download binary asset directly from github: 

```
wget -O gopeaks https://github.com/maxsonBraunLab/gopeaks/releases/download/v0.1.9/gopeaks-linux-amd64
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
    	Minimum width to be considered a peak (default 150)
  -mr int
    	Min reads per coverage bin to be considered (default 15)
  -o string
    	Output prefiex to write peaks and metrics file (default "sample")
  -pval float
    	Pvalue threshold for keeping a peak bin (default 0.05)
  -slide int
    	Slide size for coverage bins (default 50)
  -step int
    	Bin size for coverage bins (default 100)
  -version
    	Print the current gopeaks version
```


## call peaks on a bam file using an IgG control

```
./gopeaks -bam ../chr1.bam -control chr1_igg.bam -cs data/hg38.known.chrom.sizes -o chr1
```

## Output

Two output files are generated each with the output prefix ${prefix}, set to "sample" by default.

    - sample_peaks.bed
    - sample_gopeaks.json

```
head sample_peaks.bed
chr1	9950	10550
chr1	21250	22650
chr1	96050	97050
```

```
cat sample_gopeaks.json
{
	"gopeaks_version": "0.1.9",
	"date": "2021-08-06 11:4:58 AM",
	"elapsed": "1m23.43085221s",
	"prefix": "sample",
	"peak_counts": 4765
}
```
