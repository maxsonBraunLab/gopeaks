[![Anaconda-Server Badge](https://anaconda.org/bioconda/gopeaks/badges/version.svg)](https://anaconda.org/bioconda/gopeaks)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gopeaks/badges/downloads.svg)](https://anaconda.org/bioconda/gopeaks)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gopeaks/badges/license.svg)](https://anaconda.org/bioconda/gopeaks)
[![Go](https://github.com/maxsonBraunLab/gopeaks/actions/workflows/go.yml/badge.svg?branch=main)](https://github.com/maxsonBraunLab/gopeaks/actions/workflows/go.yml)
![Maintainer](https://img.shields.io/badge/maintainer-gartician-blue)

# GoPeaks

GoPeaks is a peak caller designed for CUT&TAG/CUT&RUN sequencing data. GoPeaks by default works best with narrow peaks such as H3K4me3 and transcription factors. However, broad epigenetic marks like H3K27Ac/H3K4me1 require different the step, slide, and minwidth parameters. We encourage users to explore the parameters of GoPeaks to analyze their data.

# Configure

Download the latest release using conda: 

```
conda install -c bioconda -c conda-forge gopeaks
```

Or download binary asset directly from github: 

```
wget -O gopeaks https://github.com/maxsonBraunLab/gopeaks/releases/download/v1.0.0/gopeaks-linux-amd64
chmod +x gopeaks
```

# Example Usage

```
usage: GoPeaks [-h|--help] [-b|--bam "<value>"] [-c|--control "<value>"]
               [-s|--chromsize "<value>"] [-m|--mdist <integer>] [-r|--minreads
               <integer>] [-p|--pval <float>] [-t|--step <integer>] [-l|--slide
               <integer>] [-w|--minwidth <integer>] [-o|--prefix "<value>"]
               [-v|--version] [--broad] [--verbose]

               GoPeaks is a peak caller designed for CUT&TAG/CUT&RUN sequencing
               data. GoPeaks by default works best with narrow peaks such as
               H3K4me3 and transcription factors. GoPeaks can be used with the
               "--broad" flag to call broad peaks like H3K27Ac/H3K4me1. We
               encourage users to explore the parameters of GoPeaks to analyze
               their data.

Arguments:

  -h  --help       Print help information
  -b  --bam        Input BAM file (must be paired-end reads)
  -c  --control    Input BAM file with control signal to be normalized (e.g.
                   IgG, Input)
  -s  --chromsize  Chromosome sizes for the genome if not found in the bam
                   header
  -m  --mdist      Merge peaks within <mdist> base pairs. Default: 1000
  -r  --minreads   Test genome bins with at least <minreads> read pairs..
                   Default: 15
  -p  --pval       Define significance threshold <pval> with multiple
                   hypothesis correction via Benjamini-Hochberg. Default: 0.05
  -t  --step       Bin size for coverage bins. Default: 100
  -l  --slide      Slide size for coverage bins. Default: 50
  -w  --minwidth   Minimum width (bp) of a peak. Default: 150
  -o  --prefix     Output prefix to write peaks and metrics file. Default:
                   sample
  -v  --version    Print the current GoPeaks version
      --broad      Run GoPeaks on broad marks (--step 5000 & --slide 1000)
      --verbose    Run GoPeaks in verbose mode.

```

## Call peaks on a bam file using an IgG control

```
gopeaks -b <sample>.bam -c <control>.bam -o path/to/gopeaks/<sample>
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
        "gopeaks_version": "1.0.0",
        "date": "2022-05-10 3:19:14 PM",
        "elapsed": "5m20.184082356s",
        "prefix": "K562_1_H3K4me3",
        "command": "gopeaks -b K562_1_H3K4me3.bam -c K562_1_IgG.bam",
        "peak_counts": 10329
}
```

## Recommended parameters

| Sequencing Modality                      | Recommended Parameters       |
| ---------------------------------------- | ---------------------------- |
| CUT&TAG or CUT&RUN narrow peaks          | Default parameters           |
| CUT&TAG or CUT&RUN transcription factors | Default parameters           |
| CUT&TAG or CUT&RUN broad peaks           | `--broad` and `--mdist 3000` |
| ATAC-Seq                                 | Default parameters           |

