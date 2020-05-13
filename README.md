# HiCRes

Estimating and predicting HiC library resolution

## Introduction

The purpose of this container is to estimate the resolution of your HiC library and to predict the resolution the same library will reach when sequenced at deeper level. There are two main functionality to use HiCRes: you can either start from your analyzed library by providing a bam file containing the valid read pairs and the genome index you used (much faster option) or you can start from your raw sequencing data by providing the 2 fastq files, the enzyme you used, and the species (much slower, limited to human and mouse for now, and MboI, HindIII or Arima digestions).

## Quick start

### Using docker

To start from fastq files, use

`docker run --rm -v /Path/To/Your/Fastq/Folder:/home/input:ro -v /Path/To/The/Ouput/Folder:/home/output:rw marchalc/hicres -m raw -t 40 -e MboI -1 sample_R1.fastq -2 sample_R2.fastq -s hg38`

This command will run hicres on your library "sample", using human genome, MboI digestion and using 40 threads. It will look for your fastq files within /Path/To/Your/Fastq/Folder and will stock large temporary files and output the results within /Path/To/The/Ouput/Folder.

To start from bam of valid interactions, use

`docker run --rm -v /Path/To/Your/Files/Folder:/home/input:ro -v /Path/To/The/Ouput/Folder:/home/output:rw marchalc/hicres -m bam -t 40 -c your_genome.chrom.sizes -b your_valid_interactions.bam`

This command will run hicres on bam file located in /Path/To/Your/Files/Folder and will look for the chrom.sizes file corresponding to the genome you used within /Path/To/Your/Files/Folder. It will use 40 threads and will stock large temporary files and output the results within /Path/To/The/Ouput/Folder.

### Using singularity

Will come soon.

## Usage

`hicres -m [raw,bam] [options]`

Arguments:

-m, --method [raw,bam,raw_fast,bam_fast]     The method raw starts from fastq files and output the resolution of the library versus the number of read pairs sequenced. The fast versions of each of these options (raw_fast and bam_fast) are in beta mode. See the Fast modes section bellow for more information.

-s, --species [hg38,mm10]     The species (genome name) from which the sample comes from. Either hg38 for human or mm10 for mouse. This is required for method raw and is ignored for method bam.

-e, --enzyme [HindIII,MboI,Arima]     The restriction digestion method, either HindIII for HindIII digestion, MboI for MboI digestion or Arima, for the Arima kit. This is required for method raw and is ignored for method bam.

-c, --chromsize <path to file>     The path to index of the genome used to anaylze the HiC, for method bam only.

-1, --fastq_1 <path to file>     The path to the first end of the sequenced reads, in fastq or fastq.gz format. For method raw only.

-2, --fastq_2 <path to file>     The path to the second end of the sequenced reads, in fastq or fastq.gz format. For method raw only.

-b, --bam <path to file>     The path to bam file of the valid read pairs, for method bam only.

-t, --threads <int>     The number of threads to use. Default is 1.

## Fast modes

On normal mode, HiCRes tries to keep the paired-end imformation from the valid interactions and subsamples the datasets keeping the read pairs together. Fast modes ignore the pairing of the valid interactions, resulting in a faster execution. Neverless, this a beta mode, since all the tests performed to check the accuracy of this tool have been done using normal mode. You should be cautious using these modes.
To use these modes, use methods "raw_fast" instead of "raw" and "bam_fast" instead of "bam."

## Using individual sub-programs

This tool is composed of several bash and R scripts. Among these last, two scripts within the folder scriptsR may be of interest for users who do not want to run the docker.

### script_20th_perc.R

usage: `Rscript path/to/script_20th_perc.R coverage.bed`

This script takes in input a coverage.bed file containiing the read number per window on the whole genome and output the 20th percentile. The input is a 4 column tab delimited file, with no header, for which the first column is the chr, the second, the start of the window, the third, the end of the window and the fourth, the number of reads mapping within the window.
This is script is used within the docker to output the 20th percentiles of the valid read coverage from several window and library sizes.

### script_equa.R

usage: `Rscript path/to/script_equa.R input_20th_percentiles.txt output_equation.txt`

This script takes in input a tab delimited file containing 9 data points (see bellow), and output a text file containing the equation to calculate the resolution of the HiC library for any valid reads number. This script check for the linearity between the 20th percentile of the read pairs distribution and the window size and for the linearity between the 20th percentile of the read pairs distribution and the sample size (number of valid read pairs). In case of non-linearity, the output file will contain only an error message. 

The input is the 20th percentile value for three sample sizes (20M and 2 values > 20M read pairs) and three window sizes (20kb, 50kb, 100kb). It  must have these columns, with headers:

- column 1, "SizeRaw": This column is not used by this program
- column 2, "SizeU": The number of read pairs. It must have one and only one sample containing 20000000 reads or less (try to get close to this number). And at least two samples with more read pairs.
- column 3, "Window": The window size used to calculate the read distribution in bp. It must include 20000 for each sample size and at least 2 other size larger than 20000bp (e.g. 50000 and 100000)
- column 4, "Perc_20th": The 20th percentile of the distribution of these reads 


example of input file:

| SizeRaw | SizeU | Window | Perc_20th |
| :-----: | :---: | :----: | :-------: |
| NA | 19387136 | 20000 | 132 |
| NA | 42006016 | 20000 | 290 |
| NA | 64622583 | 20000 | 447 |
| NA | 19387136 | 50000 | 346 |
| NA | 42006016 | 50000 | 753 |
| NA | 64622583 | 50000 | 1159 |
| NA | 19387136 | 100000 | 700 |
| NA | 42006016 | 100000 | 1522 |
| NA | 64622583 | 100000 | 2342 |

## Troubleshooting



## Benchmarking

Will come soon.

## References

Please cite the paper associated to this tool if you use it:

The citation will be updated here when the paper will be published.

This tool uses also:

- samtools: Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. 
- bowtie2: Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357-359 (2012).
- HiCUP: HiCUP: pipeline for mapping and processing Hi-C data. F1000Res 4, 1310 (2015).
- preseq: Daley, T. & Smith, A.D. Predicting the molecular complexity of sequencing libraries. Nat Methods 10, 325-327 (2013).
- R: R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.
- parallel: https://doi.org/10.5281/zenodo.1146014