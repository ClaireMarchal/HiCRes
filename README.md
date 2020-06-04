# HiCRes

Estimating and predicting HiC library resolution

## Introduction

The purpose of this container is to estimate the resolution of your HiC library and to predict the resolution the same library will reach when sequenced at deeper level. There are two main functionality to use HiCRes: you can either start from your analyzed library by providing a bam file containing the valid read pairs and the genome index you used (much faster option) or you can start from your raw sequencing data by providing the 2 fastq files, the enzyme you used, and the species (much slower, limited to human and mouse for now, and MboI, HindIII or Arima digestions).


## Quick start

### Using docker

#### Starting from fastq files

`docker run --rm -v /Path/To/Your/Fastq/Folder:/home/input:ro -v /Path/To/The/Ouput/Folder:/home/output:rw marchalc/hicres -m raw -t 40 -e MboI -1 sample_R1.fastq -2 sample_R2.fastq -s hg38`

This command will run hicres on your library "sample", using human genome, MboI digestion and using 40 threads. It will look for your fastq files within /Path/To/Your/Fastq/Folder and will stock large temporary files and output the results within /Path/To/The/Ouput/Folder.

Tip: You can start from a subsampled library (100M read pairs) for faster results.

#### Starting from bam of valid interactions

`docker run --rm -v /Path/To/Your/Files/Folder:/home/input:ro -v /Path/To/The/Ouput/Folder:/home/output:rw marchalc/hicres -m bam -t 40 -c your_genome.chrom.sizes -b your_valid_interactions.bam`

This command will run hicres on bam file located in /Path/To/Your/Files/Folder and will look for the chrom.sizes file corresponding to the genome you used within /Path/To/Your/Files/Folder. It will use 40 threads and will stock large temporary files and output the results within /Path/To/The/Ouput/Folder.

### Using singularity

#### Starting from fastq files

`singularity run --bind /Path/To/Your/Fastq/Folder:/home/input --bind /Path/To/The/Ouput/Folder:/home/output docker://marchalc/hicres -m raw -t 40 -e MboI -1 sample_R1.fastq -2 sample_R2.fastq -s hg38`
This command will run hicres on your library "sample", using human genome, MboI digestion and using 40 threads. It will look for your fastq files within /Path/To/Your/Fastq/Folder and will stock large temporary files and output the results within /Path/To/The/Ouput/Folder.

Tip: You can start from a subsampled library (100M read pairs) for faster results.

#### Starting from bam of valid interactions

`singularity run --bind /Path/To/Your/Files/Folder:/home/input --bind /Path/To/The/Ouput/Folder:/home/output docker://marchalc/hicres -m bam -t 40 -c your_genome.chrom.sizes -b your_valid_interactions.bam`

This command will run hicres on bam file located in /Path/To/Your/Files/Folder and will look for the chrom.sizes file corresponding to the genome you used within /Path/To/Your/Files/Folder. It will use 40 threads and will stock large temporary files and output the results within /Path/To/The/Ouput/Folder.

## Usage

`hicres -m [raw,bam,bam_fast] [options]`

Arguments:

-m, --method [raw,bam,bam_fast]     The method raw starts from fastq files and output the resolution of the library versus the number of read pairs sequenced. The fast version of the bam option (bam_fast) is in beta mode. See the Fast modes section bellow for more information.

-s, --species [hg38,mm10]     The species (genome name) from which the sample comes from. Either hg38 for human or mm10 for mouse. This is required for method raw and is ignored for method bam.

-e, --enzyme [HindIII,MboI,Arima]     The restriction digestion method, either HindIII for HindIII digestion, MboI for MboI digestion or Arima, for the Arima kit. This is required for method raw and is ignored for method bam.

-c, --chromsize <path to file>     The path to index of the genome used to anaylze the HiC, for method bam and bam_fast only.

-1, --fastq_1 <path to file>     The path to the first end of the sequenced reads, in fastq or fastq.gz format. For method raw only.

-2, --fastq_2 <path to file>     The path to the second end of the sequenced reads, in fastq or fastq.gz format. For method raw only.

-b, --bam <path to file>     The path to bam file of the valid read pairs, for method bam and bam_fast only.

-t, --threads <int>     The number of threads to use. Default is 1.

## Fast mode

On normal mode, HiCRes tries to keep the paired-end imformation from the valid interactions and subsamples the datasets keeping the read pairs together. Fast mode ignores the pairing of the valid interactions, resulting in a faster execution. Neverless, this a beta mode, since all the tests performed to check the accuracy of this tool have been done using normal mode. You should be cautious using this mode.

To use this mode, use method "bam_fast" instead of "bam."

## Using individual sub-programs (ongoing section)

This tool is composed of several bash and R scripts. Some may be of interest for users who do not want to run the docker.

### script_raw.sh


### script_bam.sh and script_bam_fast.sh


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

- "Error: too low  percentage of reads mapped"

This error is issued if less than 50%  of the read pairs are mapped. In this case, it is likely that there is a problem. Check if you selected the correct species. If you selected the correct species, several problems are possible:

     - The fastq have not only read pairs: if you pre-processed the fastq to remove adpaters for example, make sure the final fastq files contains read pairs only, both file in  the same order. These files will be splitted to use use multiple cpus. Thus, if the read pairs are not in the  same order, they will be proccessed separetly and will not be paired after the mapping.

     - There is a contamination of adapters in the fastq. Adapters are not removed before  processing. Thus, if you have a contaminatioin, it is possible that you loose too many reads. In this case clip the adpters before usinig HiCRes, and make sure to keep only read pairs, in the same order ni both fastq files.

     - There is a contamination with another species. If your samples are contaminated (by mixing bymistake to samples with the same index for example), your samples are not usable. You can seperate the species by mapping on chimeric genome your reads and selecting the one mapping on your genome of interest, but this solution is for troubleshooting only, and not recommanded to perform further analysis as you  have a risk to introduce a bias in your data.

- "Warning: low percentage of valid interactions" and "Warning: low proportioin of cis interactions."

In both cases, this indicates a poor quality for the library. HiCRes will still compute the resolution, but it is recommended to try to improve the experimental conditions before sequencing deeper. In the case of low cis interactions proportion, iit is especially important to work with the resolution computed from cis interactions only.

- "Error: preseq failed"

This should be preceded by preseq error message. Use this message for refering to preseq documentation if you want to understand why preseq failed. When preseq fails, the library yields in function of the sequencing depth cannot be predicted. In this case, HiCRes will issue a warning and predict the resolution in function of the unique valid read pairs sequenced, instead of the toatl number of the  total number of read pairs sequenced (including duplicates).

- "Error: too low number of mapped reads in the library."

HiCRes has been succesfuly tested to work with datasets having as low as 30M valid read pairs. If your dataset is smaller, HiCRes will not process it as we don't know how accurate this will be.

- "Error while generating the tag files." and "Error: bedtools failed."

This error is likely caused by  bedtools not workinig properly. If you get this error while running HiCRes on processed data, make sure the genome index (.chrom.sizes file) you gave to HiCRes corresponds to genome used to process the data.

- "Error: sorting failed."

There are several steps that require sorting of big files, using the function sort. It has been configured to use a temporary folder created in the output directory. This error could be caused by a lack of memory in this directory.

- "Error: no prediction computed."

This error occurs either if the equation.txt file has not been generated (unkown error), or if the equation has not been predicted (most likely). In this second case, an equation.txt file should be present in the output directory, containing the error details instead of the final equation. The error is the absence of linearity of the distribution of the mapped reads in function of the read number or in function of the windows size (specified in the equation.txt file). This will happen if you mixed several libraries together. This tool cannot process several libraries. The linearity of the distribution is the main condition to extract the equation to predict the resolution of the library in function of the depth. Thus, if the distribution is not linear, it unfortunately impossible to make any prediction. If you are analysing a single library, trying to start with a deeper library may improve that.

## Benchmarking (ongoing section)

Below are the times needed to process HiC data using the raw method (starting from fastq files). This times have been measured using singularity on an HPC server, allocating 40 cpus. 

| Dataset | Size (read pairs) | Species | Enzyme | Time | Resolution |
| :-----: | :---------------: | :-----: | :----: | :--: | :--------: |
| SRR1658692 | 274M | Human | HindIII | 5h12m<sup>1</sup> | |
| SRR1658573<sup>2</sup> | 161M | Human | MboI | 4h55m | |
| SRR443883 SRR443884 SRR443885 | 465M | Mouse | HindIII | 4h36m | |
| SRR9906313<sup>2</sup> | 270M | Mouse | MboI | 6h53m | |
| Unpublished | 195M | Mouse | Arima | 5h14m | |

1. Including the time to pull the docker from dockerhub.

2. Prediction for the library yields (duplicates) failed. In this case, HiCRes gives a warning and generates the predictions on the uniquely sequenced read pairs.


Below are the times needed to process HiC data	using the bam and bam_fast method (starting from bam files containing valid interactions). This times have been measured using singularity on an HPC server, allocating 40 cpus.

| Datasets | Size (valid interactions) | Species | Time "bam" | Resolution "bam" | Time "bam_fast" | Resolution "bam_fast" |
| :------: | :-----------------------: | :-----: | :--------: | :--------------: | :-------------: | :-------------------: |
| SRR1658692 | 165M | Human | 44m<sup>1</sup> | | 22m | |
| SRR1658573 | 105M | Human | 25m | | 15m | |
| SRR443883 SRR443884 SRR443885 | 145M | Mouse | 29m | | 19m | |
| SRR9906313 | 183M | Mouse | error | | 26m | |
| Unpublished | 120M | Mouse | 29m | | 16m | |

1. Including the time to pull the docker from dockerhub.

## References

Please cite the paper associated to this tool if you use it:

The citation will be updated here when the paper will be published.

This tool uses also:

- samtools: Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. 
- bowtie2: Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357-359 (2012).
- HiCUP: HiCUP: pipeline for mapping and processing Hi-C data. F1000Res 4, 1310 (2015).
- preseq: Daley, T. & Smith, A.D. Predicting the molecular complexity of sequencing libraries. Nat Methods 10, 325-327 (2013).
- R: R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.