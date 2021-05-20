#! /bin/bash

processors=1
species=""
enzyme=""
fastq1=""
fastq2=""
method=""
bam=""
chromsize=""

## Checking if there are arguments
if (( $# % 2 != 0 && $# != 1)); then
    echo -e "Error: bad arguments\nPlease use -h or --help to print the help screen." >&2
    exit 1
fi

## Printing welcome message
echo -e "Welcome to HiCRes!\n\
Please use -h or --help to print the help screen.\n\
If you use this tool, please refer to the GitHub page to cite it: http://github.com/ClaireMarchal/HiCRes\n\
\n\
*********************************\n\
" >&2

## Reading the arguments
while (( "$#" )); do
  case "$1" in
    -h|--help)
    echo -e  "\nThe purpose of this container is to estimate the resolution of your HiC library and to predict the resolution the same library will reach when sequenced at deeper level.\n\
There are two main functionality to use HiCRes: you can either start from your analyzed library by providing a bam file containing the valid read pairs and the genome index you used (much faster option) or you can start from your raw sequencing data by providing the 2 fastq files, the enzyme you used, and the species (much slower, limited to human and mouse for now, and MboI, HindIII or Arima digestions).\n\
Usage:\n\
hicres -m [raw,bam] [options]\n\
\n\
Arguments:\n\
-m, --method [raw,bam,bam_fast]\t The method raw starts from fastq files and output the resolution of the library versus the number of read pairs sequenced. The fast version of the bam option (bam_fast) is in beta mode. See the GitHub page http://github.com/ClaireMarchal/HiCRes for more information. (required)\n\
-s, --species [hg38,mm10,ce10,dm3,TAIR10]\t The species (genome name) from which the sample comes from. Either hg38 for human, mm10 for mouse, ce10 for C. elegans, TAIR10 for A. thaliana or dm3 for D. melanogaster. This is required for method raw and is ignored for method bam.\n\
-e, --enzyme [HindIII,MboI,Arima]\t The restriction digestion method, either HindIII for HindIII digestion, MboI for MboI digestion or Arima, for the Arima kit. This is required for method raw and is ignored for method bam.\n\
-c, --chromsize <path to file>\t The path to index of the genome used to anaylze the HiC, for method bam only.\n\
-1, --fastq_1 <path to file>\t The path to the first end of the sequenced reads, in fastq or fastq.gz format. For method raw only.\n\
-2, --fastq_2 <path to file>\t The path to the second end of the sequenced reads, in fastq or fastq.gz format. For method raw only.\n\
-b, --bam <path to file>\t The path to bam file of the valid read pairs, for method bam only.\n\
-t, --threads <int>\t The number of threads to use. Default is 1.\n\
\n\
Files access help:\n\
Permission to read and write can be difficult to set up with containers. If have errors involving reading / writing, please look at our GitHub (http://github.com/ClaireMarchal/HiCRes) or consult singulariy or docker online helps." >&2
    exit 0
    ;;
    -m|--method)
      method=$2
      shift 2
      ;;
    -s|--species)
      species=$2
      shift 2
      ;;
   -t|--threads)
     processors=$2
     shift 2
     ;;
    -e|--enzyme)
     enzyme=$2
     shift 2
     ;;
   -1|--fastq_1)
     fastq1=$2
     shift 2
     ;;
   -2|--fastq_2)
     fastq2=$2
     shift 2
     ;;
    -b|--bam)
     bam=$2
     shift 2
     ;;
    -c|--chromsize)
     chromsize=$2
     shift 2
     ;;
    --)
      shift
      break
      ;;
    -*|--*=)
      echo "Error: Unsupported flag \"$1\"." >&2
      exit 1
      ;;
    *)
      echo "Error: Unknown option \"$1\"." >&2
      exit 1
      ;;
  esac
done

## Checking if necessary options are here
if [[ $method != "raw" ]] && [[ $method != "bam" ]] && [[ $method != "bam_fast" ]]; then
    echo "Error: method should be either raw, bam or bam_fast (beta version)." >&2
    exit 1
fi

if [[ $species != "hg38" ]] && [[ $species != "mm10" ]] && [[ $species != "ce10" ]] && [[ $species != "dm3" ]] && [[ $species != "TAIR10" ]] && [[ $method == "raw" ]]; then
    echo "Error: species should be either hg38 (for human), mm10 (for mouse), ce10 (for C. elegans), TAIR10 (for A. thaliana) or dm3 (for D. melanogaster)." >&2
    exit 1
fi

if [[ $enzyme != "MboI" ]] && [[ $enzyme != "HindIII" ]]  && [[ $enzyme != "Arima" ]] && [[ $method == "raw" ]]; then
    echo "Error: enzyme should be MboI, HindIII or Arima." >&2
    exit 1
fi

if [[ -z $fastq1 || -z $fastq2 ]] && [[ $method == "raw" ]]; then
    echo "Error: fastq files are not specified." >&2
    exit 1
fi

if  [[ ! -r /home/input/$fastq1 || ! -r /home/input/$fastq2 ]]  && [[ $method == "raw" ]]; then
    echo "Error: fastq files cannot be read." >&2
    exit 1
fi

if [[ -z $bam ]] && [[ $method == "bam" || $method == "bam_fast" ]]; then
    echo "Error: bam file is not specified." >&2
    exit 1
fi

if [[ ! -r /home/input/$bam ]] && [[ $method == "bam" || $method == "bam_fast" ]]; then
    echo "Error: cannot be read." >&2
    exit 1
fi

if [[ $chromsize == "" ]] && [[ $method == "bam" || $method == "bam_fast" ]]; then
    echo "Error: please specify a genome index (.chom.sizes file)." >&2
    exit 1
fi

## Setting up some environmental variables that R needs to avoid R startup warnings
export LC_TIME="C"
export LC_MESSAGES="C"
export LC_MONETARY="C"
export LC_PAPER="C"
export LC_MEASUREMENT="C"
export LC_COLLATE="C"
export LC_CTYPE="C"

## Running the script corresponding to the method choosen by the user
if [[ $method == "raw" ]]; then
    echo "Here are your options for analysing raw data:"
    echo -e "\tspecies: "$species
    echo -e "\tthreads: "$processors
    echo -e "\tenzyme: "$enzyme
    echo -e "\tfastq 1: "$fastq1
    echo -e "\tfastq 2: "$fastq2
    echo -e "Running the script.\n\
Please wait, depending on the number of threads allocated and the size of the library, it may take a long time."
    /home/src/script_raw.sh $species $processors $enzyme /home/input/$fastq1 /home/input/$fastq2
    exit 0
fi

if [[ $method == "bam" ]]; then
    echo "Here are your options for analysing processed data:"
    echo -e "\tthreads: "$processors
    echo -e "\tbam: "$bam
    echo -e "\tchromsize: "$chromsize
    echo -e "Running the script.\n\
Please wait, depending on the number of threads allocated and the size of the library, it may take a long time."
    /home/src/script_bam.sh $processors /home/input/$bam /home/input/$chromsize
    exit 0
fi

if [[ $method == "bam_fast" ]]; then
    echo "Here are your options for analysing processed data in fast mode:"
    echo -e "\tthreads: "$processors
    echo -e "\tbam: "$bam
    echo -e "\tchromsize: "$chromsize
    echo -e "Running the script.\n\
Please wait, depending on the number of threads allocated and the size of the library, it may take a long time.\n\
Please note that the fast mode is a beta version, be carreful when using the result. For more information, please see the GitHub page."
    /home/src/script_bam_fast.sh $processors /home/input/$bam /home/input/$chromsize
    exit 0
fi

