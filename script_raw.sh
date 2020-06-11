#! /bin/bash


workdir=/home/output/workdir
hicupdir=${workdir}/hicup
preseqdir=${workdir}/preseq
hicup_path=/home/src/HiCUP-master/hicup
preseq=/home/src/preseq
script20thPerc=/home/src/script_20th_perc.R
finalpredictdir=/home/output
script_equaR=/home/src/script_equa.R
script_calcR=/home/src/script_predict.R
script_calcR_no_preseq=/home/src/script_predict_no_preseq.R
figuresdir=${workdir}/figures
tmp=/home/output/tmp

# ./script_raw.sh $species $processors $enzyme $fastq1 $fastq2

processors=$2
species=$1
enzyme=$3
fastq1=$4
fastq2=$5

chromsizes=/home/genomes/${species}/${species}.chrom.sizes
hicuptemplate=/home/templates/hicup_${species}_${enzyme}_template.conf


#################################################################################

## Preparation of the directories

if [ ! -d $workdir ]; then
    mkdir $workdir
fi

cd $workdir

for d in $hicupdir $preseqdir $finalpredictdir $figuresdir $tmp; do
    if [ ! -d $d ]; then
	mkdir $d
    fi
done

#################################################################################

## Mappping and filtering of the reads

cd $hicupdir

# Splitting the fastq files to maximise the cpus use by HiCUP

>&2 echo -e "\tSplitting the read files"
name1=${fastq1##*/}
name2=${fastq2##*/}
if [[ $fastq1 =~ .*.gz ]]; then
    gunzip -c $fastq1 | split -a 5 -l 10000000 - ${name1%.fastq*}_split_
else
    split -a 5 -l 10000000 $fastq1 ${name1%.fastq*}_split_
fi
if [[ $fastq2 =~ .*.gz ]]; then
    gunzip -c $fastq2 | split -a 5 -l 10000000 - ${name2%.fastq*}_split_
else
     split -a 5 -l 10000000 $fastq2 ${name2%.fastq*}_split_
fi

# Running HiCUP

>&2 echo -e "\tRunning HiCUP"
if [ -d hicup_out ]; then
    rm -r hicup_out
fi
mkdir hicup_out
cat $hicuptemplate | sed "s/Threads: 40/Threads: $processors/g"  > hicup.conf
for file in ${name1%.fastq*}_split_*; do
    echo $file >> hicup.conf
    echo ${name2%.fastq*}_split_${file#*_split_} >> hicup.conf
    echo "" >> hicup.conf
done
if [ -f hicup_log.txt ]; then
    rm hicup_log.txt
fi
${hicup_path} --conf hicup.conf >> hicup_log.txt 2>> hicup_log.txt

if (( `ls ${hicupdir}/hicup_out/*.hicup.bam | wc -l` < 1 )); then
    >&2 echo -e "Error: HiCUP execution failed.\n See HiCUP log bellow for troubleshooting:\n"
    >&2 cat hicup_log.txt
    rm -r $tmp
    rm -r $workdir
    exit 1
fi

# Extracting the  mapping and valid reads percentages

echo -e "Raw_pairs\tPerc_map\tPerc_valid\tPerc_cis" > ${figuresdir}/hicup_stats.txt
perc_map=`cat ${hicupdir}/hicup_out/HiCUP_summary_report_*.txt | awk 'BEGIN{n=0;l=0} !/^File/{n=n+$32; l=l+1} END{print n/l}'`
perc_valid=`cat ${hicupdir}/hicup_out/HiCUP_summary_report_*.txt | awk 'BEGIN{n=0;l=0} !/^File/{n=n+$33; l=l+1} END{print n/l}'`
perc_cis=`cat ${hicupdir}/hicup_out/HiCUP_summary_report_*.txt | awk 'BEGIN{n=0;l=0} !/^File/{n=n+$35; l=l+1} END{print 100-(n/l)}'`
echo -e "Total\t"$perc_map"\t"$perc_valid"\t"$perc_cis >> ${figuresdir}/hicup_stats.txt

if (( $(echo "$perc_map<50" | bc -l) )); then
    >&2 echo "Error: too low percentage of reads mapped ("$perc_map"%)."
    rm -r $tmp
    rm -r $workdi
    exit 1
fi

if (( $(echo "$perc_valid<60" | bc -l) )); then
    >&2 echo "Warning: low percentage of valid interactions ("$perc_valid"%)."
fi

if (( $( echo "$perc_cis<70" | bc -l) )); then
    >&2 echo "Warning: low proportion of cis interactions ("$perc_cis"%)."
fi

################################################################################# 

## Unique read pairs prediction for different library size

>&2 echo "Computing duplicates stats"

cd ${preseqdir}

# Starting from the non-filtered mapped read pairs

if [ -f tmp.txt ]; then
    rm tmp.txt
fi
>&2 echo -e "\tprediction of unique read pairs using preseq"
for file in ${hicupdir}/hicup_out/*.pair.bam; do
    samtools view $file | awk '//{if($7=="="){if($4<$8){print $3,$4,$3,$8} else {print $3,$8,$3,$4}} else {if($3<$7){print $3,$4,$7,$8} else print $7,$8,$3,$4}}' OFS='paired' >> tmp.txt
    if [ ! -f tmp.txt ]; then
	>&2 echo "Error: samtools failed"
	rm -r $tmp
	rm -r $workdir
	exit 1
    fi
done

# Sorting the pair by position 

sort -k1,1 -k2,2n -k3,3 -k4,4n --parallel=$processors tmp.txt > tmp_srt.txt
if [ ! -f tmp_srt.txt ]; then
    >&2 echo "Error: sorting failed"
    rm -r $tmp
    rm -r $workdir
    exit 1
fi
rm tmp.txt

# Calculating the number of copies for each read pair
## NB:n/2 because two lines per copy

cat tmp_srt.txt | awk 'BEGIN{line="";n=1} //{if(line!=$0){if(NR>1){print n/2}; line=$0; n=1} else{n=n+1}} END{print n/2}' > tmp_stats.txt
rm tmp_srt.txt

# Running preseq

preseq_out=0
$preseq lc_extrap -o ${figuresdir}/preseq_100M_raw_pairs.txt -V tmp_stats.txt
rm tmp_stats.txt
if [ ! -f ${figuresdir}/preseq_100M_raw_pairs.txt ]; then
    >&2 echo -e "Error: preseq failed.\n\
\t\tUsing mapped interactions only to compute the result.\n\
\t\tThis means that predictions will be done on unique sequenced pairs, not on total sequenced reads."
    preseq_out=1
fi



#################################################################################

## Calculating the coefficiants of the equation to calculate the resolution

>&2 echo "Calculating the resolution equation"

#  Starting from the filtered mapped read pairs, removing manually the duplicates

>&2 echo -e "\tExtracting positions from bam file"
for file in ${hicupdir}/hicup_out/*.filt.bam; do
    samtools view $file | awk '//{if($7=="="){if($4<$8){print $3,$4,$3,$8} else {print $3,$8,$3,$4}} else {if($3<$7){print $3,$4,$7,$8} else print $7,$8,$3,$4}}' OFS='+' >> tmp.txt
    if [ ! -f tmp.txt ]; then
	>&2 echo "Error: samtools failed"
	rm -r $tmp
	rm -r $workdir
	exit 1
    fi
done

>&2 echo -e "\tSorting read pairs"
sort -k1,1 -k2,2n -k3,3 -k4,4n --parallel=$processors tmp.txt > tmp_srt.txt
if [ ! -f tmp_srt.txt ]; then
    >&2 echo "Error: sorting failed"
    rm -r $tmp
    tm -r $workdir
    exit 1
fi
rm tmp.txt
cat tmp_srt.txt | awk 'BEGIN{line="";m=0} //{if(line!=$0){if(NR>1){print line > "tmp_pair1.txt"; m=m+1}; line=$0}} END{print m > "all_count.txt"}'
rm tmp_srt.txt

# Checking if there are enough valid read pairs

n=`cat all_count.txt`
if (( $n<30000000 )); then
    >&2 echo "Error, too low number of mapped reads in the library."
    rm -r $tmp
    rm -t $workdir
    exit 1
fi

# Splitting the read file to maximise cpus use

>&2 echo -e "\tSubsampling read pairs"
split -a 5 -l 10000000 tmp_pair1.txt tmp_pair1_split_

# Generating subsamples tag files for 20M (low) valid read pairs, intermediate (mid) (mean between 20M and the total) and total mapped pairs (all)

mid=`echo "(${n}+20000000)/2" | bc`
s1=`echo "scale=2; 20000000/${n}" | bc`
s2=`echo "scale= 2; ${mid}/${n}" | bc`

for file in tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]; do
    echo "cat $file | awk -v name=${file} -v s1=$s1 -v s2=$s2 'BEGIN{FS=\"+\"; l=0; m=0; srand()} //{print \$1\"\t\"\$2\"\t\"\$2+100 >> (name \"_all.bed\"); print \$3\"\t\"\$4\"\t\"\$4+100 >> (name \"_all.bed\"); if(rand()<s1){print \$1\"\t\"\$2\"\t\"\$2+100 >> (name \"_low.bed\"); print \$3\"\t\"\$4\"\t\"\$4+100 >> (name \"_low.bed\"); l=l+1}; if(rand()<s2){print \$1\"\t\"\$2\"\t\"\$2+100 > (name \"_mid.bed\"); print \$3\"\t\"\$4\"\t\"\$4+100 >> (name \"_mid.bed\"); m=m+1}} END{print l > (name \"_low_count.txt\"); print m > (name \"_mid_count.txt\")}'" >> parallel.sh
done

cat parallel.sh | xargs -P $processors -L 1 -0 bash -c
rm parallel.sh

# Merging the sub-sampled splitted files

for file in tmp_all.bed tmp_low.bed tmp_mid.bed; do
    if [ -f $file ]; then
	rm $file
    fi
done
cat tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]_all.bed >> tmp_all.bed
cat tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]_low.bed >> tmp_low.bed
cat tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]_mid.bed >> tmp_mid.bed
rm tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]
rm tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]_all.bed
rm tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]_low.bed
rm tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]_mid.bed

if [ ! -f tmp_all.bed ] || [ ! -f tmp_low.bed ] || [ ! -f tmp_mid.bed ]; then
    >&2 echo "Error while generating tags files"
    rm -r $tmp
    rm -r $workdir
    exit 1
fi

# Sorting the tag files to use them with bedtools

>&2 echo -e "\tSorting reads"
sort -T $tmp --parallel=$processors -k1,1 -k2,2n --output=tmp_all_srt.bed tmp_all.bed
sort -T $tmp --parallel=$processors -k1,1 -k2,2n --output=tmp_low_srt.bed tmp_low.bed
sort -T $tmp --parallel=$processors -k1,1 -k2,2n --output=tmp_mid_srt.bed tmp_mid.bed
if [ ! -f tmp_all_srt.bed ] || [ ! -f tmp_low_srt.bed ] || [ ! -f tmp_mid_srt.bed ]; then
     >&2 echo "Error: sorting failed"
     rm -r $tmp
     rm -r $workdir
     exit 1
fi

# Extracting the real read count for each subsample
## NB: the subsampling is based on a random function and output an approximate number of reads.

cat tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]_low_count.txt | awk 'BEGIN{n=0} //{n=n+$1} END{print n > "low_count.txt"}'
cat tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]_mid_count.txt | awk 'BEGIN{n=0} //{n=n+$1} END{print n > "mid_count.txt"}'
rm tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]_low_count.txt
rm tmp_pair1_split_[a-z][a-z][a-z][a-z][a-z]_mid_count.txt

# Using bedtools to assess the coverage of the subsamples on 20kb, 50kb and 100kb windows

>&2 echo -e "\tCalculating coverage for 20kb, 50kb and 100kb windows"
sort --parallel=$processors -T $tmp -k1,1 -k2,2n $chromsizes > genome_sorted.chrom.sizes
if [ ! -f genome_sorted.chrom.sizes ]; then
    >&2 echo "Error: sorting failed"
    rm -r $tmp
    rm -r $workdir
    exit 1
fi
echo -e "SizeRaw\tSizeU\t\tWindow\tPerc_20th" > ${figuresdir}/res_data.txt
for i in 20 50 100; do
    j=`echo ${i}*1000 | bc`
    bedtools makewindows -w $j -g genome_sorted.chrom.sizes > genome_windows_${i}kb.bed
    if [ ! -f genome_windows_${i}kb.bed ]; then
	>&2 echo "Error: bedtools failed"
	rm -r $tmp
	rm -r $workdir
	exit 1
    fi
done

for i in 20 50 100; do
    for name in low mid all; do
	file=tmp_${name}_srt.bed
	echo "bedtools intersect -sorted -c -b $file -a genome_windows_${i}kb.bed | awk '{print \$1,\$2,\$3,\$4}' OFS='\t' > ${file%_srt.bed}_${i}kb.bedGraph" >> parallel.sh
    done
done
cat parallel.sh | xargs -P $processors -L 1 -0 bash -c
rm parallel.sh

# Calulating  the 20th percentile for each subsample - window size pair

>&2 echo -e "\tCalculating the 20th percentile of coverage for each subsample / window size"
for i in 20 50 100; do
    for name in low mid all; do
	file=tmp_${name}_srt.bed
	if [ ! -f ${file%_srt.bed}_${i}kb.bedGraph ]; then
            >&2 echo "Error: bedtools failed"
	    rm -r $tmp
	    rm -r $workdir
	    exit 1
	fi
	sizeU=`cat ${name}_count.txt`
	perc=`Rscript $script20thPerc ${file%_srt.bed}_${i}kb.bedGraph`
	echo -e "NA\t"$sizeU"\t"$((i*1000))"\t"$perc >> ${figuresdir}/res_data.txt
    done
done

# Using the 6 20th percentiles calculated to check the linearity and solve the resolution equation

>&2 echo -e "\tExtracting the coefficiants of the resolution equation"
Rscript $script_equaR ${figuresdir}/res_data.txt ${finalpredictdir}/equation.txt


#################################################################################

## Prediction table and graph generation

cd $finalpredictdir

>&2 echo -e "Prediction table generation"
touch ${finalpredictdir}/equation.txt
n=`wc -l equation.txt | cut -d" " -f 1`
if (( $n==5 )); then
    i=0
    while read line; do
	if [[ $i > 0 ]]; then
	    declare "coef_$i"=`echo $line | cut -d" " -f 3`
	fi
	i=$((i+1))
    done < equation.txt
    perc_mapp=`cat ${figuresdir}/hicup_stats.txt | awk '/^Total/{print $2}'`
    perc_valid=`cat ${figuresdir}/hicup_stats.txt | awk '/^Total/{print $3}'`
    perc_cis=`cat ${figuresdir}/hicup_stats.txt | awk '/^Total/{print $4}'`
    s=`cat ${preseqdir}/all_count.txt`
    res=`echo "( 1000 - $coef_4 - ( $coef_3 * $s ) ) / ( ( $coef_1 * $s ) - $coef_2 )" | bc`
    echo "************************************************************"
    echo ""
    echo "The resolution of your libary is "$res" bp."
    echo ""
    if (( $preseq_out==0)); then
	Rscript $script_calcR $coef_1 $coef_2 $coef_3 $coef_4 $perc_mapp $perc_valid $perc_cis ${figuresdir}/preseq_100M_raw_pairs.txt ${finalpredictdir}/prediction.txt
    else
	Rscript $script_calcR_no_preseq $coef_1 $coef_2 $coef_3 $coef_4 $perc_mapp $perc_valid $perc_cis ${finalpredictdir}/prediction.txt
    fi
else
    >&2 echo "error, no prediction computed"
    rm -r $tmp
    rm -r $workdir
    exit 1
fi

################################################################################# 

## Removing temporary files

rm -r $workdir
rm -r $tmp
