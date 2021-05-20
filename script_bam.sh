#! /bin/bash


workdir=/home/output/workdir
script20thPerc=/home/src/script_20th_perc.R
finalpredictdir=/home/output
script_equaR=/home/src/script_equa.R
script_calcR=/home/src/script_predict_bam.R
figuresdir=${workdir}/figures
tmp=/home/output/tmp

# ./script_bam.sh $species $processors $bam $chromsize

processors=$1
bam=$2
chromsizes=$3

if [ ! -d $workdir ]; then
    mkdir $workdir
fi

cd $workdir

for d in $finalpredictdir $figuresdir $tmp; do
    if [ ! -d $d ]; then
	mkdir $d
    fi
done


#################################################################################

## Calculating the coefficiants of the equation to calculate the resolution

>&2 echo "Calculating the resolution equation"

#  Starting from the unique filtered mapped read pairs

>&2 echo -e "\tExtracting positions from bam file"
samtools view $bam | awk '//{if($7=="="){if($4<$8){print $3,$4,$3,$8} else {print $3,$8,$3,$4}} else {if($3<$7){print $3,$4,$7,$8} else print $7,$8,$3,$4}}' OFS='+' >> tmp.txt
if [ ! -f tmp.txt ]; then
    >&2 echo "Error: samtools failed"
    rm -r $workdir
    rm -r $tmp
    exit 1
fi


>&2 echo -e "\tSorting read pairs"
sort -k1,1 -k2,2n -k3,3 -k4,4n --parallel=$processors tmp.txt > tmp_srt.txt
if [ ! -f tmp_srt.txt ]; then
    >&2 echo "Error: sorting failed"
    rm -r $workdir
    rm -r $tmp
    exit 1
fi
rm tmp.txt
cat tmp_srt.txt | awk 'BEGIN{line="";m=0} //{if(line!=$0){if(NR>1){print line > "tmp_pair1.txt"; m=m+1}; line=$0}} END{print m > "all_count.txt"}'
rm tmp_srt.txt

# Checking if there are enough valid read pairs

n=`cat all_count.txt`
if (( $n<30000000 )); then
    >&2 echo "Error, too low number of mapped reads in the library."
    rm -r $workdir
    rm -r $tmp
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
    rm -r $workdir
    rm -r $tmp
    exit 1
fi

# Sorting the tag files to use them with bedtools

>&2 echo -e "\tSorting reads"
sort -T $tmp --parallel=$processors -k1,1 -k2,2n --output=tmp_all_srt.bed tmp_all.bed
sort -T $tmp --parallel=$processors -k1,1 -k2,2n --output=tmp_low_srt.bed tmp_low.bed
sort -T $tmp --parallel=$processors -k1,1 -k2,2n --output=tmp_mid_srt.bed tmp_mid.bed
if [ ! -f tmp_all_srt.bed ] || [ ! -f tmp_low_srt.bed ] || [ ! -f tmp_mid_srt.bed ]; then
     >&2 echo "Error: sorting failed"
     rm -r $workdir
     rm -r $tmp
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
    rm -r $workdir
    rm -r $tmp
    exit 1
fi
echo -e "SizeRaw\tSizeU\t\tWindow\tPerc_20th" > ${figuresdir}/res_data.txt
for i in 20 50 100; do
    j=`echo ${i}*1000 | bc`
    bedtools makewindows -w $j -g genome_sorted.chrom.sizes > genome_windows_${i}kb.bed
    if [ ! -f genome_windows_${i}kb.bed ]; then
	>&2 echo "Error: bedtools failed"
	rm -r $workdir
	rm -r $tmp
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
	    rm -r $workdir
	    rm -r $tmp
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

## Prediction table generation

cd $finalpredictdir

>&2 echo -e "Prediction table generation"
touch  ${finalpredictdir}/equation.txt
n=`wc -l equation.txt | cut -d" " -f 1`
if (( $n==5 )); then
    i=0
    while read line; do
	if [[ $i > 0 ]]; then
	        declare "coef_$i"=`echo $line | cut -d" " -f 3`
		fi
	i=$((i+1))
    done < equation.txt
    s=`cat ${workdir}/all_count.txt`
    res=`echo "( 1000 - $coef_4 - ( $coef_3 * $s ) ) / ( ( $coef_1 * $s ) - $coef_2 )" | bc`
    echo "************************************************************"
    echo ""
    echo "The resolution of your libary is "$res" bp (using all interactions)."
    echo ""
    Rscript $script_calcR $coef_1 $coef_2 $coef_3 $coef_4 ${finalpredictdir}/prediction.txt
else
    >&2 echo "Error, no prediction computed."
    rm -r $workdir
    rm -r $tmp
    exit 1
fi

################################################################################# 

## Removing temporary files

rm -r $workdir
rm -r $tmp
