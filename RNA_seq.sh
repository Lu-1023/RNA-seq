echo "usage:"
echo "bash rnaseq.sh quality_cutoff(20/30) species(hg38. ce11) feature(exon, transcript) name(gene_id, transcript_id)"
species=$2  
feature=$3
name=$4

if [ ! -n "$3" ];then
	echo “no enough parameters
	exit 8
fi

source ~/.bashrc
source activate qc
if [[ $1 == 20 ]] || [[ $1 == 30 ]]
then
    echo $1
else
    echo "without quality_cutoff(20/30)"
    exit 8
fi
# need all fastq file in fq files
mkdir qc
mkdir fastqc
cd fq 
for i in *.fq.gz
do
    ii=`echo $i | sed "s/.*_L0._//g"`
    # remove the sequencing mark (FXXXX_L01_sample.fq.gz)
    mv $i $ii
done


cd ..
##trimadapter
mkdir trimmed
for line in $x
do
    trim_galore --paired -q $1 fq/${line}_1.fq.gz  fq/${line}_2.fq.gz  -o trimmed &
done
cp *.txt fastqc
cd ..
# remove the rRNA reads
mkdir de_rRNA
touch rRNA_report
for line in $x
do
    echo 'rRNA_sample:'$line >> rRNA_report
    bowtie -p 16 --un de_rRNA/${line}.fq -x ~/reference/${species}/ribosome/rDNA -1 trimmed/${line}_1.fq.gz -2 trimmed/${line}_2.fq.gz -S useless.sam 2>>rRNA_report
done
rm *.sam
# get the quality of fastq

for line in $x
do
    fastqc de_rRNA/${line}_1.fq de_rRNA/${line}_2.fq -o fastqc &
done
wait

mkdir aligned
for line in $x
do
    STAR --genomeDir ~/reference/hg38/star_index/ --runThreadN 20 \
    --readFilesIn de_rRNA/${line}_1.fq de_rRNA/${line}_2.fq  --outFileNamePrefix aligned/${line} --outSAMtype SAM \
    --outSAMunmapped Within --outSAMattributes Standard
    mv aligned/${line}Aligned.out.sam aligned/${line}.sam
    mv aligned/${line}Log* fastqc
    samtools sort -@ 16 aligned/${line}.sam -o aligned/${line}.bam
    rm aligned/${line}.sam
    samtools index -@ 16 aligned/${line}.bam
done

mkdir ../bw
# generate none normalized bw file
for i in *.bam
do
    bamCoverage -b $i -o ../bw/${i%.*}_fw.bw --normalizeUsing None --binSize 10 --filterRNAstrand forward -p 16
    bamCoverage -b $i -o ../bw/${i%.*}_rev.bw --normalizeUsing None --binSize 10 --filterRNAstrand reverse -p 16
    # calculate the total reads wiht samtools
done

    # generate featureCounts.txt for DESeq2 analysis
    featureCounts -p -T 16 -s 2 -a ~/reference/$species/ref_nmr.gtf -t $feature -g $name -o featureCounts.txt $bam
    # analyze with DESeq2 in R

