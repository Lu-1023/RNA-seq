**<font color="grey"><font size=10>UPSTREAM ANALYSIS of RNA-Seq </font></font>**
<font size=5><font color="grey"><p align="right">2020.10.27</p></font></font>
# <font color="steelblue">Pipe for RNA-seq(polyA enrichment,mRNA) </font>

[TOC]

***    
##  <font size=4>1   remove adapter(trim_galore)</font>
 Trim the adapter of the reads.
```shell
trim_galore --paired ${line}_1.fastq.gz ${line}_2.fastq.gz -o trimmed_fastq  #pair end
trim_galore ${line}.fastq.gz -o trimmed_fastq  #single end
```

##  <font size=4>2   quality control (fastQC)</font>
Get the basic and quality information of the library.

```shell
fastqc -o ./fastqc -t 16 ${line}/${line}_1.fq.gz
```

##  <font size=4>3   mapping reads(RNA-seq mappers)</font>
Using the RNA-seq mappers , such as <kbd>hisat2</kbd> , <kbd>bowtie</kbd> , <kbd>bowtie2</kbd> , <kbd>STAR</kbd> or another , mapping the reads against the genome reference and identifying their genomic positions.
```shell
hisat2 -p 16 --dta --rna-strandness RF -x /media/hp/disk1/song/Genomes/${species}/Sequence/WholeGenomeFasta/hisat2/genome -1 ${line}/${line}_1_val_1.fq.gz -2 ${line}/${line}_2_val_2.fq.gz -S align/${line}.sam 2>> align/mapping_report.txt
```
```shell
bowtie -p 16 -S /media/hp/disk1/song/Genomes/${species}/${species}_ref/bowtie/longest ${line}_timmed.fastq align/${line}.sam 2>> align/mapping_report.txt
```
```shell
bowtie2 -S ${input}.sam -p 16 -5 5 -3 10 -x /media/hp/disk1/song/Genomes/NC10/Sequences/WholeGenomeFasta/bowtie2/NC10 -1 ${input}_1_val_1.fq.gz -2 ${input}_2_val_2.fq.gz 2>> align/mapping_report.txt
```
```shell
mkdir /media/hp/disk1/song/Genomes/${species}/Sequence/STAR
STAR --runMode genomeGenerate --runThreadN 16 --genomeFasta /media/hp/disk1/song/Genomes/${species}/Sequence/WholeGenomeFasta/genome.fa --genomeDir /media/hp/disk1/song/Genomes/${species}/Sequence/STAR --sjdbGTFfile /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf --sjdOverhang 149 
# assign length of reads
    --sjdOverhang represents the length of reads
STAR --runThreadN 20 --genomeDir /media/hp/disk1/song/Genomes/${species}/Sequence/STAR --readFilesCommand zcat --readFilesIn ${line}_1_val_1.fq ${line}_2_val_2.fq --sjdbOverhang 149 --outSAMtype BAM SortedBycoordinate   --outFileNamePrefix align/${line} #pair end
STAR --runThreadN 20 --genomeDir /media/hp/disk1/song/Genomes/${species}/Sequence/STAR --readFilesIn ${line}_trimmed.fq --sjdbOverhang 149 --outSAMtype BAM SortedBycoordinate  --outFileNamePrefix align/${line} #single end
#if your input file is compressed (*.fq.gz) add the option
    --readFilesCommand zcat
    --readFilesCommand gzip -c 
#assign the format of  output file, such as BAM or SAM
    --outSAMtype BAM Unsorted #sort by name and can link to HTSeq
    --outSAMtype BAM SortedByCoordinate
    --outSAMtype BAM Unsorted SortedByCoordinate #output two formats mentioned above
mv align/${line}Aligned.sortedByCoord.out.bam align/${line}.bam
mv align/${line}Log.fianl.out fastqc
```
##  <font size=4>4   integrated information(multiqc)</font>
Use multiqc to merge the quality message(fastqc) , trim adapter reports(trim_galore) and mapping results(STAR) of all samples into one report file.
```shell
multiqc fastqc -o ./fastqc
```

##  <font size=4>5   sorting alignment and converting(samtools)</font>
Sorting the alignment by the genomic positions or names .
```shell
samtools view -@ 16 -Sb ${line}.sam > ${line}.bam
```
Converting .bam to .sam saves the storage.
```shell
samtools sort -@ 16 ${line}.bam -o ${line}.sort.bam
stringtie -e -B -p 16 -G /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o $line/$line.FPKM.gtf ../align/${line}.sort.bam
```
Creating the index for the .bam.
```shell
samtools index ${line}.sort.bam
```
##  <font size=4>4   assembly and quantitate(stringtie)</font>
Reconstrcting all the isoforms from each genes and estimating their relative abundance.
```shell
stringtie -e -B -p 16 -G /media/hp/disk1/song/Genomes/${species}/Genes/genes.gtf -o $line/$line.FPKM.gtf ../align/${line}.sort.bam
``` 
<!--/TOC-->
