# Variant Analysis
## Quality Assessment
### FastQC
Download Page: <http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc>

Unzip FastQC: `unzip fastqc_v0.11.9.zip`

Set environmental variables: `echo export PATH=$PATH:/data/notebook/Jerry/Tools/FastQC/fastqc >> ~/.bashrc`

`source ~/.bashrc`

Usage of FastQC: `fastqc -o ./ -t 6 V300035135_L03_531_1.clean.fq.gz ...`

Or you can process multiple .fq.gz files in batch with a .sh script:
```javascript
for id in *fastq
do
echo $id
/data/notebook/Jerry/Test/Data/Data20200323 $id
Done
```

Explanation: 
-o --outdir     
> Create all output files in the specified output directory. Please note that this directory must exist as the program will not create it. If this option is not set then the output file for each sequence file is created in the same directory as the sequence file which was processed.

-t --threads    
> Specifies the number of files which can be processed simultaneously. Each thread will be allocated 250MB of memory so you shouldn't run more threads than your available memory will cope with, and not more than 6 threads on a 32 bit machine.

After the execution, an .html file will be generated. This file is the fastqc report of the associate fq data. The explanation of this report can be found at: 
<http://www.bio-info-trainee.com/95.html>

### Trimmomatic
Download Page: <http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip>

Unzip Trimmomatic: `unzip Trimmomatic-0.39.zip`

The tutorial of using Trimmomatic can be found at: 
<http://www.bio-info-trainee.com/1958.html>

<http://www.usadellab.org/cms/?page=trimmomatic>

Installation of Java: `yum install [java-1.8.0-openjdk.x86_64](java-1.8.0-openjdk.x86_64) `

For more details about using Trimmomatic: `java -jar Trimmomatic-0.39.jar -h`

Usage of Trimmomatic (Paired End): `java -jar trimmomatic-0.39.jar PE -phred33 /data/notebook/Jerry/Test/Input/Data20200323/V300035135_L03_531_1.clean.fq.gz /data/notebook/Jerry/Test/Input/Data20200323/V300035135_L03_531_2.clean.fq.gz /data/notebook/Jerry/Test/Reference/Output_V300035135_L03_531_1_paired.clean.fq.gz /data/notebook/Jerry/Test/Reference/Output_V300035135_L03_531_1_unpaired.clean.fq.gz /data/notebook/Jerry/Test/Reference/Output_V300035135_L03_531_2_paired.clean.fq.gz /data/notebook/Jerry/Test/Reference/Output_V300035135_L03_531_2_unpaired.clean.fq.gz /data/notebook/Jerry/Tools/Trimmomatic-0.39/adapters/ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

Output_paired: Usable data

Output_unpaired: Removed adapters, leading low quality, and trailing low quality

Usage of Trimmomatic (Single End): `java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

After the execution, we could obtain paired data for next step.

## Read Alignment
### BWA-mem
Download Page:
<https://sourceforge.net/projects/bio-bwa/files/>

Tutorial: <https://www.jianshu.com/p/1552cc6ac3be>

Unzip bwa-0.7.17.tar.bz2: `tar -jxvf bwa-0.7.17.tar.bz2`

Install bwa:`cd bwa-0.7.17` `make`

After `make`, execute bwa file: `./bwa`

Set environmental variables: `export PATH=$PATH:/data/notebook/Jerry/Tools/bwa-0.7.17/`

`echo export PATH=$PATH:/data/notebook/Jerry/Tools/bwa-0.7.17 >> ~/.bashrc` `source ~/.bashrc`

Create a new folder for mm10.fasta: `cd /data/notebook/Jerry/Test/Reference`

Create a new folder at Output: `mkdir bwa_index`

`mkdir mm10` `cd mm10`

Download mm10.fast: `wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz`

cd to: `/data/notebook/Jerry/Test/Reference/mm10/Mus_musculus/UCSC/mm10/Sequence/BWAIndex`

Copy genome.fa to the Output folder: `cp genome.fa /data/notebook/Jerry/Test/Output/bwa_index/`

Generate index sequence: `bwa index genome.fa`

After 613 iterations, five new files are generated: genome.fa.amb, genome.fa.ann, genome.fa.bwt, genome.fa.pac, and genome.fa.sa.

Use BWA-mem to obtain .sam file: `bwa mem genome.fa /data/notebook/Jerry/Test/Input/Data20200323/V300035135_L03_531_1.clean.fq.gz /data/notebook/Jerry/Test/Input/Data20200323/V300035135_L03_531_2.clean.fq.gz > aln-pe.sam`

After 11142.574 sec, aln-pe.sam file is generated.

## Variant Identification

### Samtools

First download samtools: `wget -c https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2`

Then unzip samtools-1.9: `tar jxvf samtools-1.9.tar.bz2`

Set Configure: `./configure --prefix=/data/notebook/Jerry/Tools/samtools-1.9`

If it shows a bug like: "fatal error: curses.h: No such file or directory", please install libncurses5-dev (ubuntu) or curse-devel (centos).

Then `make` `make install`

If it shows a bug on "htslib-1.9", install this package.

Generate a new folder at Output: `mkdir samtools_bam`

Copy the generated .sam file to samtools_bam: `cp aln-pe.sam ~/samtools_bam`

Generate .bam file (15 min): `samtools view -bS aln-pe.sam > aln-pe.bam`

Sort the generated .bam file (25 min): `samtools sort -n aln-pe.bam -o aln-pe.sort.bam`

### GATK

Tutorial: <http://www.bio-info-trainee.com/3144.html>

First download GATK: `cd /data/notebook/Jerry/Tools`

`wget https://github.com/broadinstitute/gatk/releases/download/4.0.2.1/gatk-4.0.2.1.zip`

Then unzip gatk-4.0.2.1.zip: `unzip gatk-4.0.2.1`

Set environmental variables: `echo export PATH=$PATH:/data/notebook/Jerry/Tools/gatk-[4.0.2.1](4.0.2.1) >> ~/.bashrc ` `source ~/.bashrc`

Use GATK to mark duplicates (18 min): `java -jar gatk-package-4.0.2.1-local.jar MarkDuplicates \-I /data/notebook/Jerry/Test/Output/samtools_bam/aln-pe.sort.bam -O aln-pe.sort.markdup.bam -M aln-pe.sort.markdup.bam.metrics`

Sort the generated .bam file again: `samtools sort aln-pe.sort.markdup.bam -o aln-pe.sort1.markdup.bam`

Create index for sorted marked .bam file to generate aln-pe.sort1.markdup.bam.bai: `samtools index aln-pe.sort1.markdup.bam`

Genreate a new folder at Reference, named GATK. Then download .vcf files from: `wget ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz`

Switch to /data/notebook/Jerry/Test/Output/bwa_index, generate .fai file: `samtools faidx genome.fa`

Generate .dict file based on genome.fa: `samtools dict genome.fa > genome.dict`

Do recalibration at the folder of gate: `java -jar gatk-package-4.0.2.1-local.jar BaseRecalibrator -R /data/notebook/Jerry/Test/Output/bwa_index/chromosomes.fa -I /data/notebook/Jerry/Test/Output/samtools_bam/aln-pe-0406-sort1-markdup.bam -O /data/notebook/Jerry/Test/Output/samtools_bam/aln-pe-0406-sort1-markdup.recal.table --known-sites /data/notebook/Jerry/Test/Reference/mm10SNP/mgp.v3.snps.rsIDdbSNPv137.vcf --known-sites /data/notebook/Jerry/Test/Reference/mm10SNP/mgp.v3.indels.rsIDdbSNPv137.vcf`

Use IndexFeatureFile to generate vcd.idx file (This step solves the problem >= 1 but = 0): `java -jar gatk-package-4.0.2.1-local.jar IndexFeatureFile -F /data/notebook/Jerry/Test/Reference/GATK/00-All.vcf `

Check out this tutorial: <https://www.jianshu.com/p/fe0c876563b0>

Check out the detail on .bam file: `samtools view -H aln-pe-0406-sort1-markdup.bam`

Merge all the .fa file at Chromosomes to one .fa file: `cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chrM.fa chrX.fa chrY.fa > chromosomes.fa`