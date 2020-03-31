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

After 8489.8 sec, aln-pe.sam file is generated.

## Variant Identification

### GATK

Tutorial: <http://www.bio-info-trainee.com/3144.html>

First download GATK: `cd /data/notebook/Jerry/Tools`

`wget https://github.com/broadinstitute/gatk/releases/download/4.0.2.1/gatk-4.0.2.1.zip`

Then unzip gatk-4.0.2.1.zip: `unzip gatk-4.0.2.1`

Set environmental variables: `echo export PATH=$PATH:/data/notebook/Jerry/Tools/gatk-[4.0.2.1](4.0.2.1) >> ~/.bashrc ` `source ~/.bashrc`



### Samtools

First download samtools: `wget -c https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2`

Then unzip samtools-1.9: `tar jxvf samtools-1.9.tar.bz2`

Set Configure: `./configure --prefix=/data/notebook/Jerry/Tools/samtools-1.9`

Then `make` `make install`



