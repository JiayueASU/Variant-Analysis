# Variant Analysis
## Quality Assessment
### FastQC
Download Page: <http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc>

Unzip FastQC: `unzip fastqc_v0.11.9.zip`

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

For more details about using Trimmomatic: `java -jar Trimmomatic-0.39.jar -h`

## Read Alignment
### BWA-mem
Download Page:
<https://sourceforge.net/projects/bio-bwa/files/>

Unzip bwa-0.7.17.tar.bz2: `tar -jxvf bwa-0.7.17.tar.bz2`



