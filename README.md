# RNAseq analysis for Homo sapiens

## Login Midway2

```
ssh yourID@midway2.rcc.uchicago.edu
quota
rcchelp balance
rcchelp usage --byjob
module load avail
```

## Transfer files: Mac OS X Connect to Server utility

https://rcc.uchicago.edu/docs/data-transfer/index.html

*If you are on campus*
home:            smb://midwaysmb.rcc.uchicago.edu/homes
project2:        smb://midwaysmb.rcc.uchicago.edu/project2
scratch/midway:  smb://midwaysmb.rcc.uchicago.edu/midway2-scratch

ADLOCAL\yourID
Registered User

*If you are off campus - use Globus*
https://rcc.uchicago.edu/docs/data-transfer/index.html#globus-online

https://globus.rcc.uchicago.edu 
Sign in
-> enter ucrcc#midway for at the "Search all endpoints" box (Upper right corner)
-> File Manager (Top left panel)
-> Tansfer or sync file
-> Globus Connect Personal (links to your local computer)

## Establish your PATH -> /home/yourID/local/bin

```
mkdir ~/local/bin
```
```
## Installations
module avail 
module list
	
**Then build your PATH **
* echo 'export PATH=~/local/bin:$PATH' >> ~/.bashrc
* source ~/.bashrc
```

## Data processing

### 0. Merge replicate fastq files

```
cat file*.fastq > bigfile.fastq
cat file1.fastq file2.fastq > bigfile.fastq

*21
cat ~/scratch-midway2/fastq_EC-CC1/EC-CC-21_S1_R1_001.fastq.gz ~/scratch-midway2/fastq_EC-CC2/EC-CC-21_S1_R1_001.fastq.gz > ~/scratch-midway2/fastq_EC-CC3/EC-CC-21_S1_R1_001.fastq.gz | cat ~/scratch-midway2/fastq_EC-CC1/EC-CC-21_S1_R2_001.fastq.gz ~/scratch-midway2/fastq_EC-CC2/EC-CC-21_S1_R2_001.fastq.gz > ~/scratch-midway2/fastq_EC-CC3/EC-CC-21_S1_R2_001.fastq.gz 
```

### 2. Map RNA-seq reads to a reference genome (STAR)
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4631051/
https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/master/lessons/03_alignment.md

#### 2.1: Generate genome indices

```
*sinteractive --exclusive --partition=broadwl --nodes=2 --time=05:00:00*

**ensembl reference genome **
https://useast.ensembl.org/Homo_sapiens/Info/Index

Homo_sapiens.GRCh38.99.gtf.gz
Homo_sapiens.GRCh38.dna.toplevel.fa.gz
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gunzip ~/scratch-midway2/RefGenome/Homo_sapiens.GRCh38.99.gtf.gz
gunzip ~/scratch-midway2/RefGenome/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
gunzip ~/scratch-midway2/RefGenome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

*NOT RUN - examine the number of references in the genome file*
grep "^>" ~/scratch-midway2/RefGenome/Homo_sapiens.GRCh38.dna.toplevel.fa | wc -l
639
a genome with a large > 5,000 number of references (chrosomes/scaffolds), reduce the --genomeChrBinNbits
--genomeChrBinNbits 15 *use when genome is big*
--genomeSAindexNbases 10 *use when genome is small*


*run the following*
STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir ~/scratch-midway2/EC-CC/1.star_index/ensembl \
--genomeFastaFiles ~/scratch-midway2/RefGenome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile ~/scratch-midway2/RefGenome/Homo_sapiens.GRCh38.99.gtf \
--sjdbOverhang 99 

```

#### 2.2: Align reads 

```
*sinteractive --exclusive --partition=broadwl --nodes=2 --time=05:00:00*
cd ~/scratch-midway2/EC-CC/2.star/star2

mkdir star_EC-CC-All | mkdir star_EC-CC-All01 | mkdir star_EC-CC-All02 | mkdir star_EC-CC-All03 | mkdir star_EC-CC-All04 | mkdir star_EC-CC-All05 | mkdir star_EC-CC-All06 | mkdir star_EC-CC-All07 | mkdir star_EC-CC-All08 |
*done*
```

```
STAR --genomeDir ~/scratch-midway2/EC-CC/1.star_index/ensembl  \
--runThreadN 6 \
--readFilesIn ~/scratch-midway2/EC-CC/0.fastq/fastq_EC-CC3/EC-CC-21_R1.fastq.gz ~/scratch-midway2/EC-CC/0.fastq/fastq_EC-CC3/EC-CC-21_R2.fastq.gz \
--outFileNamePrefix ~/scratch-midway2/EC-CC/2.star/star_EC-CC-All01/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--limitBAMsortRAM 1232875179 \
--readFilesCommand zcat 
*21 done*

*Examine results*
cd ~/scratch-midway2/EC-CC/2.star/star2/star_EC-CC-All01
less Log.final.out
cat Log.progress.out


cd ~/scratch-midway2/EC-CC/2.star
cp star_EC-CC-All01/Aligned.sortedByCoord.out.bam star_EC-CC-All/Aligned.sortedByCoord.EC-CC-21.out.bam
```


### 3. Count transcripts 

```
mkdir ~/scratch-midway2/EC-CC/3.featureCounts/


featureCounts -T 4 \
-a ~/scratch-midway2/RefGenome/Homo_sapiens.GRCh38.99.gtf \
-o ~/scratch-midway2/EC-CC/3.featureCounts/featureCounts_EC-CC.txt \
-g gene_id \
~/scratch-midway2/EC-CC/2.star/star_EC-CC-All/*.out.bam

featureCounts -T 4 \
-a ~/scratch-midway2/RefGenome/Homo_sapiens.GRCh38.99.gtf \
-o ~/scratch-midway2/EC-CC/3.featureCounts/featureCounts_new2.txt \
-g gene_name \
~/scratch-midway2/EC-CC/2.star/star_EC-CC-All/*.out.bam

An example of attributes included in your GTF annotation is 'gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; transcript_support_level "1";' 
```


## MultiQC Pipeline

https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/master/lessons/06_multiQC.md

### Files
.zip files from FastQC
.Log.final.out files from STAR
.summary file from featurecounts

cp ~/scratch-midway2/EC-CC/2.star/star_EC-CC-All01/Log.final.out ~/scratch-midway2/EC-CC/2.star/star_EC-CC-log/Log.final-21.out 

### Install MultiQC
https://multiqc.info/docs/
pip install multiqc
pip --update multiqc

### Run

```
cd ~/"Google Drive"/UC/Projects/scRNAseq_EC-CC/MultiQC

multiqc -n multiqc_report_rnaseq \
~/"Google Drive"/UC/Projects/scRNAseq_EC-CC/MultiQC/*zip \
~/"Google Drive"/UC/Projects/scRNAseq_EC-CC/MultiQC/*Log.final.out \
~/"Google Drive"/UC/Projects/scRNAseq_EC-CC/MultiQC/featureCounts_EC-CC.txt.summary
```






