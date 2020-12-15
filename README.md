## STAR vs HISAT2 

### 1. Rawdata 

- GEO number: [GSE157852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157852)
- SRA number: [SRP282132](https://www.ncbi.nlm.nih.gov/sra?term=SRP282132)
- Sample description: Neurological complications are common in patients with COVID-19. While SARS-CoV-2, the causal pathogen of COVID-19, has been detected in some patient brains, its ability to infect brain cells and impact their function are not well understood, and experimental models using human brain cells are urgently needed. Here we investigated the susceptibility of human induced pluripotent stem cell (hiPSC)-derived monolayer brain cells and region-specific brain organoids to SARS-CoV-2 infection. We found modest numbers of infected neurons and astrocytes, but greater infection of choroid plexus epithelial cells. We optimized a protocol to generate choroid plexus organoids from hiPSCs, which revealed productive SARS-CoV-2 infection that leads to increased cell death and transcriptional dysregulation indicative of an inflammatory response and cellular function deficits. Together, our results provide evidence for SARS-CoV-2 neurotropism and support use of hiPSC-derived brain organoids as a platform to investigate the cellular susceptibility, disease mechanisms, and treatment strategies for SARS-CoV-2 infection. Bulk RNA-seq of choroid plexus organoids (CPOs) was performed on mock 72 hours post-infection (hpi), SARS-CoV-2 24 hpi, and SARS-CoV-2 72 hpi samples. All conditions were profiled in triplicate.
- Rawdata download: [SRAtoolkit](https://github.com/Mira0507/using_SRA/blob/master/README.md)

#### 1-1. Downloading SRA files 

- rawdata_download.sh

```bash

#!/bin/bash


# Numbers of SRA data:
# SRR12626040: SARS-CoV-2 72 hpi rep1	
# SRR12626041: SARS-CoV-2 72 hpi rep2
# SRR12626042: SARS-CoV-2 72 hpi rep3
# SRR12626034: Mock 72 hpi rep1
# SRR12626035: Mock 72 hpi rep2
# SRR12626036: Mock 72 hpi rep3

mkdir rawdata

cd rawdata

SRR_number=(34 35 36 40 41 42)

for x in ${SRR_number[*]}
do
    prefetch SRR126260$x

done


cd ..
```

#### 1-2. Converting to FASTQ files 

- sra_fastq.sh

```bash
#!/bin/bash

# Numbers of SRA data:
# SRR12626040: SARS-CoV-2 72 hpi rep1	
# SRR12626041: SARS-CoV-2 72 hpi rep2
# SRR12626042: SARS-CoV-2 72 hpi rep3
# SRR12626034: Mock 72 hpi rep1
# SRR12626035: Mock 72 hpi rep2
# SRR12626036: Mock 72 hpi rep3


cd rawdata

SRR_number=(34 35 36 40 41 42)

for x in ${SRR_number[*]}
do
     fastq-dump --split-files SRR126260$x/SRR126260$x.sra 

done


cd ..

```

#### 1-3. Renaming 

- name_change.sh

```bash

#!/bin/bash

# SRR12626040: SARS-CoV-2 72 hpi rep1	
# SRR12626041: SARS-CoV-2 72 hpi rep2
# SRR12626042: SARS-CoV-2 72 hpi rep3
# SRR12626034: Mock 72 hpi rep1
# SRR12626035: Mock 72 hpi rep2
# SRR12626036: Mock 72 hpi rep3

cd rawdata

mv SRR12626040_1.fastq SARS-CoV-2-rep1.fastq
mv SRR12626041_1.fastq SARS-CoV-2-rep2.fastq
mv SRR12626042_1.fastq SARS-CoV-2-rep3.fastq
mv SRR12626034_1.fastq Mock-rep1.fastq
mv SRR12626035_1.fastq Mock-rep2.fastq
mv SRR12626036_1.fastq Mock-rep3.fastq

cd ..
```

#### 1-4. Gunzip 

- fastq_gz.sh 

```bash
#!/bin/bash



cd rawdata

gzip -v -k *.fastq 



cd ..
```

### 2. Conda environment 

- Reused from [this project](https://github.com/Mira0507/seqc_comparison/blob/master/README.md)

### 3. Reference files 

- Same references as [this project](https://github.com/Mira0507/seqc_comparison/blob/master/README.md)
- GENCODE GRCh38 (hg38) release 36 (v36)

### 3. STAR alignment

#### 3-1. Indexing

- Reused from [this project](https://github.com/Mira0507/seqc_comparison/blob/master/README.md)



#### 3-2. Alignment

- star_alignment.sh

```bash
#!/bin/bash

# Define directory names 
outdir=star_output   # directory storing output files
indir=../rawdata        # input directory (fastq files)
refdir=/home/mira/Documents/programming/Bioinformatics/SEQC/reference_GENCODE     # reference directory (absolute path needed!)
indexdir=star_index    # index directory
genome=*genome.fa      # reference file
gtf=*.gtf              # GTF file


samples=(Mock-rep{1..3} SARS-CoV-2-rep{1..3})

mkdir $outdir 
cd $outdir

for read in ${samples[*]} 

do 
    STAR --runThreadN 16 --runMode alignReads --genomeDir $refdir/$indexdir --readFilesCommand zcat --sjdbGTFfile $refdir/$gtf -sjdbOverhang 100 --readFilesIn $indir/${read}.fastq.gz --outFileNamePrefix $read --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic --chimOutType Junctions

done

cd ..

# --readFilesIn: For paired-end reads, use comma separated list for read1, followed by space, followed by comma separated list for read2
# --readFilesCommand: required when the input files are .gzip format (e.g. --readFilesCommand zcat, --readFilesCommand gunzip -c, or --readFilesCommand bunzip2)
# Note: the output files are generated in the current directory
```

### 4. HISAT2 Alignment

#### 4-1. Indexing 

- Reused from [this project](https://github.com/Mira0507/seqc_comparison/blob/master/README.md)

#### 4-2. Alignment 

- hisat2_align.sh

```bash

```
