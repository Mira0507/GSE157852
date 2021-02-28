## STAR vs HISAT2 vs Salmon

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

#### 2-1. Tools for alignment/mapping

- HISAT2: http://daehwankimlab.github.io/hisat2
- Samtools: http://www.htslib.org/doc/#manual-pages
- STAR: https://github.com/alexdobin/STAR
- Salmon: https://salmon.readthedocs.io/en/latest
- bedtools: https://bedtools.readthedocs.io/en/latest
- gawk: https://www.gnu.org/software/gawk/manual/gawk.html

#### 2-2. Tools for counting and differential expression (DE) analysis

- DESeq2: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- Tximport: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
- Rsubread: https://bioconductor.org/packages/release/bioc/html/Rsubread.html


### 3. Reference files 

- Same references as [this project](https://github.com/Mira0507/seqc_comparison/blob/master/README.md)
- GENCODE GRCh38 (hg38) release 36 (v36)

### 4. STAR alignment

#### 4-1. Indexing

- Reused from [this project](https://github.com/Mira0507/seqc_comparison/blob/master/README.md)


#### 4-2. Alignment

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

### 5. HISAT2 Alignment

#### 5-1. Indexing 

- Reused from [this project](https://github.com/Mira0507/seqc_comparison/blob/master/README.md)

#### 5-2. Alignment 

- hisat2_align.sh

```bash


#!/bin/bash 

# Define directory and sample names
refdir=../SEQC/reference_GENCODE/hisat2_index/index   # reference directory
samples=(Mock-rep{1..3} SARS-CoV-2-rep{1..3})                 # sample names
outdir=hisat2_output                          # output directory
indir=rawdata                                 # input directory

mkdir $outdir 


for read in ${samples[*]} 

do 

    hisat2 -q -p 16 --seed 23 -x $refdir -U $indir/${read}.fastq.gz -S $outdir/$read.sam 

done 


```


#### 5-3. Converting SAM to BAM 

- uses samtools
- hisat2_samtobam.sh

```bash
#!/bin/bash

cd hisat2_output

input=(Mock-rep{1..3} SARS-CoV-2-rep{1..3})  

for f in ${input[*]}

do
    samtools view -bS -@ 16 $f.sam > $f.bam 
done 


cd ..
```

#### 5-4. Sorting 

- uses samtools
- hisat2_sort.sh

```bash
#!/bin/bash

cd hisat2_output

input=(Mock-rep{1..3} SARS-CoV-2-rep{1..3})  

for f in ${input[*]}

do
    samtools sort -@ 16 $f.bam -o $f.sorted.bam

done 


cd ..

# Delete SAM files after samtools run
```


### 6. Salmon mapping

#### 6-1. Indexing 

- Reused from [this project](https://github.com/Mira0507/seqc_comparison/blob/master/README.md)

#### 6-2. Mapping

- salmon_map.sh

```bash
#!/bin/bash

# Define name of input/output directory
in=rawdata
out=salmon_output


# Define index file directory
ind=../../SEQC/reference_GENCODE/salmon_index/gencode_index

# Define file names 
samples=(Mock-rep{1..3} SARS-CoV-2-rep{1..3})

mkdir $out

cd $in

for read in ${samples[*]}
do
    salmon quant -i $ind -l A --gcBias --seqBias -r ${read}.fastq.gz -p 16 --validateMappings -o ../$out/${read}.salmon_quant
done

cd ..
```

### 7. Counting and downstream analysis in R 

- [DE_EnsDb.Rmd](https://github.com/Mira0507/GSE157852/blob/master/DE_EnsDb.Rmd): DE analysis (STAR vs HISAT2 vs Salmon) with EnsDb 

- [DE_EnsDb.html](https://github.com/Mira0507/GSE157852): output of DE_EnsDb.Rmd

- [DE_OrDb.Rmd](https://github.com/Mira0507/GSE157852/blob/master/DE_OrgDb.Rmd): DE analysis (STAR vs HISAT2 vs Salmon) with OrgDb

- [DE_OrDb.html](https://github.com/Mira0507/GSE157852/blob/master/DE_OrgDb.html): output of DE_OrDb.Rmd

- [Ranking_FDR.Rmd](https://github.com/Mira0507/GSE157852/blob/master/Ranking_FDR.Rmd): Exploring factors affecting gene ranking difference in FDR, input from EnsDb 

- [Ranking_FDR.html](https://github.com/Mira0507/GSE157852/blob/master/Ranking_FDR.html): output of Ranking_FDR.Rmd

- [Ranking_LFC.Rmd](https://github.com/Mira0507/GSE157852/blob/master/Ranking_LFC.Rmd): Exploring factors affecting gene ranking difference in LFC, input from EnsDb 

- [Ranking_LFC.html](https://github.com/Mira0507/GSE157852/blob/master/Ranking_LFC.html): output of Ranking_LFC.Rmd

- [salmon_tpm_count.Rmd](https://github.com/Mira0507/GSE157852/blob/master/salmon_tpm_count.Rmd): Salmon TPM vs Count DE analysis by shrinkage

- [salmon_tpm_count.html](https://github.com/Mira0507/GSE157852/blob/master/salmon_tpm_count.html): output of salmon_tpm_count.Rmd

- [Ranking_tpmVScount_FDR.Rmd](https://github.com/Mira0507/GSE157852/blob/master/Ranking_tpmVScount_FDR.Rmd): TPM vs Count input in Salmon/DESeq2 DE analysis (no shrinkage)

- [Ranking_tpmVScount_FDR.html](https://github.com/Mira0507/GSE157852/blob/master/Ranking_tpmVScount_FDR.html): output of Ranking_tpmVScount_FDR.Rmd

- [Ranking_tpmVScount_LFC.Rmd](https://github.com/Mira0507/GSE157852/blob/master/Ranking_tpmVScount_LFC.Rmd): TPM vs Count input in Salmon/DESeq2 DE analysis (no shrinkage)

- [Ranking_tpmVScount_LFC.html](https://github.com/Mira0507/GSE157852/blob/master/Ranking_tpmVScount_LFC.html): output of Ranking_tpmVScount_LFC.Rmd

- [OrDb_vs_EnsDb.Rmd](https://github.com/Mira0507/GSE157852/blob/master/OrgDb_vs_EnsDb.Rmd): Comparing OrgDb vs EnsDb in AnnotationHub

- [OrDb_vs_EnsDb.html](https://github.com/Mira0507/GSE157852/blob/master/OrgDb_vs_EnsDb.html): output of OrDb_vs_EnsDb

- [h_FDR_RankDiff.csv]
