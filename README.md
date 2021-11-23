# nanotei
A pipeline for detection of non-reference transposon insertions using Nanopore data

## Installation
nanotei depends on several programs and python modules:
```
  - python =3.7.5
  - pysam =0.17.0
  - samtools =1.14
  - minimap2 =2.22
  - biopython =1.76
  - statistics =1.0.3.5
  - pandas =0.24.2

``` 

The best way to install nanotei is to use nanotei.yaml file to create conda environment
```
 conda env create --file nanotei.yaml
```

## How to run
nanotei is written in python and can be run from command line.

```
  usage: nanotei.py [-h] [-q MAP_Q] [-mlc MIN_LEN_CLIPPED]
                  [-md MIN_MERGE_DISTANCE] [-samtp SAMTOOLS_PATH]
                  [-mm2 MINIMAP2_PATH] [-bt BEDTOOLS]
                  [-ovt OVERLAP_WITH_ORIGTE] [-mpv MINPVALUE] [--fpv]
                  [-mrs MIN_READ_SUPPORT] [--bed]
                  bam_file fastq genome_fasta outfolder bed_TE outtab

```

Script to find regions of insertions using Nanopore reads

### positional arguments:
  **bam_file**              path to the sorted bam file with mapped reads to the genome
  
  **fastq**                 path to fastq file of reads
  
  **genome_fasta**          path to fasta file of genome reference
  
  **outfolder**             path to the output folder
  
  **bed_TE**                path to the bed file with annotated TE
  
  **outtab**                name of the output table
  

### optional arguments:
  *-h*, *--help*            show this help message and exit
  
  *-q* MAP_Q, *--map_q* MAP_Q  minimum mapping quality
  
  *-mlc* MIN_LEN_CLIPPED, *--min_len_clipped* MIN_LEN_CLIPPED
                        minimum length of the clipped part
                        
  *-md* MIN_MERGE_DISTANCE, *--min_merge_distance* MIN_MERGE_DISTANCE
                        minimum distance between two reference positions of
                        clipping to be merged into one
                        
  *-samtp* SAMTOOLS_PATH, *--samtools_path* SAMTOOLS_PATH
                        path to samtools program
                        
  *-mm2* MINIMAP2_PATH, *--minimap2_path* MINIMAP2_PATH
                        path to minimap2 program
                        
  *-bt* BEDTOOLS, *--bedtools* BEDTOOLS
                        path to minimap2 program
                        
  *-ovt* OVERLAP_WITH_ORIGTE, *--overlap_with_origTE* OVERLAP_WITH_ORIGTE
                        minimum portion of origin TE overlapped with clipped
                        repa alignemnt
                        
  *-mpv* MINPVALUE, *--minpvalue* MINPVALUE
                        minimum pvalue to filter (it is used with --fpv). 
                        pvalue here reflects the probability that the TEI is OK (so usually it should be > 0.05) 
  *--fpv*                 
  
  *-mrs* MIN_READ_SUPPORT, *--min_read_support* MIN_READ_SUPPORT
                        minimum read support for filtering
                        
  *--bed*                 output bed file with TEIs
  
  ## input files
  
  **bam file** can be obtained by read mapping to the reference genome by minimap2 (`minimap2 -ax map-ont -t 150  genome.fasta reads.fastq > file.sam`) followed by sam to bam conversion (`samtools view -Sb file.sam > file.bam`) and sorting (`samtools sort -o sorted_file.bam file.bam`)
  
  **bed_TE** nanotei requires 4-column bed file without header: chromosome_id, start, end, TEid. 
  
  
  ## example of the command
### get copy of githup repository
### change directory to test
```
git clone https://github.com/Kirovez/nanotei.git
cd ./nanotei/test

```
### 1. download TAIR10 genome from NCBI and move it to test directory (it must be current directory)
TAIR10.1 repository: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001735.4#/def_asm_Primary_Assembly

unpack tar and unzip .gz in the folder. Then move GCF_000001735.4_TAIR10.1_genomic.fna file to test directory (current directory). Examples of the commands:
    
```
tar -xvf genome_assemblies_genome_fasta.tar

gunzip ./ncbi-genomes-2021-11-23/GCF_000001735.4_TAIR10.1_genomic.fna.gz

mv ./ncbi-genomes-2021-11-23/GCF_000001735.4_TAIR10.1_genomic.fna .

```
        
### 2. download Nanopore reads of ddm1 plant

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR147/094/SRR14765194/SRR14765194_1.fastq.gz 
## gunzip
gunzip SRR14765194_1.fastq.gz

```
### 3. TE annotation file that was build by Panda and Slotkin, 2020 must be already in the test directory. The chromosome ids were changed according to the NCBI chromosome ids. The file is available on GitHub:
```
wget https://raw.githubusercontent.com/Kirovez/nanotei/master/test/TE_31189_Slotkin_NCBI_id.bed
```

### 4. map the reads to the genome using minimap2, then convert sam file to bam, sort and index:
#### mapping
minimap2 -ax map-ont -t 100 GCF_000001735.4_TAIR10.1_genomic.fna SRR14765194_1.fastq > SRR14765194_vs_TAIR10.sam
#### convert sam to bam

```
samtools view -Sb SRR14765194_vs_TAIR10.sam > SRR14765194_vs_TAIR10.bam 
```

#### sort bam file

``` 
samtools sort -o sorted_SRR14765194_vs_TAIR10.bam -@ 100 SRR14765194_vs_TAIR10.bam 
```
#### index
```
samtools index -@ 100 sorted_SRR14765194_vs_TAIR10.bam
```
### 5. run nanotei (please, use absolute path to input files or ./ before file name ). The resulted files will appear in the current directory (it will be changed in later versions).
```
python3 ../nanotei.py ./sorted_SRR14765194_vs_TAIR10.bam ./SRR14765194_1.fastq GCF_000001735.4_TAIR10.1_genomic.fna ./out_nanotei ./TE_31189_Slotkin_NCBI_id.bed SRR14765194.nanotei --fpv --bed --minpvalue 0.05 -ovt 0.3

```

#### an example ot output table and bed file can be found in out_example folder. Note: TEIs in output folder for ddm1 includes ddm1 specific TEIs and TEIs that resulted from assembly artefacts (~44 TEIs, see nanotei paper for details)  
  
  ## Output table
  The output nanotei table contains the following columns:
  
  **TEIcoord** - the coordinates of TEI: chromosome:startTEI..endTEI 
  
  **supp.reads** - the number of reads supporting this TEI (reads with clipped ends + reads with insertions)
  
  **pvalue** - the probability that this insertion is real based on the number of reads supporting TEI and the genome coverage distribution (so 1 is absolutely OK and 0 - is bad TEI).
  
  **isTEIregioIsOK** - it is always TRUE because nanotei is filtered out all the regions with too high coverage based on expectation from genome coverage distribution
  
  **TE_reads** - it shows the number of mapped sequences of clipped parts and insertion parts of the reads overlaped with TEs from bed file.
  
  **chosenTE** - a TE with the highest number of mapped clipped parts and insertion parts of the reads
  
  **TE.chr** - chromosome name of chosen TE
  
  **TE.start** - chromosome start of chosen TE
  
  **TE.end** - chromosome end of chosen TE
  
  ## nanotei citation
  
  We are currently preparing the manuscript.
  
