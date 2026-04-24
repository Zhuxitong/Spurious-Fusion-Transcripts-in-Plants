# 01. Sequencing Data Preprocessing

## Purpose
Preprocess long-read sequencing datasets, including PacBio Iso-Seq and Nanopore direct RNA/direct cDNA data, for downstream fusion transcript analyses.

## Related Figure
- Supplementary Figure S1

---

## A. PacBio Iso-Seq Preprocessing

### Purpose
Generate full-length non-chimeric (FLNC) reads from raw PacBio subreads for downstream fusion transcript identification and breakpoint analyses.

### Input
- Raw PacBio subreads BAM file: `isoseq.subreads.bam`
- Associated index/metadata files: `.pbi`, `.xml`
- Primer file: `test_data/primers.fasta`

### Output
- CCS reads: `isoseq.ccs.bam`
- Full-length reads: `isoseq.fl.bam`
- FLNC reads: `isoseq.flnc.bam`
- FLNC FASTA: `isoseq.flnc.fasta`

### Workflow

#### Step A1. Generate Circular Consensus Reads (CCS)
```bash
ccs isoseq.subreads.bam isoseq.ccs.bam \
-j 3 \
--min-length 50 \
--report-file isoseq.ccs.report \
--min-rq 0.9
````

#### Step A2. Remove Primers and Identify Full-Length Reads
```bash
lima isoseq.ccs.bam test_data/primers.fasta isoseq.fl.bam \
--isoseq -j 5
```

#### Step A3. Refine Full-Length Non-Chimeric Reads (FLNC)
```bash
isoseq refine isoseq.fl.primer_5p--primer_3p.bam \
test_data/primers.fasta \
isoseq.flnc.bam \
--require-polya -j 5
```

#### Step A4. Convert FLNC BAM to FASTA
```bash
samtools fasta isoseq.flnc.bam > isoseq.flnc.fasta
```

---

## B. Nanopore Direct RNA / Direct cDNA Preprocessing
### Purpose
Perform quality control, error correction, and genome alignment for Nanopore direct RNA or direct cDNA reads prior to downstream fusion transcript analysis.

### Input
* Raw Nanopore FASTQ file: `dRNA.fastq.gz`
* Matched short-read RNA-seq FASTQ file for correction: `short_reads.fq.gz`
* Reference genome: `Omh63.genomic.fa`

### Output
* Quality-filtered reads: `dRNA.qc.fastq.gz`
* Corrected reads: `dRNA.correct.fasta`
* Genome-aligned BAM file: `dRNA.2mh.bam`

### Workflow
#### Step B1. Summarize Raw Read Statistics
```bash
seqkit stats -b dRNA.fastq.gz
```

#### Step B2. Perform Quality Filtering
```bash
chopper -l 150 -q 7 -i dRNA.fastq.gz | pigz - > dRNA.qc.fastq.gz
```

#### Step B3. Build the Short-Read BWT Index for FMLRC2 Correction
```bash
pigz -c -d short_reads.fq.gz | \
awk 'NR % 4 == 2' | \
sort | \
tr NT TN | \
ropebwt2 -LR | \
tr NT TN | \
msbwt2-convert short_reads_fmlrc2.npy
```

#### Step B4. Correct Nanopore Reads Using FMLRC2
```bash
fmlrc2 -C 10 -t 2 \
short_reads_fmlrc2.npy \
dRNA.qc.fastq.gz \
dRNA.correct.fasta
```

#### Step B5. Align Corrected Reads to the Reference Genome
```bash
minimap2 -ax splice -uf -k14 -t 2 Omh63.genomic.fa dRNA.correct.fasta | \
samtools sort --threads 2 - | \
samtools view --threads 2 -b - > dRNA.2mh.bam
```

### Notes
* `short_reads.fq.gz` refers to matched short-read RNA-seq data from the same tissue, which are used to improve Nanopore read correction.
* `Omh63.genomic.fa` is the MH63 reference genome and can be downloaded from: http://rice.hzau.edu.cn/rice_rs3/download_ext/MH63RS3.fasta.gz



# 02. Fusion Transcript Identification Using the Optimized JAFFAL Pipeline
## Purpose
Identify putative fusion transcripts from long-read sequencing data using a modified JAFFAL (v2.3) pipeline adapted for plant genomes and annotation formats.
## Related Figures
- Figure 1b–d
- Supplementary Figures S2, S3, S4a–b, S14b

---

## A. Install the Optimized JAFFAL Package
### Purpose
Install the customized JAFFAL package included in the `scripts/` directory. This version contains modifications to several core C++ components to improve compatibility with plant genomes and GFF/GTF annotations.

### Reference
- Original JAFFAL repository:  
  https://github.com/Oshlack/JAFFA
- Installation guide:  
  https://github.com/Oshlack/JAFFA/wiki/HowToSetUpJAFFA

### Workflow
#### Step A1. Unpack and install JAFFAL

```bash
unzip JAFFA-version-2.3.zip
cd JAFFA-version-2.3
./install_linux64.sh
````

---

## B. Build a Plant-Compatible JAFFAL Reference Database
### Purpose
JAFFAL provides prebuilt references for human and selected animal genomes, but not for plants. Therefore, a custom database must be generated using a plant reference genome FASTA file and corresponding annotation file.

### Input
* Reference genome FASTA: `MH63RS3.fasta.gz`
* Annotation file: `MH63RS3_nonTE.gff3.tar.gz`

### Workflow
#### Step B1. Create database directory and download reference files
```bash
mkdir JAFFAL_DB_MH63RS3
cd JAFFAL_DB_MH63RS3
wget http://rice.hzau.edu.cn/rice_rs3/download_ext/MH63RS3.fasta.gz
wget http://rice.hzau.edu.cn/rice_rs3/download_ext/MH63RS3_nonTE.gff3.tar.gz
```

#### Step B2. Rename genome file and convert annotation format
```bash
mv MH63RS3.fasta Omh63.fa
gffread MH63RS3_nonTE.gff3 -T -o Omh63.gtf
```

#### Step B3. Generate transcript FASTA
```bash
gffread Omh63.gtf -g Omh63.fa -w Omh63_RS3.fasta
```

#### Step B4. Reformat transcript FASTA headers
```bash
cat Omh63_RS3.fasta | \
awk '{if(/^>/){print $1}else{print}}' | \
~/PATH_to_JAFFAL/JAFFA-version-2.3/tools/bin/reformat \
fastawrap=0 in=stdin.fa out=stdout.fa > Omh63_RS3.fa
```

#### Step B5. Generate GenePred annotation table
```bash
gtfToGenePred Omh63.gtf Omh63.gpd -genePredExt

awk '{print "1\t"$0}' Omh63.gpd | \
csvtk_tH add-header \
-n '#bin,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,score,name2,cdsStartStat,cdsEndStat,exonFrames' \
> Omh63_RS3.tab
```

#### Step B6. Generate exon BED and masked genome
```bash
grep -P "\texon\t" Omh63.gtf | \
awk '{print $1"\t"$4-1"\t"$5"\t"$10"\t.\t"$7}' | \
sortBed > Omh63_RS3.bed

bedtools maskfasta \
-fi Omh63.fa \
-fo Masked_Omh63.fa \
-bed Omh63_RS3.bed
```

#### Step B7. Build Bowtie2 indices
```bash
bowtie2-build Omh63_RS3.fa Omh63_RS3
bowtie2-build Masked_Omh63.fa Masked_Omh63
```

#### Step B8. Build BLAST database
```bash
makeblastdb -in Omh63_RS3.fa -dbtype nucl -out Omh63_RS3_blast
```

#### Step B9. Create placeholder known fusion file

```bash
touch known_fusions.txt
```

### Output
A complete JAFFAL-compatible plant reference database analogous to the human `hg38_genCode22` package.

---

## C. Validate the Custom JAFFAL Database
### Purpose
Test whether the customized JAFFAL database functions correctly using simulated fusion transcript sequences.
### Input
* Simulated fusion transcripts generated by FusionTranscriptSimulator

### Reference
[https://github.com/FusionSimulatorToolkit/FusionSimulatorToolkit/wiki](https://github.com/FusionSimulatorToolkit/FusionSimulatorToolkit/wiki)

### Workflow
#### Step C1. Generate simulated fusion transcripts
```bash
~/PATH_to_FusionTranscriptSimulator/FusionSimulatorToolkit-master/FusionTranscriptSimulator \
MH63RS3_nonTE.gff3 Omh63.fa 20 > test.fasta
```

#### Step C2. Run JAFFAL on test data
```bash
$JAFFAL_HOME/tools/bin/bpipe run -n 1 \
-p genome=Omh63 \
-p annotation=RS3 \
-p refBase=~/PATH_to_JAFFAL_DB/ \
-p outputName=test \
refBase=~/PATH_to_JAFFAL_DB/JAFFAL_DB_MH63RS3 \
~/PATH_to_JAFFAL_DB/JAFFA-version-2.3/JAFFAL.groovy \
test.fasta
```

### Expected Output
* `test.csv` — predicted fusion transcript table
* `test.fasta` — fusion sequences
* `test/` — JAFFAL intermediate files

---

## D. Identify Putative Fusion Transcripts from Real Data
### Purpose
Run the optimized JAFFAL pipeline on real long-read sequencing datasets.

### Input
* Example FASTA file: `test_data/Isoseq.fasta`

### Workflow
#### Step D1. Run JAFFAL
```bash
$JAFFAL_HOME/tools/bin/bpipe run -n 10 \
-p refBase=~/PATH_to_JAFFAL_DB/JAFFAL_DB_MH63RS3 \
-p genome=Omh63 \
-p annotation=RS3 \
-p outputName=Isoseq \
~/PATH_to_JAFFAL_DB/JAFFA-version-2.3/JAFFAL.groovy \
Isoseq.fasta
```

### Expected Output
* `Isoseq.csv`
* `Isoseq.fasta`
* `Isoseq/`

---

## E. Reformat JAFFAL Output for Downstream Analysis
### Purpose
Convert raw JAFFAL output into formats suitable for statistical analyses, visualization, and genome browser inspection.

### Workflow
#### Step E1. Generate ranked fusion candidate table
```bash
sed '1d' Isoseq.csv | \
perl scripts/reformat_jaffal.pl - | \
csvtk -tH sort -k 5:u -L 5:ranks > mh63_flagleaf.pc
```

#### Step E2. Generate BEDPE file for breakpoint visualization
```bash
sed '1d' Isoseq.csv | \
perl scripts/convert_jaffal_to_bedpe.pl - | \
sortBed -i - > mh63_flagleaf.pc.bedpe
```

#### Step E3. Assign unique IDs to each fusion candidate
```bash
awk '{print $0"\t"$10"|"$1":"$2"|"$11}' \
mh63_flagleaf.pc > mh63_flagleaf.pc.uniqueID
```

### Output
* `mh63_flagleaf.pc`
* `mh63_flagleaf.pc.bedpe`
* `mh63_flagleaf.pc.uniqueID`

### Notes
Example output files are available in the `examples/` directory.



# 03. Visualization of Fusion Transcripts and Fusion Breakpoints
## Purpose
Generate IGV-compatible files for visual inspection of putative fusion transcripts and breakpoint structures. This workflow enables visualization of parental genes, gene boundaries, annotations, SNP positions, protein domains, transposable elements, and alignments of fusion-supporting reads.

## Related Figures
- Supplementary Figures S6a–b, S7, S9, S12, S13b, S15, S16, S19b–c, S21b, S22
- Figure 2c

---

## A. IGV Visualization of Putative Fusion Transcripts
### Purpose
Construct custom pseudo-fusion reference sequences and annotation tracks for direct visualization of candidate fusion transcripts in IGV.

### Input Files
The following reference files can be generated from the MH63 genome/annotation or are available in the `examples/` directory:

- `Omh63_RS3.txt` — genomic coordinates of MH63 genes  
- `Omh63.transSNPs.bed` — SNP positions in transcript coordinates  
- `Omh63.geneLocus.fa.length` — gene length table  
- `Omh63.transcriptAnno.gtf` — transcript-coordinate annotation file  
- `Omh63_RS3.fasta` — gene locus FASTA file (available from the JAFFAL reference database)

### Workflow
#### Step A1. Convert candidate fusion list into visualization input format

```bash
cut -f 1-4 mh63_flagleaf.pc | \
sed "s/:/\t/g" | \
cut -f 1-3,5 | \
sed "s/\t/:/" | \
sed "s/\t/#/g" | \
sort | uniq > mh63_flagleaf.pc.visual.input
````

Example output:
```text
OsMH_01G0005000:OsMH_01G0005100#Chr01#Chr01
OsMH_01G0006900:OsMH_08G0138000#Chr01#Chr08
```

#### Step A2. Generate pseudo-fusion reference files
```bash
cat mh63_flagleaf.pc.visual.input | \
while read line; do
sh scripts/03_fusion_visualization/create_fake_fusionTrans.refMH63RS3.sh \
$line mh63_flagleaf.pc
done
```

### Output Files
* `mh63_flagleaf.pc.visual.fa` — pseudo-fusion reference FASTA
* `mh63_flagleaf.pc.visual.gtf` — fusion annotation file
* `mh63_flagleaf.pc.visual.gene_boundary.bed` — gene boundary track
* `mh63_flagleaf.pc.visual.SNPs.bed` — SNP track
* `mh63_flagleaf.pc.visual.breakpoint.bed` — predicted breakpoint track

---

#### Step A3. Align JAFFAL fusion reads to pseudo-fusion references
```bash
cp ~/PATH_to_JAFFAL_results/Isoseq.fusions.fa \
mh63_flagleaf.jaffal.fusions.fa

minimap2 -ax splice:hq -uf \
mh63_flagleaf.pc.visual.fa \
mh63_flagleaf.jaffal.fusions.fa | \
samtools view -b | \
samtools sort > mh63_flagleaf.jaffal.fusions.bam

samtools index mh63_flagleaf.jaffal.fusions.bam
```

### Output Files
* `mh63_flagleaf.jaffal.fusions.bam`
* `mh63_flagleaf.jaffal.fusions.bam.bai`

### Example Figure
* `figures/fake_fusion_IGVexample.png`

---

## B. Visualization of Fusion Breakpoint Types
### Purpose
Visualize local breakpoint structures by aligning parental gene sequences to fusion transcript sequences. This enables classification of breakpoint types such as SHS-mediated overlaps.

### Related Figure
* Figure 2c

### Input Files
* `Isoseq.fasta`
* `mh63_flagleaf.pc`

### Workflow
#### Step B1. Reformat fusion transcript FASTA and IDs
```bash
awk '{print $10"|"$1":"$2"|"$11}' \
mh63_flagleaf.pc > mh63_flagleaf.pc.reads.id

cat Isoseq.fasta | \
fasta_formatter -w 80 > mh63_flagleaf.pc.reads.fa
```

#### Step B2. Extract parental gene sequences
```bash
cut -f 1,2 mh63_flagleaf.pc | \
sed "s/\t/\n/" | \
sort | uniq > mh63_flagleaf.pc.reads.geneID

seqtk subseq Omh63.geneLocus.fa \
mh63_flagleaf.pc.reads.geneID | \
fasta_formatter -w 80 > mh63_flagleaf.pc.reads.gene.fa
```

#### Step B3. Align parental genes to fusion transcripts
```bash
blastn \
-subject mh63_flagleaf.pc.reads.fa \
-query mh63_flagleaf.pc.reads.gene.fa \
-outfmt '6 qseqid qlen qstart qend sstrand sseqid slen sstart send nident length bitscore' \
-perc_identity 95 \
> mh63_flagleaf.pc.reads.gene.blastn
```

#### Step B4. Convert BLAST results into IGV BED tracks
```bash
awk '{
if($5=="minus"){s=$9-1;e=$8}
else{s=$8-1;e=$9};
print $6"\t"s"\t"e"\t"$1
}' mh63_flagleaf.pc.reads.gene.blastn | \
sortBed > mh63_flagleaf.pc.reads.gene.bed
```

### Output Files
* `mh63_flagleaf.pc.reads.fa`
* `mh63_flagleaf.pc.reads.gene.bed`

### Example Figure
* `figures/Breakpoints_SHS_example.png`


# 04. Classification of Fusion Breakpoint Types
## Purpose
Classify the breakpoints of putative fusion transcripts identified by JAFFAL into four categories: **short homologous sequence (SHS)**, **Joined**, **Unknown**, and **Misaligned**. This analysis helps infer the likely origins of false-positive fusion transcripts.

## Related Figures
- Figure 2a–b
- Supplementary Figures S4c, S5a

---

## A. Convert Breakpoints from Genomic Coordinates to Transcript Coordinates
### Purpose
Transform fusion breakpoint coordinates from genomic space into transcript space for downstream breakpoint classification.

### Input
- `mh63_flagleaf.pc`
- `Isoseq.csv`
- `Omh63.geneLocus.bed`
- `Omh63.gtf`

### Workflow
#### Step A1. Extract breakpoint coordinates in genomic space

```bash
cut -f 1,2,3,4,6,7 mh63_flagleaf.pc | \
sed "s/:/\t/g" | \
perl scripts/04_breakpoint_classification/get_breakpoints_genomic_coords.pl - | \
sortBed -i - > mh63_flagleaf.pc.breakpoints.genomic.bed
````

#### Step A2. Map genomic breakpoints to transcript coordinates in R
```r
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(readr)

geneLoc <- import.bed("Omh63.geneLocus.bed")
names(geneLoc) <- geneLoc$name

genomeGtf <- import.gff("Omh63.gtf", format = "gtf")
fusionMap <- import.bed("mh63_flagleaf.pc.breakpoints.genomic.bed")

fusionMap_trans <- mapToTranscripts(x = fusionMap, transcripts = geneLoc)

fusionMap %>%
  as.data.frame() %>%
  dplyr::select(seqnames, start, end, name) %>%
  add_column(order = 1:length(fusionMap)) %>%
  merge(
    x = fusionMap_trans %>% as.data.frame() %>% dplyr::select(seqnames, start, end, strand, xHits),
    y = .,
    by.x = "xHits",
    by.y = "order"
  ) %>%
  write_tsv(file = "mh63_flagleaf.pc.breakpoints.trans.txt")
```

#### Step A3. Convert transcript-space breakpoints to BED format
```bash
sed '1d' mh63_flagleaf.pc.breakpoints.trans.txt | \
awk '{print $2"\t"$3-1"\t"$4"\t"$9}' | \
sortBed > mh63_flagleaf.pc.breakpoints.trans.bed
```

### Output
* `mh63_flagleaf.pc.breakpoints.genomic.bed`
* `mh63_flagleaf.pc.breakpoints.trans.txt`
* `mh63_flagleaf.pc.breakpoints.trans.bed`

---

## B. Prepare Breakpoint and Coverage Files on Fusion Reads
### Purpose
Generate breakpoint and parental-gene coverage tracks on fusion transcript sequences.

### Input
* `mh63_flagleaf.pc`
* `Isoseq.csv`
* `mh63_flagleaf.pc.reads.gene.bed`
* `mh63_flagleaf.pc.reads.fa`

### Workflow
#### Step B1. Extract breakpoint coordinates on fusion reads
```bash
perl scripts/04_breakpoint_classification/get_breakpoints_on_read.pl \
mh63_flagleaf.pc Isoseq.csv > mh63_flagleaf.pc.reads.breakpoints.bed
```

#### Step B2. Merge parental-gene alignment intervals
```bash
mergeBed -d -30 -c 4 -o distinct \
-i mh63_flagleaf.pc.reads.gene.bed > mh63_flagleaf.pc.reads.gene.merged.bed
```

#### Step B3. Build FASTA index and length table
```bash
samtools faidx mh63_flagleaf.pc.reads.fa
cut -f 1,2 mh63_flagleaf.pc.reads.fa.fai > mh63_flagleaf.pc.reads.fa.length
```

### Output
* `mh63_flagleaf.pc.reads.breakpoints.bed`
* `mh63_flagleaf.pc.reads.gene.merged.bed`
* `mh63_flagleaf.pc.reads.fa.length`

---

## C. Classify Breakpoint Types Based on Parental-Gene Coverage
### Purpose
Classify fusion breakpoints according to how the aligned parental genes cover the breakpoint region on the fusion read. SHS events are further separated by overlap length.

### Workflow
#### Step C1. Compute parental-gene coverage across fusion reads
```bash
genomeCoverageBed \
-i mh63_flagleaf.pc.reads.gene.merged.bed \
-g mh63_flagleaf.pc.reads.fa.length \
-bga > mh63_flagleaf.pc.reads.gene.coverage
```

#### Step C2. Identify SHS breakpoints with overlap >3 bp
```bash
intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
awk '$7==2{print $1"\t"$5"\t"$6"\t"$6-$5}' | \
awk '$4>3' > mh63_flagleaf.pc.reads.out.shsGreat3.bed
```

#### Step C3. Generate temporary support files
```bash
cut -f 1 mh63_flagleaf.pc.reads.out.shsGreat3.bed > mh63_flagleaf.pc.reads.out.shsGreat3.id

intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
awk '$7==2{print $1"\t"$5"\t"$6"\t"$6-$5}' \
> mh63_flagleaf.pc.reads.out.interBed

awk '$4<=3' mh63_flagleaf.pc.reads.out.interBed > mh63_flagleaf.pc.reads.out.join1.bed

cut -f 1 mh63_flagleaf.pc.reads.out.interBed > tmp.id

awk '$4==0' mh63_flagleaf.pc.reads.gene.coverage > mh63_flagleaf.pc.reads.out.coverage0
awk '$4==2' mh63_flagleaf.pc.reads.gene.coverage > mh63_flagleaf.pc.reads.out.coverage2
```

#### Step C4. Identify Joined breakpoints
```bash
intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
awk '$7<2' | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4=="."' | \
cut -f 1-3 > mh63_flagleaf.pc.reads.out.join2.bed
```

```bash
intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
awk '$7<2' | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4!="."' | \
awk '{print $1"\t"$2"\t"$3"\t"$6-$5"\t"$8}' | \
awk '$5>=30' | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage2 -d | \
awk '{print $1"\t"$2"\t"$3"\t"$8-$7"\t"$10}' | \
awk '$4<=3 || $5>10' > mh63_flagleaf.pc.reads.out.join3.bed
```

```bash
intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4!="."' | \
awk '{print $1"\t"$2"\t"$3"\t"$6-$5"\t"$8}' | \
awk '$5<30 && $4<30' | \
awk '$4==1 || $4==2' > mh63_flagleaf.pc.reads.out.join4.bed
```

#### Step C5. Identify additional SHS, Unknown, and Misaligned breakpoints
```bash
intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4!="."' | \
awk '{print $1"\t"$2"\t"$3"\t"$6-$5"\t"$8}' | \
awk '$5>=30' | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage2 -d | \
awk '{print $1"\t"$7"\t"$8"\t"$8-$7"\t"$10}' | \
awk '$4>3 && $5<=10' > mh63_flagleaf.pc.reads.out.shsGreat3_1.bed
```

```bash
intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4!="."' | \
awk '{print $1"\t"$2"\t"$3"\t"$6-$5"\t"$8}' | \
awk '$5<30 && $4<30' | \
awk '$4!=1 && $4!=2' > mh63_flagleaf.pc.reads.out.unknown.bed
```

```bash
intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4!="."' | \
awk '{print $1"\t"$2"\t"$3"\t"$6-$5"\t"$8}' | \
awk '$4>=30 && $5<30' > mh63_flagleaf.pc.reads.out.error.bed
```

---

## D. Merge Final Classification Results
### Purpose
Combine all intermediate outputs into the four final breakpoint-type files.

### Workflow
#### Step D1. Merge SHS breakpoints
```bash
cat mh63_flagleaf.pc.reads.out.shsGreat3.bed \
mh63_flagleaf.pc.reads.out.shsGreat3_1.bed | \
cut -f 1-4 > mh63_flagleaf.pc.reads.out.shs.bed
```

#### Step D2. Merge Joined breakpoints
```bash
cat mh63_flagleaf.pc.reads.out.join?.bed | \
cut -f 1-3 | \
awk '{print $0"\t0"}' > mh63_flagleaf.pc.reads.out.join.bed
```

### Final Output
* `mh63_flagleaf.pc.reads.out.shs.bed` — SHS-type breakpoints
* `mh63_flagleaf.pc.reads.out.join.bed` — Joined-type breakpoints
* `mh63_flagleaf.pc.reads.out.unknown.bed` — Unknown-type breakpoints
* `mh63_flagleaf.pc.reads.out.error.bed` — Misaligned-type breakpoints

### Notes
These example output files are available in the `test_data/` directory.



# 05. Fusion Network Analysis
## Purpose
Construct gene-level fusion networks from putative fusion transcript pairs and identify network modules for downstream topological analyses.

## Related Figures
- Figure 4a–e

## Script
- `scripts/05_fusion_network/network.R`

## Input
- `mh63_flagleaf.pc` — list of putative fusion transcripts (available in the `examples/` directory)

## Software / R package
- `igraph`

## Description
Fusion network analysis was performed using the R script `network.R`. In this analysis, genes involved in putative fusion transcripts were represented as nodes, and fusion relationships between gene pairs were represented as edges. The resulting network was further partitioned into modules for downstream analysis of network structure and gene connectivity.

## Output
- `isoseq_gene_modules.txt` — gene-to-module assignment table
- `isoseq_gene_modules_stats.txt` — module summary statistics, including module ID, number of genes per module, and average degree


# 06. Fusion Transcript Identification from Short-Read RNA-seq Data Using Arriba
## Purpose
Identify putative fusion transcripts from short-read RNA-seq data using the Arriba pipeline. These results were used for analyses of adjacent-gene fusions and for comparison with long-read RNA sequencing datasets.

## Related Figures
- Figure 5a
- Supplementary Figures S11, S13a

## Description
Fusion transcripts were identified from paired-end short-read RNA-seq data using STAR + Arriba. Because very small test datasets may not yield sufficient chimeric signals, we recommend using a complete RNA-seq dataset.

The short-read datasets used in this study are publicly available in NCBI under accession:
- `PRJNA1291274`

---

## A. Build STAR Genome Index
### Purpose
Generate the STAR genome index using the MH63 reference genome and annotation.

### Input
- `Omh63.fa`
- `Omh63.gtf`

### Workflow
#### Step A1. Build STAR index

```bash
STAR \
--runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles Omh63.fa \
--sjdbGTFfile Omh63.gtf \
--runThreadN 4 \
--genomeSAindexNbases 13
````

---

## B. Align Short Reads to the MH63 Reference Genome
### Purpose
Align paired-end RNA-seq reads using STAR with chimeric alignment settings required for Arriba fusion calling.

### Input
* `data1.fq.gz`
* `data2.fq.gz`

### Workflow
#### Step B1. Run STAR alignment
```bash
STAR \
--runMode alignReads \
--runThreadN 5 \
--genomeDir mh63 \
--genomeLoad NoSharedMemory \
--readFilesIn data1.fq.gz data2.fq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outFilterMultimapNmax 50 \
--peOverlapNbasesMin 10 \
--alignSplicedMateMapLminOverLmate 0.5 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--chimSegmentMin 10 \
--chimOutType WithinBAM HardClip \
--chimJunctionOverhangMin 10 \
--chimScoreDropMax 30 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreSeparation 1 \
--chimSegmentReadGapMax 3 \
--chimMultimapNmax 50 \
--outFileNamePrefix output_pre
```

### Output
* `output_preAligned.sortedByCoord.out.bam`

---

## C. Call Fusion Transcripts Using Arriba
### Purpose
Detect candidate fusion transcripts from STAR-aligned BAM files.

### Workflow
#### Step C1. Run Arriba
```bash
arriba \
-x output_preAligned.sortedByCoord.out.bam \
-g Omh63.gtf \
-a Omh63.fa \
-o output_pre.arribaFusion.tsv \
-O output_pre.arribaFusion.discarded \
-f blacklist \
-i Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12
```

### Output
* `output_pre.arribaFusion.tsv` — candidate fusion transcript table

---

## D. Merge and Filter Fusion Calls Across Samples
### Purpose
Reformat individual Arriba outputs, merge biological replicates, and generate a final filtered fusion set for downstream analyses.

### Workflow
#### Step D1. Reformat each Arriba result file
```bash
for i in *.tsv; do
pre=$(echo $i | cut -f 1 -d '.')
perl reformat_tsv.pl $i | \
awk -v n=$pre '{print $0"\t"n}' > $pre.basic
done
```

#### Step D2. Merge replicate samples
```bash
cat *.basic | perl merge_rep.pl - > flagleaf.arriba
```

#### Step D3. Filter repetitive or low-confidence calls
```bash
for i in *.arriba; do
pre=$(basename $i | cut -f 1 -d '.')
cat $i | perl filter_repeats.pl - > $pre.arriba
done
```

#### Step D4. Merge all flag leaf treatment groups
```bash
cat flagleaf_HTLD.arriba \
flagleaf_HTSD.arriba \
flagleaf_LTLD.arriba \
flagleaf_LTSD.arriba | \
perl merge_flagleaf.pl - | \
csvtk -tH sort -k 13:u -L 13:ranks > flagleaf.result.filter1_rep
```

### Final Output
* `flagleaf.result.filter1_rep`

### Example Output
```text
OsMH_06G0046600 OsMH_06G0046500 Chr06:2654636 Chr06:2643216 -/- -/- intra N/S CDS 5'UTR deletion/read-through out-of-frame high 36 5 279

OsMH_05G0422200 OsMH_09G0065600 Chr05:28230500 Chr09:4496942 -/- +/+ inter S/S CDS CDS translocation out-of-frame high 109 8318
```



# 07. Identification of Allelic and Non-Allelic Fusion Transcripts in Hybrid Rice Using the Trans-Mate Pair Strategy
## Purpose
Identify candidate **allelic** and **non-allelic** fusion transcripts in the hybrid rice line SY63 using a trans-mate pair strategy based on paired-end short-read RNA-seq data.

## Description
This workflow uses paired-end RNA-seq reads from the hybrid SY63 and mixed-parent controls. Read pairs are classified according to their parental mapping patterns against the MH63 and ZS97 genomes. Candidate fusion-supporting read pairs are then used to infer possible allelic or non-allelic fusion genes.

## Recommended Input
Complete raw paired-end RNA-seq FASTQ datasets are recommended for optimal sensitivity.

## Related Figures
- Supplementary Figure S18a–g

---

## A. Phase Reads and Assess Data Quality
### Purpose
Phase RNA-seq reads to MH63 or ZS97 haplotypes and evaluate alignment quality.

### Input
- `toMH63.fs.bam`
- `toZS97.fs.bam`
- `MH63vsZS97.snp`

### Workflow
#### Step A1. Phase reads using PP2PG
```bash
phasing.py \
--pt=RNAseq \
--g1=toMH63.fs.bam \
--g2=toZS97.fs.bam \
--snp=MH63vsZS97.snp \
--gop1=MH63 \
--gop2=ZS97 \
--st=5
````

#### Step A2. Sort phased BAM files
```bash
samtools sort -@ 10 mh.bam > mh.sort.bam
samtools index mh.sort.bam

samtools sort -@ 10 zs.bam > zs.sort.bam
samtools index zs.sort.bam
```

#### Step A3. Evaluate BAM correlations
```bash
multiBamSummary bins \
-b *sort.bam \
-o rna.multiBamSummary \
-bs 100000 \
--outRawCounts rna.multiBamSummary.rawCounts \
--scalingFactors rna.multiBamSummary.scalingFactors \
-p 30
```

#### Step A4. Generate coverage tracks
```bash
bamCoverage \
-b mh.sort.bam \
-o mh.sort.bw \
-of bigwig \
-bs 1 \
-p 10
```

---

## B. Extract Read Pairs with Differential Parental Matching Patterns
### Purpose
Identify paired reads in which one mate maps perfectly while the other contains mismatches, enabling parental-origin classification.

### Workflow
#### Step B1. Name-sort BAM files
```bash
samtools view -@ 4 -b -F 256 toMH63.bam | \
samtools sort -@ 4 -n - > toMH63.nameSorted.bam

samtools view -@ 4 -b -F 256 toZS97.bam | \
samtools sort -@ 4 -n - > toZS97.nameSorted.bam
```

#### Step B2. Convert BAM to BEDPE
```bash
bamToBed -bedpe -ed -mate1 -i toMH63.nameSorted.bam > toMH63.nameSorted.bedpe

bamToBed -bedpe -ed -mate1 -i toZS97.nameSorted.bam > toZS97.nameSorted.bedpe
```

#### Step B3. Extract mismatch statistics (NM tag)
```bash
samtools view toMH63.nameSorted.bam | \
perl scripts/parse_bam_for_NM.pl - > toMH63.nameSorted.stats

samtools view toZS97.nameSorted.bam | \
perl scripts/parse_bam_for_NM.pl - > toZS97.nameSorted.stats
```

---

## C. Merge Parental Alignment Statistics
### Purpose
Combine mismatch information from both parental alignments to determine read-pair matching categories.

### Workflow
#### Step C1. Merge NM-tag summaries
```bash
paste \
<(cut -f 1,4,8 toMH63.nameSorted.stats) \
<(cut -f 4,8 toZS97.nameSorted.stats) \
> sy63.tagNM
```

#### Step C2. Convert NM tags into match categories
```bash
cat sy63.tagNM | \
perl nm_to_match.pl - > sy63.parentalAlign
```

Example:
```text
E250031905L1C001R00200000011 Y Y Y Y
E250031905L1C001R00200000049 Y Y Y N
```

#### Step C3. Summarize category proportions
```bash
cat sy63.parentalAlign | perl make_category.pl -
```

---

## D. Extract Cis- and Trans-Mate Pair Reads
### Purpose
Label each read pair as cis- or trans-mate and retain informative candidate pairs.

### Workflow
#### Step D1. Assign read-pair categories
```bash
cat sy63.parentalAlign | \
perl get_category.pl - > sy63.cistrans
```

Example:
```text
E250031905L1C001R00200000974 MHcis
E250031905L1C001R00200001211 ZScis
```

#### Step D2. Extract selected BEDPE pairs
```bash
perl get_bedpe.pl sy63.cistrans toMH63.nameSorted.bedpe > sy63.mh.bedpe.filter
perl get_bedpe.pl sy63.cistrans toZS97.nameSorted.bedpe > sy63.zs.bedpe.filter
```

#### Step D3. Convert BEDPE to BED format
```bash
perl bed12Tobed6_mh.pl sy63.cistrans sy63.mh.bedpe.filter > sy63.mh.bedpe.filter.bed
perl bed12Tobed6_zs.pl sy63.cistrans sy63.zs.bedpe.filter > sy63.zs.bedpe.filter.bed
```

---

## E. Identify Candidate Allelic and Non-Allelic Fusion Genes
### Purpose
Intersect cis/trans mate pairs with exon annotations and infer candidate hybrid fusion genes.
### Related Figures
* Supplementary Figure S18b–g

### Workflow
#### Step E1. Retain longest isoforms
```bash
agat_sp_keep_longest_isoform.pl \
-f Omh63.maker.gff \
-o Omh63.longest.gff

agat_sp_keep_longest_isoform.pl \
-f Ozs97.maker.gff \
-o Ozs97.longest.gff
```

#### Step E2. Add exon order
```bash
perl add_exon_order.pl Omh63.longest.gff | \
sort -k1,1 -k4,4n -k5,5 > Omh63.exon.gff

perl add_exon_order.pl Ozs97.longest.gff | \
sort -k1,1 -k4,4n -k5,5 > Ozs97.exon.gff
```

#### Step E3. Intersect read pairs with exons
```bash
intersectBed \
-a sy63.mh.bedpe.filter.bed \
-b Omh63.exon.gff \
-wao > sy63.mh.bedpe.filter.bed.exon

intersectBed \
-a sy63.zs.bedpe.filter.bed \
-b Ozs97.exon.gff \
-wao > sy63.zs.bedpe.filter.bed.exon
```

#### Step E4. Determine overlapping genes
```bash
perl reads_overlap_exon.pl \
Omh63.gene2trans.map \
sy63.mh.bedpe.filter.bed.exon > sy63.mh.readStatus

perl reads_overlap_exon.pl \
Ozs97.gene2trans.map \
sy63.zs.bedpe.filter.bed.exon > sy63.zs.readStatus
```

#### Step E5. Identify allelic fusion candidates
```bash
cat sy63.mh.readStatus sy63.zs.readStatus | \
perl cis_category.pl - | \
grep 'Diff_exon' > sy63.merged.cis.diffexon

cat sy63.merged.cis.diffexon | \
perl cis_map_gene.pl - | \
sort -k1,1 -k2,2 > sy63.merged.cis.diffexon.exonLevel

cat sy63.merged.cis.diffexon.exonLevel | \
perl convert_cis_exon_to_gene.pl - | \
sort -k1,1 > sy63.merged.cis.diffexon.geneLevel
```

#### Step E6. Identify non-allelic fusion candidates
```bash
cat sy63.mh.readStatus sy63.zs.readStatus | \
awk '$2~"trans"' | \
perl merge_trans_pair.pl - > sy63.merged.trans

cat sy63.merged.trans | \
perl trans_categroy.pl - | \
awk '$1=="Pass"' | \
cs-th cut -f -1 | \
perl trans_find_properPair.pl - > sy63.trans.diffexon

cat sy63.trans.diffexon | \
grep -v 'multie' | \
perl trans_map_gene.pl - > sy63.trans.diffexon.exonLevel

cat sy63.trans.diffexon.exonLevel | \
perl convert_trans_exon_to_gene.pl - | \
sort -k1,1 > sy63.trans.diffexon.geneLevel
```

### Final Output
* `sy63.merged.cis.diffexon.geneLevel` — candidate allelic fusion genes
* `sy63.trans.diffexon.geneLevel` — candidate non-allelic fusion genes



# 08. Additional Analyses
## Purpose
This section includes three additional analyses used in the study:
1. overlap between fusion transcripts and 3D interaction data,
2. primary genome assembly of MH63 using HiFi reads,
3. comparison of different genome annotation versions.

---

## A. Overlap Between Fusion Transcripts and 3D Interaction Data
### Purpose
Evaluate whether putative fusion transcripts overlap with chromatin or RNA interaction signals, including DNA–DNA, DNA–RNA, and RNA–RNA interaction datasets.

### Related Figure
- Figure 1f

### Data Sources
The rice 3D interaction datasets used in this study were obtained from the following publications:

#### DNA–DNA interactions
- **H3K4me3 ChIA-PET**  
  https://www.nature.com/articles/s41467-019-11535-9  
  Accession: `GSM3767544`
- **RNAPII ChIA-PET**  
  https://www.nature.com/articles/s41467-019-11535-9  
  Accession: `GSM3767545`
 
#### DNA–RNA interactions
- **ChRD–PET**  
  https://www.nature.com/articles/s41477-021-01089-4#data-availability  
  Accession: `GSE163119`
- **TaDRIM-seq**  
  https://www.nature.com/articles/s41467-024-53534-5  
  Accession: `GSE234228`
#### RNA–RNA interactions
- **TaDRIM-seq**  
  https://www.nature.com/articles/s41467-024-53534-5  
  Accession: `GSE234228`

### Example Test Data
- `test_data/dna2rna_H3K4me3_rs3.bedpe`
- `test_data/isoseq.fusions.promoter.bedpe`

### Workflow
#### Step A1. Intersect fusion transcript pairs with 3D interaction data
```bash
pairToPair \
-a test_data/isoseq.fusions.promoter.bedpe \
-b test_data/dna2rna_H3K4me3_rs3.bedpe \
-is > isoseq.dna2rna_H3K4me3.observed
````

#### Step A2. Generate matched random gene pairs and compute overlaps
```bash
sh random_gene_pair.sh

pairToPair \
-a test_data/random_expressed_pair.promoter.bedpe \
-b test_data/dna2rna_H3K4me3_rs3.bedpe \
-is > dna2rna_H3K4me3.random
```

### Output
* `isoseq.dna2rna_H3K4me3.observed`
* `dna2rna_H3K4me3.random`

---

## B. Assembly of the MH63 Primary Genome Using HiFi Reads
### Purpose
Assemble the MH63 primary genome from PacBio HiFi sequencing data.

### Related Figure
* Supplementary Figure S6b

### Data Source
The MH63 HiFi sequencing data used for assembly are available from NCBI:
* [https://www.ncbi.nlm.nih.gov/sra/SRR10238608](https://www.ncbi.nlm.nih.gov/sra/SRR10238608)

### Workflow
#### Step B1. Run hifiasm

```bash
hifiasm -o MH63.asm -t 32 SRR10238608.fastq.gz
```

### Output
* `MH63.asm.*`

---

## C. Comparison of Different Genome Annotation Versions
### Purpose
Compare different annotation versions of the same genome at the gene and transcript levels, with a particular focus on:
* **split-type annotations**, in which one gene in one annotation corresponds to two or more genes in the other annotation,
* **fused-type annotations**, in which two or more genes in one annotation are merged into one gene in the other annotation.

### Related Figure
* Figure 13c–e

### Example Comparison
The example below compares `MH63RS3` and `MH63RS3_man`.

### Input
* `Omh63.gff3`
* `Omh63_man.gff3`
* `Omh63.genomic.fa`

### Additional Annotation Source
* `MH63RS3_man` can be downloaded from:
  [http://rice.hzau.edu.cn/rice_rs3/download_ext/mh63-2025_annotation.gff3.gz](http://rice.hzau.edu.cn/rice_rs3/download_ext/mh63-2025_annotation.gff3.gz)

### Workflow
#### Step C1. Compare overall annotation statistics
```bash
agat_sp_statistics.pl --gff Omh63.gff3 -g Omh63.genomic.fa -o mh63.stats
agat_sp_statistics.pl --gff Omh63_man.gff3 -g Omh63.genomic.fa -o mh63_man.stats
```

#### Step C2. Convert exon annotations to BED format
```bash
agat_convert_sp_gff2bed.pl --gff Omh63.gff3 --sub exon -o Omh63.bed
agat_convert_sp_gff2bed.pl --gff Omh63_man.gff3 --sub exon -o Omh63_man.bed
```

#### Step C3. Compare annotation structures using ParsEval
```bash
parseval -d -o parseval.out Omh63.bed Omh63_man.bed &> parseval.log
```

#### Step C4. Identify split-type gene annotations
```bash
intersectBed \
-a <(awk '$3=="gene"' Omh63.gff3) \
-b <(awk '$3=="gene"' Omh63_man.gff3) \
-wao -s | \
sed 's/ID=//g' | \
awk '$18!="."' | \
csvtk_tH fold -f 9 -v 18 -s ";" | \
awk '{n=sub(/;/,";",$2); split($2,a,";"); if(n>=1){print $0"\t"length(a)}}' \
> v1_split_in_v2.id
```

#### Step C5. Identify fused-type gene annotations
```bash
intersectBed \
-a <(awk '$3=="gene"' Omh63_man.gff3) \
-b <(awk '$3=="gene"' Omh63.gff3) \
-wao -s | \
sed 's/ID=//g' | \
awk '$18!="."' | \
csvtk_tH fold -f 9 -v 18 -s ";" | \
awk '{n=sub(/;/,";",$2); split($2,a,";"); if(n>=1){print $0"\t"length(a)}}' \
> v1_fuse_in_v2.id
```

#### Step C6. Summarize new transcripts and detailed locus comparisons
```bash
perl resolve_parseval_out.newTranscripts.pl parseval.out | \
awk '$NF!=0' > parseval.out.newTrans

perl resolve_parseval_out.pl parseval.out > parseval.out.locusNum

perl resolve_parseval_out.comparison.pl parseval.out > parseval.out.comparison
```

### Output

* `mh63.stats`
* `mh63_man.stats`
* `parseval.out`
* `parseval.log`
* `v1_split_in_v2.id`
* `v1_fuse_in_v2.id`
* `parseval.out.newTrans`
* `parseval.out.locusNum`
* `parseval.out.comparison`