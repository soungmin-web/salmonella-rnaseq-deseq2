#!/usr/bin/env bash
# From raw reads to BAM files: download, QC, trimming, QC after trimming, mapping
#
# Usage: bash preprocessing.sh   (Linux, run in an empty folder)
# Tools: SRA Toolkit, FastQC, MultiQC, Trimmomatic 0.39, STAR 2.7.10a, samtools
set -euo pipefail

TRIMMOMATIC_JAR=/usr/share/java/trimmomatic-0.39.jar
ADAPTERS=/usr/share/java/Trimmomatic-0.39/adapters/TruSeq3-PE.fa

GLYCEROL="SRR10320667 SRR10320668 SRR10320669"
SUCCINATE="SRR10320670 SRR10320671 SRR10320672"
SAMPLES="$GLYCEROL $SUCCINATE"

# 1. Download the paired-end reads (first use: run "vdb-config --interactive" once)
for s in $SAMPLES; do
  fastq-dump --split-files "$s"
done

# 2. Quality control of the raw reads
fastqc *.fastq
multiqc .

# 3. Trim adapters and low-quality bases
for s in $SAMPLES; do
  java -jar "$TRIMMOMATIC_JAR" PE -threads 2 -phred33 \
    "${s}_1.fastq" "${s}_2.fastq" -baseout "${s}_trimmed.fastq" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 HEADCROP:1
done

# 4. Quality control of the trimmed reads (paired output only)
fastqc *_trimmed_*P.fastq
multiqc ./*trimmed*

# 5. Reference genome: S. Typhimurium str. 14028S
curl -OJX GET \
  "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000022165.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF&filename=GCF_000022165.1.zip" \
  -H "Accept: application/zip"
unzip GCF_000022165.1.zip
mv ncbi_dataset/data/GCF_000022165.1/* .

# 6. Index the genome
mkdir genome_dir
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir genome_dir \
  --genomeFastaFiles GCF_000022165.1_ASM2216v1_genomic.fna \
  --sjdbGTFfile genomic.gff --genomeSAindexNbases 10 \
  --sjdbGTFtagExonParentTranscript Parent

# 7. Map the trimmed reads (alignIntronMax 1: bacterial genes have no introns)
for s in $SAMPLES; do
  STAR --runThreadN 4 \
    --readFilesIn "${s}_trimmed_1P.fastq" "${s}_trimmed_2P.fastq" \
    --genomeDir genome_dir --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "$s" --outSAMunmapped Within --alignIntronMax 1
  cat "${s}Log.final.out"
done

# 8. Index the BAM files for IGV
for s in $SAMPLES; do
  samtools index "${s}Aligned.sortedByCoord.out.bam"
done
