#!/bin/bash
cp samtools ~/bin

cd ~
mkdir samtools-03713
cd samtools-03713

#download the gzipped sample SAM obtained from Bowtie

curl https://_____/sample.sam.gz > sample.sam.gz
gzip -d sample.sam.gz

# convert SAM to BAM and view (optional):
samtools view -S -b sample.sam > sample.bam
samtools view sample.bam | head


#R-SEM
echo "\n\n----------------R-SEM starts---------------\n\n"

#unzip the RNA-seq
unzip -d data data/SRR_____.zip

#quantize expression compared to req-genome
#options to pick between buid ref-seq or use prebuilt ref
software/RSEM-1.2.25/rsem-calculate-expression -p 8 --paired-end \
					--bowtie2 --bowtie2-path software/bowtie2-2.2.6 \
					--estimate-rspd \
					--append-names \
					--output-genome-bam \
					data/SRR937564_1.fastq data/SRR937564_2.fastq\
					ref/mouse_ref exp/LPS_6h


software/RSEM-1.2.25/rsem-calculate-expression -p 8 --paired-end \
					--bam \
					--estimate-rspd \
					--append-names \
					--output-genome-bam \
					exp/LPS_6h.bam \
					ref/human_ref exp/LPS_6h

