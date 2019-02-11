#!/bin/bash
mkdir project/rsap_files
gunzip project/data/genome/Homo_sapiens.GRCh38.95.gtf.gz
R-SAP --in1 $1 --in2 project/data/genome/Homo_sapiens.GRCh38.95.gtf --outDir project/rsap_files --rf GTF --tNum 4
rm project/data/genome/Homo_sapiens.GRCh38.95.gtf

#1 Reference genome alignment file (sorted on RNA‐Seq IDs) of the RNA‐Seq reads in psl format.
#Reference transcript annotation file (UCSC table browser, BED or GTF format). 
#We are looking for rsap_files/ChimericTranscriptAnnotation.out
