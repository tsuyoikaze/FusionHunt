#!/bin/bash
echo "Running: mkdir $3/rsap_files"
mkdir "$3/rsap_files"

FLAG=false
FILENAME=$2
EMPTY_STR=""

if [ ${file: -3} == ".gz" ]; then
  echo "GZIP file detected; unzipping..."
  gunzip "$2"
  FLAG=true
  FILENAME="${2/.gz/$EMPTY_STR}"
fi

R-SAP --in1 $1 --in2 $FILENAME --outDir "$3/rsap_files" --rf GTF --tNum 4

if [ $FLAG == true ]; then
  echo "Re-zipping GZIP file..."
  gzip "$FILENAME"

#1 Reference genome alignment file (sorted on RNA‐Seq IDs) of the RNA‐Seq reads in psl format.
#Reference transcript annotation file (UCSC table browser, BED or GTF format). 
#We are looking for rsap_files/ChimericTranscriptAnnotation.out
