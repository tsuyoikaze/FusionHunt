#!/bin/bash
echo "Running: mkdir $3/rsap_files"

if [ ! -d "$3/rsap_files" ]; then
  mkdir "$3/rsap_files"
fi

FLAG=false
FILENAME=$2
EMPTY_STR=""
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. $DIR/config.ini

if [ ${FILENAME: -3} == ".gz" ]; then
  echo "GZIP file detected; unzipping..."
  gunzip "$2"
  FLAG=true
  FILENAME="${2/.gz/$EMPTY_STR}"
fi

$RSAP_PATH --in1 $1 --in2 $FILENAME --outDir "$3/rsap_files" --rf GTF --tNum "$4"

if [ $FLAG == true ]; then
  echo "Re-zipping GZIP file..."
  gzip "$FILENAME"

#1 Reference genome alignment file (sorted on RNA‐Seq IDs) of the RNA‐Seq reads in psl format.
#Reference transcript annotation file (UCSC table browser, BED or GTF format). 
#We are looking for rsap_files/ChimericTranscriptAnnotation.out
