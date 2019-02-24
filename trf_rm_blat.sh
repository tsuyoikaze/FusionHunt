#!/bin/bash

# File paths and current directory
GENOME_FILE=$1
FASTA_FILE=$2
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. $DIR/config.ini

echo "Using genome file: $GENOME_FILE, fasta file: $FASTA_FILE"

# Tandem Repeat Finder
echo "/n/n----------------Tandem Repeat Finder starts---------------/n/n"

TRF_CMD="$TRF_PATH $FASTA_FILE 2 7 7 80 10 50 500 -f -d -m"
echo "Running $TRF_CMD"
$TRF_CMD
rm *.html
rm *.dat

# Repeat Masker
echo "/n/n----------------Repeat Masker starts---------------/n/n"
EXTENSION=".2.7.7.80.10.50.500.mask"
FASTA_FILE=$(basename "$FASTA_FILE")
FASTA_FILE="$FASTA_FILE$EXTENSION"

RM_CMD="$REPEAT_MASKER_PATH -pa "$4" -q -a $FASTA_FILE"
echo "Running $RM_CMD"
$RM_CMD
rm *.mask
rm *.cat
rm *.out
rm *.mask.cat.all

# BLAT
echo "/n/n----------------BLAT starts---------------/n/n"
FASTA_FILE="$FASTA_FILE.masked"
BLAT_COMMAND="$BLAT_PATH $GENOME_FILE $FASTA_FILE $3/result.psl"
printf "Running $BLAT_COMMAND"
if [ ! -d $3 ]; then
  echo "Creating output folder..."
  mkdir $3
fi
$BLAT_PATH $GENOME_FILE $FASTA_FILE $3/result.psl
