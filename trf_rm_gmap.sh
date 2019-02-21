#!/bin/bash

# For testing the program, run source trf_rm_gmap.sh genome.fasta test.fasta

# Usage: 
# First, load modules using the following command: (in your terminal)
# source project/src/load_module.sh

# Second, run this file in the format:
# source YOUR_DIR_TO_THIS_FILE/trf_rm_gmap.sh YOUR_DIR_TO_GENOME/GENOME_FILE YOUR_DIR_TO_FASTA/FASTA_FILE

# The result will be in result.psl, in the same directory of tf_rm_gmap.sh

# ----------------------------------------------------------

# File paths and current directory
GENOME_FILE=$1
FASTA_FILE=$2
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. $DIR/config.ini

echo "\n\nUsing genome file: $GENOME_FILE, fasta file: $FASTA_FILE\n\n"

# Tandem Repeat Finder
echo "\n\n----------------Tandem Repeat Finder starts---------------\n\n"

TRF_CMD="$TRF_PATH $FASTA_FILE 2 7 7 80 10 50 500 -f -d -m"
echo "Running $TRF_CMD"
$TRF_CMD
rm *.html
rm *.dat

# Repeat Masker
echo "\n\n----------------Repeat Masker starts---------------\n\n"
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
echo "\n\n----------------BLAT starts---------------\n\n"
FASTA_FILE="$FASTA_FILE.masked"
GMAP_COMMAND="$GMAP_PATH -D $GENOME_FILE -f 1 $FASTA_FILE > $3/result.psl"
echo "Running $GMAP_COMMAND"
if [ ! -d $3 ]; then
  echo "Creating output folder..."
  mkdir $3
fi
$GMAP_PATH -g $GENOME_FILE -f 1 $FASTA_FILE > $3/result.psl

# # GMAP
# echo "\n\n----------------GMAP starts---------------\n\n"
# FASTA_FILE="$FASTA_FILE.masked"
# GMAP_COMMAND="$GMAP_PATH -D $GENOME_FILE -f 1 $FASTA_FILE > $3/result.psl"
# echo "Running $GMAP_COMMAND"
# if [ ! -d $3 ]; then
#   echo "Creating output folder..."
#   mkdir $3
# fi
# $GMAP_PATH -g $GENOME_FILE -f 1 $FASTA_FILE > $3/result.psl

