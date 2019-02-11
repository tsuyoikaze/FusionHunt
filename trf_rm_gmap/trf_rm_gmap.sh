#!/bin/bash

# For testing the program, run source trf_rm_gmap.sh genome.fasta test.fasta

# Usage: 
# First, load modules using the following command: (in your terminal)
# source project/src/load_module.sh

# Second, run this file in the format:
# source YOUR_DIR_TO_THIS_FILE/trf_rm_gmap.sh YOUR_DIR_TO_GENOME/GENOME_FILE YOUR_DIR_TO_FASTA/FASTA_FILE

# The result will be in result.psl, in the same directory of tf_rm_gmap.sh

# ----------------------------------------------------------
GENOME_FILE=$1
FASTA_FILE=$2
printf "\n\nUsing genome file: $GENOME_FILE, fasta file: $FASTA_FILE\n\n"
printf "\n\nLet's run the tests --------------------------------------\n\n"

# Tandem Repeat Finder
printf "\n\n----------------Tandem Repeat Finder starts---------------\n\n"
trf $FASTA_FILE 2 7 7 80 10 50 500 -f -d -m
rm *.html
rm *.dat
printf "\n\nTandem Repeat Find for $FASTA_FILE Complete.\n\n"

# Repeat Masker
printf "\n\n----------------Repeat Masker starts---------------\n\n"
EXTENSION=".2.7.7.80.10.50.500.mask"
FASTA_FILE=$(basename "$FASTA_FILE")
FASTA_FILE="$FASTA_FILE$EXTENSION"
RepeatMasker -pa 10000 -q -a $FASTA_FILE
rm -r *RM_*
rm *.mask
rm *.cat
rm *.tbl
rm *.align
rm *.out
rm -r *output*
printf "\n\nInterspersed Repeat Find (Repeat Masker) for $FASTA_FILE Complete.\n\n"


# GMAP
printf "\n\n----------------GMAP starts---------------\n\n"
FASTA_FILE="$FASTA_FILE.masked"
gmap -g $GENOME_FILE -f 1 $FASTA_FILE > project/trf_rm_gmap/result.psl
rm *.masked
printf "\n\nGMAP finishes. The result is in result.psl.\n\n"
