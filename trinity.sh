#!/bin/bash
# Basic script to run trinity with trimmomatic
# Usage: ./trinity.sh left_seq right_seq output_dir cpu memory trimmomatic_quality_parameter

echo "Running Trinity --seqType fq --trimmomatic --left $1 --right $2 --CPU $4 --max_memory $5 --output $3/trinity --quality_trimming_params $6 --full_cleanup"

Trinity --seqType fq --trimmomatic --left "$1" --right "$2" --CPU "$4" --max_memory "$5" --output "$3/trinity" --quality_trimming_params "$6" --full_cleanup
