#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH -p RM
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=26
#SBATCH --output=sbatch.log
#SBATCH --chdir=/pylon5/mc5frap/siweixu/project

source src/load_module.sh

#
# Controller script for the pipeline
#

# Default general parameters
DEF_CPU=`grep -c ^processor /proc/cpuinfo`
DEF_CPU=`expr $DEF_CPU - 2`
DEF_MEMORY=`awk '/MemTotal/ {print $2}' /proc/meminfo`
DEF_MEMORY=`expr $DEF_MEMORY / 1024 / 1024 / 2`
DEF_MEMORY+="G"
DEF_OUTPUT_DIR="output"

# Default trinity parameters
DEF_TRIMMOMATIC_PARAM="ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/TruSeq3-PE.fa:2:30:12"

# Initialize variables
LEFT="/dev/null"
RIGHT="/dev/null"
GENOME="/dev/null"
COORD="/dev/null"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SECONDS=0
GFUSION_OPT=""
GFUSION_BOWTIE_INDEX=""
RSEM_OPT=""

# Load software path from ini file
. src/config.ini

# String literals
VERSION_TEXT="Version: 1.0.0"
HELP_TEXT="Usage: ./main.sh -l|--left <left_seq> -r|--right <right_seq> -g|--genome <genome_seq> -c|--coord <ref_coord> [-h|--help] [-v|--version] [-o|--output-dir <output_dir>] [--cpu <number_of_cores>] [--max-memory <max_memory>] [--trimmomatic-param <param>] [--GFusiion] [--GFusion-bowtie-index <bowtie_index_prefix>] [--GFusion-opt <gfusion_options>] [--rsem] [--rsem-opt <rsem_options>]
  
  Example:
    ./main.sh -l left.fq -r right.fq -g genome.fasta -c coord.gtf --cpu 8 --max-memory 10G --GFusion --GFusion-bowtie-index data/index/hg31
  
  Required:
    -l|--left <left_seq>:         The left sequence
    -r|--right <right_seq>:       The right sequence
    -g|--genome <genome_seq>:     The reference genome sequence
    -c|--coord <ref_coord>:       The reference transcript coordinate file in 
                                  GTF format or zipped GTF files 

  Optional:
    -h|--help:                    Show the help information for this pipeline
    -v|--version:                 Show the version information for this pipeline
    -o|--output-dir <output_dir>: The output directory (default is output/)
    --cpu <number_of_cores>:      The number of CPU cores to work (default is
                                  total number of cores minus 2)
    --max-memory <max_memory>:    The max number of memory to be used (default
                                  is half of the total memory)
    --GFusion:                    Perform genome alignment-based fusion gene
                                  detection powered by GFusion (default is 
                                  disabled) 
    --rsem:                       Perform RSEM estimation of gene and isoform 
                                  expression (default is disabled)
  
  Subroutine-specific Optional Parameters:
    --trimmomatic-param <param>:  Adapters and specifications for trimmomatic
                                  to perform automatic quality trim (default
                                  is \"ILLUMINACLIP:\$TRIMMOMATIC_DIR/
                                  adapters/TruSeq3-PE.fa:2:30:12\")
    --GFusion-bowtie-index <idx>: The path and prefix of bowtie index file. 
                                  NOTICE: This is REQUIRED if GFusion is enabled.
    --GFusion-opt <options>:      Optional arguments passed into GFusion. 
                                  'options' must be en entire string such 
                                  as \"-r 0 -n 1\"
    --rsem-opt <options>:         Optional arguments passed into RSEM. 
                                  'options' must be en entire string such 
                                  as \"-r 0 -n 1\"
"

# Parse options
OPTS=`getopt -o vhl:r:g:c:o: --long version,help,left:,right:,genome:,coord:,output-dir:,cpu:,max-memory:,GFusion,rsem,trimmomatic-param:,GFusion-bowtie-index:,GFusion-opt:,rsem-opt: -n 'parse-options' -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

eval set -- "$OPTS"

VERSION=false
HELP=false
GFUSION_FLAG=false
RSEM_FLAG=false

while true; do
  case "$1" in
    -v | --version ) VERSION=true; shift ;;
    -h | --help )    HELP=true; shift ;;
    -l | --left ) LEFT="$2"; shift; shift ;;
    -r | --right ) RIGHT="$2"; shift; shift ;;
    -g | --genome ) GENOME="$2"; shift; shift ;;
    -c | --coord ) COORD="$2"; shift; shift ;;
    -o | --output-dir ) DEF_OUTPUT_DIR="$2"; shift; shift ;;
    --cpu ) DEF_CPU="$2"; shift; shift ;;
    --max-memory ) DEF_MEMORY="$2"; shift; shift ;;
    --GFusion ) GFUSION_FLAG=true; shift ;;
    --rsem ) RSEM_FLAG=true; shift ;;
    --trimmomatic-param ) DEF_TRIMMOMATIC_PARAM="$2"; shift; shift ;;
    --GFusion-opt ) GFUSION_OPT=$2; shift; shift ;;
    --GFusion-bowtie-index ) GFUSION_BOWTIE_INDEX=$2; shift; shift ;;
    --rsem-opt ) RSEM_OPT=$2; shift; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done


# Display error messages
if [ $VERSION == true ]; then
  echo "$VERSION_TEXT"
  exit 1
fi

if [ $HELP == true ]; then
  echo "$HELP_TEXT"
  exit 1
fi

if [ $LEFT == "/dev/null" ] || [ $RIGHT == "/dev/null" ] || [ $GENOME == "/dev/null" ] || [ $COORD == "/dev/null" ]; then
  echo "Invalid arguments!

"
  echo "$HELP_TEXT"
  exit 127
fi

if [ $GFUSION_FLAG == true ] && [ $GFUSION_BOWTIE_INDEX == "" ]; then
  echo "Must provide bowtie index if GFusion is enabled!"
fi

# Step 1: Run Trinity with Trimmomatic
echo "Step 1: Run Trinity with Trimmomatic"
TRINITY_COMMAND="src/trinity.sh $LEFT $RIGHT $DEF_OUTPUT_DIR $DEF_CPU $DEF_MEMORY $DEF_TRIMMOMATIC_PARAM"
$TRINITY_COMMAND

# Step 2: Run BLAT and TRF/Repeat Masker
echo "Step 2: Run BLAT and Tandem Repeat Finder with RepeatMasker"
source "src/trf_rm_blat.sh" "$GENOME" "$DEF_OUTPUT_DIR/trinity/Trinity.fasta" "$DEF_OUTPUT_DIR/blat" "$DEF_CPU"

# Step 3: Run R-SAP
echo "Step 3: Run R-SAP"
source "src/rsap.sh" "$DEF_OUTPUT_DIR/blat/result.psl" "$COORD" "$DEF_OUTPUT_DIR" $DEF_CPU

echo "Finished! Results are in $DEF_OUTPUT_DIR/rsap_files/ChimericTranscriptAnnotation.out"
echo "Total runtime: $SECONDS seconds."

# Step 4 (optional): Run GFusion
if [ $GFUSION_FLAG == true ]; then
  echo "GFusion enabled"
  GFUSION_COMMAND="perl $GFUSION_PATH -o $DEF_OUTPUT_DIR/GFusion -1 $LEFT -2 $RIGHT -p $DEF_CPU -i $GFUSION_BOWTIE_INDEX -g $COORD $GFUSION_OPT"
  echo "Running $GFUSION_COMMAND"
  $GFUSION_COMMAND
fi

# Step 5 (optional): Run RSEM
if [ $RSEM_FLAG == true ]; then
  echo "RSEM enabled"
  BOWTIE_PATH=`which bowtie`
  RSEM_PREP_COMMAND="rsem-prepare-reference --gtf $COORD --bowtie2 --bowtie2-path $BOWTIE_PATH $GENOME -p $DEF_CPU"
  echo "RUnning $RSEM_PREP_COMMAND"
  $RSEM_PREP_COMMAND
  RSEM_CALC_COMMAND="rsem-calculate-expression --bowtie2 --bowtie2-path $BOWTIE_PATH --paired-end --append-names --output-genome-bam -p $DEF_CPU --calc-ci --ci-memory $DEF_MEMORY $LEFT $RIGHT"
  echo "Running $RSEM_CALC_COMMAND"
  $RSEM_CALC_COMMAND
fi
