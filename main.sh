#!/bin/bash
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

# Load software path from ini file
. $DIR/config.ini

# String literals
VERSION_TEXT="Version: 1.0.0"
HELP_TEXT="Usage: ./main.sh -l|--left <left_seq> -r|--right <right_seq> -g|--genome <genome_seq> -c|--coord <ref_coord> [-h|--help] [-v|--version] [-o|--output-dir <output_dir>] [--cpu <number_of_cores>] [--max-memory <max_memory>] [--trimmomatic-param <param>]
  
  Example:
    ./main.sh -l left.fq -r right.fq -g genome.fasta -c coord.gtf --cpu 8 --max-memory 10G
  
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
  
  Subroutine-specific Optional Parameters:
    --trimmomatic-param <param>:  Adapters and specifications for trimmomatic
                                  to perform automatic quality trim (default
                                  is \"ILLUMINACLIP:\$TRIMMOMATIC_DIR/
                                  adapters/TruSeq3-PE.fa:2:30:12\")
"

# Parse options
OPTS=`getopt -o vhl:r:g:c:o: --long version,help,left:,right:,genome:,coord:,output-dir:,cpu:,max-memory:,trimmomatic-param: -n 'parse-options' -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

eval set -- "$OPTS"

VERSION=false
HELP=false

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
    --trimmomatic-param ) DEF_TRIMMOMATIC_PARAM="$2"; shift; shift ;;
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

# Step 1: Run Trinity with Trimmomatic
echo "Step 1: Run Trinity with Trimmomatic"
TRINITY_COMMAND="$DIR/trinity.sh $LEFT $RIGHT $DEF_OUTPUT_DIR $DEF_CPU $DEF_MEMORY $DEF_TRIMMOMATIC_PARAM"
$TRINITY_COMMAND

# Step 2: Run GMAP and TRF/Repeat Masker
echo "Step 2: Run BLAT and Tandem Repeat Finder with RepeatMasker"
source "$DIR/trf_rm_gmap.sh" "$GENOME" "$DEF_OUTPUT_DIR/trinity/Trinity.fasta" "$DEF_OUTPUT_DIR/gmap" "$DEF_CPU"

# Step 3: Run R-SAP and Bowtie
echo "Step 3: Run R-SAP and Bowtie"
source "$DIR/rsap.sh" "$DEF_OUTPUT_DIR/gmap/result.psl" "$COORD" "$DEF_OUTPUT_DIR" $DEF_CPU
