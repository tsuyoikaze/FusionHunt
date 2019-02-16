#!/usr/bin/bash

# This script adds the required modules into the PSC system
# This is not required for final submission. Only for PSC machines!

# Manually add Bioperl to perl library and ActivePerl
export PATH=/pylon5/mc5frap/siweixu/project/lib/ActivePerl-5.26/bin:$PATH
export PERL5LIB=/pylon5/mc5frap/siweixu/project/lib/r-sap-v1.1

# Get current directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Restore software path from backup
cp "$DIR/config.ini.psc" "$DIR/config.ini"


# PSC built-in software
module load trinity/2.8.4
module load tophat/2.1.1
module load bowtie/1.1.2
module load bowtie2/2.2.7
module load samtools/1.3.1
module load jellyfish2/2.2.6
module load salmon/0.9.1
module load gcc/8.2.0

# add variables
export TRIMMOMATIC_DIR=/opt/packages/trinity/2.8.4/trinity-plugins/Trimmomatic

