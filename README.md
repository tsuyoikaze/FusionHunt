# FusionHunt: A Fusion Gene Detection Software
Bioinformatics Data Integration Practicum Project - Team 2

## Description

This software, the FusionHunt, is a pipeline that detects possible fusion genes based on *de novo* transcriptome assembly from human RNA-seq data with benchmark result and functionality as existing software [GFusion](https://github.com/xiaofengsong/GFusion) for clean and efficient fusion gene detection.


The pipeline runs by first assembling the reads by *de novo* transcriptome assembly using Trinity, removing tandem repeats and interspersed repeats using Tandem Repeat Finder and RepeatMasker, then aligning to the genome using BLAT, detecting possible chimeric transcripts using R-SAP, and filtering out for more likely fusion genes by aligning back to the original reads using Bowtie. Then, the pipeline will compare its results directly to those from GFusion, a pipeline that uses reference-based assembly to detect fusion genes.

## Prerequisites

* Linux environment with bash
* [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) v2.8.4 or higher with trimmomatic embedded
* [Tophat](https://ccb.jhu.edu/software/tophat/index.shtml) v2.1.1 or higher
* [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) v1.1.2
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) v2.2.7
* [Samtools](http://samtools.sourceforge.net/) v1.3.1
* [Jellyfish](https://www.cbcb.umd.edu/software/jellyfish/) v2.2.6 or higher
* [Salmon](https://combine-lab.github.io/salmon/) v0.9.1
* [GCC 8.2.0](https://gcc.gnu.org/) or higher
* [Tandem Repeat Finder](https://tandem.bu.edu/trf/trf.html) v4.09 or higher
* [RepeatMasker](http://www.repeatmasker.org/) v1.332 or higher
* [RMBlast](http://www.repeatmasker.org/RMBlast.html) v2.60+ or higher
* [ActivePerl](https://www.activestate.com/products/activeperl/) v5.8.0 or higher
* [Perl Bioperl](https://bioperl.org/) module v1.007 or higher
* [R-SAP](http://www.mcdonaldlab.biology.gatech.edu/r-sap.htm) v1.1 or higher
* [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)
* [GFusion](https://github.com/xiaofengsong/GFusion) v1.0 (optional)
* [RSEM](https://deweylab.github.io/RSEM/) v1.1.17 or higher (optional)

## Configuration

Although some software being used in this pipeline are by default installed and added to `$PATH`, some require manual specification of path information. 

Please make sure that the following is true

* Trinity is accessible via `Trinity` command
* Tophat is accessible via `tophat` command
* Bowtie is accessible via `bowtie` command
* Samtools is accessible via `samtools` command
* RSEM is accessible in command line directly (optional)
* ActivePerl with BioPerl module installed is accessible via `perl` command. (You may need to add ActivePerl to your `PATH` in order to do that)
* Trimmomatic installation directory is accessible via the `$TRIMMOMATIC_DIR` environment variable

*Then, please fill the `config.ini` file with the installation paths of the remaining software in order to get the pipeline working as intended.*

## Usage

You should call `main.sh` with the following command line arguments:


`Usage: ./main.sh -l|--left <left_seq> -r|--right <right_seq> -g|--genome <genome_seq> -c|--coord <ref_coord> [-h|--help] [-v|--version] [-o|--output-dir <output_dir>] [--cpu <number_of_cores>] [--max-memory <max_memory>] [--trimmomatic-param <param>] [--GFusiion] [--GFusion-bowtie-index <bowtie_index_prefix>] [--GFusion-opt <gfusion_options>] [--rsem] [--rsem-opt <rsem_options>]
`


```  

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
                                  is "ILLUMINACLIP:$TRIMMOMATIC_DIR/
                                  adapters/TruSeq3-PE.fa:2:30:12")
    --GFusion-bowtie-index <idx>: The path and prefix of bowtie index file.
                                  NOTICE: This is REQUIRED if GFusion is enabled.
    --GFusion-opt <options>:      Optional arguments passed into GFusion.
                                  'options' must be en entire string such
                                  as "-r 0 -n 1"
    --rsem-opt <options>:         Optional arguments passed into RSEM.
                                  'options' must be en entire string such
                                  as "-r 0 -n 1"


```

## Output
The output will be in `<output_dir>/rsap_files/ChimericTranscriptAnnotation.out`.
## FASTQ File formatting

Please use a file format that is supported by Trinity. If you are downloading data from SRA, make sure to use this command:
`fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files <accession_ID>`

## Genome downloads:

* Human Genome (hg38/GRCh38) from UCSC in fasta/fa format can be downloaded [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) for the -g flag
* Human Genome (hg38/GRch38) transcript file in gtf format can be downloaded [here](https://genome.ucsc.edu/cgi-bin/hgTables) for the -c flag

## Known bugs and issues

* **Bugs from prerequisite software**:

There is a bug in Trinity v2.2.8-v2.3.4.1 that puts the entire program into an infinite loop. Please use Trinity v2.2.7

* **Bugs within this pipeline**:

There are currently no known bugs at this point. However, should there be any issues please submit it to the issues page. 

## The group

* **Program Manager**: Sarah Hsu
* **Technical Lead**: Siwei Xu
* **Lead Technical Writer**: Martin Ma
* **Q&A Lead**: Juntao Chen
* **Communications Lead**: Saina Mahera

