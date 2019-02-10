# 03-713 Bioinformatics Data Integration Practicum Project - Team 2

## Description
A pipeline that finds fusion genes base on *De novo transcriptome assembly* with benchmark functionality to compare results with existing software: [GFusion](https://github.com/xiaofengsong/GFusion) for refined and efficient fusion gene detection. 

## Prerequisites
* Linux environment with bash
* Trinity v2.8.4 or higher with trimmomatic embedded
* Tophat v2.1.1 or higher
* Bowtie v1.1.2
* Samtools v0.1.19
* Tandem Repeat Finder v4.09 or higher
* RepeatMasker v1.332 or higher
* RMBlast v2.60+ or higher
* ActivePerl v5.8.0 or higher (can be obtained [here](https://www.activestate.com/products/activeperl/))
* Perl Bioperl module v1.007 or higher
* R-SAP v1.1 or higher
* GMAP v2019-01-24 or higher
* GFusion v1.0

## Configuration
As this is a pipeline written as pure shell script, there is no need to configurate the software itself. However, please make sure that the following command works:

* Trinity is accessible via `Trinity` command
* Tophat is accessible via `tophat` command
* Bowtie is accessible via `bowtie` command
* Samtools is accessible via `samtools` command
* RepeatMasker is accessible via `RepeatMasker` command and is configurated correctly to link with Tandem Repeat Finder and RMBlast
* ActivePerl with BioPerl module installed is accessible via `perl` command. (You may need to add ActivePerl to your `PATH` in order to do that)
* R-SAP is accessible via `R-SAP` command
* GMAP is accessible via `gmap` command
* GFusion is accessible via `GFusion` command. (Since the program itself is a perl script, you may need to add `alias GFusion="perl /path/to/GFusion.pl"`)

## Usage

## Known bugs

## The group
* **Program Manager**: Sarah Hsu
* **Technical Lead**: Siwei Xu
* **Lead Technical Writer**: Martin Ma
* **Q&A Lead**: Juntao Chen
* **Communications Lead**: Saina Mahera
