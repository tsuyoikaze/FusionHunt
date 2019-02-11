## Tandem Repeat Finder, RepeatMakser and GMAP

This folder contains trf_rm_gmap.sh, which process the fasta file to Tandem Repeat Finder, RepeatMakser, and GMAP.    

The Tandem Repeat Finder masks the tandem repeats to Ns, while the RepeatMasker masks the interspersed repeats to Ns. GMAP then aligns the result to genome file, outputing a PSL format file.    

To learn more about those three software, check the links below:  
[Tandem Repeat Finder](https://tandem.bu.edu/trf/trf.html)  
[RepeatMasker](http://www.repeatmasker.org/)   
[GMAP](http://research-pub.gene.com/gmap/)    


## Usage
First, load modules using the following command: (in your terminal)  
`source project/src/load_module.sh`  

Second, run this file in the format:  
`source YOUR_DIR_TO_THIS_FILE/trf_rm_gmap.sh YOUR_DIR_TO_GENOME/GENOME_FILE YOUR_DIR_TO_FASTA/FASTA_FILE`  

The result would be in `result.psl`.  

To test the program with sample data, simply run:   
`source trf_rm_gmap genome.fasta test.fasta` in this directory.

