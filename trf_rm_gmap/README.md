## Tandem Repeat Finder, RepeatMakser and GMAP

This folder contains trf_rm_gmap.sh, which process the fasta file to Tandem Repeat Finder, RepeatMakser, and GMAP.  
The Tandem Repeat Finder masks the tandem repeats to Ns, while the RepeatMasker masks the interspersed repeats to Ns. GMAP then aligns the result to genome file, outputing a PSL format file.  
To learn more about those three software, check the links below:
[Tandem Repeat Finder](https://tandem.bu.edu/trf/trf.html)  
[RepeatMasker](http://www.repeatmasker.org/)   
[GMAP](http://research-pub.gene.com/gmap/)  
The result would be in `result.psl`.  
To test the program, simply run:   
`source trf_rm_gmap genome.fasta test.fasta` in this directory.
