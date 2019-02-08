#!/bin/bash
mkdir rsap_files
perl R-SAP --in1 $1 --in2 $2 --outDir rsap_files --rf GTF --tNum 4

#1 Reference genome alignment file (sorted on RNA‐Seq IDs) of the RNA‐Seq reads in psl format.
#2 Reference transcript annotation file (UCSC table browser, BED or GTF format). 
#We are looking for rsap_files/ChimericTranscriptAnnotation.out
