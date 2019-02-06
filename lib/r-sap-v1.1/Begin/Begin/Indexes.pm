use strict;
use warnings;

package Begin::Indexes;

	sub pslx_indx{
	my %pslxIndx = (
	"matches" => 0,
	"misMatch" => 1,
	"repMatch" => 2,
	"nCount" => 3,
	"qNumIn" => 4,
	"qBsIn" => 5,
	"tNumIn" => 6,
	"tBsIn" => 7,
	"strand" => 8,
	"qID" => 9,
	"qLen" => 10,
	"qSt" => 11,
	"qEn" => 12,
	"tID" => 13,
	"tLen" => 14,
	"tSt" => 15,
	"tEn" => 16,
	"blNum" => 17,
	"blLen" => 18,
	"qBlSt" => 19,
	"tBlSt" => 20,
	"qSeq" => 21,
	"tSeq" => 22,
	);
	return(%pslxIndx);  ## returning the index containing hash
	}  ## function end
	###########################	

	sub pslx2_indx{
	my %pslx2Indx;  ## hash containing array indexes of different columns in a pslx2 file
	%pslx2Indx = (
		"score" => 0,
		"coverage" => 1,
		"identity" => 2,
		"matches" => 3,
		"misMatch" => 4,
		"repMatch" => 5,
		"nCount" => 6,
		"qNumIn" => 7,
		"qBsIn" => 8,
		"tNumIn" => 9,
		"tBsIn" => 10,
		"strand" => 11,
		"qID" => 12,
		"qLen" => 13,
		"qSt" => 14,
		"qEn" => 15,
		"tID" => 16,
		"tLen" => 17,
		"tSt" => 18,
		"tEn" => 19,
		"blNum" => 20,
		"blLen" => 21,
		"qBlSt" => 22,
		"tBlSt" => 23,
		"qSeq" => 24,
		"tSeq" => 25,
		);
	return(%pslx2Indx);
	}  ## function ends
	###########################
	sub ucsc_known_gene_indexes{
	my %knownGeneIndx;
	%knownGeneIndx = (
		"transcriptID" => 0,
		"chr" => 1,
		"strand" => 2,
		"txSt" => 3,
		"txEn" => 4,
		"cdsSt" => 5,
		"cdsEn" => 6,
		"exCount" => 7, 
		"exSt" => 8,
		"exEn" => 9,
		"geneSymbol" => 10,
		);

	return(%knownGeneIndx);
	}  ## function ends
	###########################
	sub chimer_indx{
	my %chimerIndexes;
	%chimerIndexes = (
		"qID" => 0,
		"qLen" => 1,
		"chrs" => 2,
		"identities" => 3,
		"qSts" => 4,
		"gapLen" => 5,
		"overlap" => 6,
		"indx" => 7,
		"qCov1" => 8,
		"qCov2" => 9,
		"score1" => 10,
		"coverage1" => 11,
		"strand1" => 12,
		"tLen1" => 13,
		"tSt1" => 14,
		"tEn1" => 15,
		"blNum1" => 16,
		"blLen1" => 17,
		"qBlSt1" => 18,
		"tblSt1" => 19,
		"score2" => 20,
		"coverage2" => 21,
		"strand2" => 22,
		"tLen2" => 23,
		"tSt2" => 24,
		"tEn2" => 25,
		"blNum2" => 26,
		"blLen2" => 27,
		"qBlSt2" => 28,
		"tblSt2" => 29,
		);
	return(%chimerIndexes);
	}  ## function ends
	###########################
	sub ref_gene_array_info{  ## to store reference gene info (read counts, read characterization info)
	my %geneFileIndx = (
		"transcriptID" => 0,	   ## transcript ID (each reference transcript ID may be related to multiple genomic loci. So the unique should be different than just the transcriptID).
		"chr" => 1,		   ## transcripte chromosome
		"txSt" => 2,		   ## transcript ends
		"txEn" => 3,		   ## transcript end genomic position (transcription 
		"transcriptLen" => 4,      ## length of the transcript
		"codingFlag" => 5,         ## coding or non-coding gene
		"readCount" => 6,          ## total reads mapped within the exons
		"exonReadCounts" => 7,	   ## readCount for each exon
		"exonLengths" => 8,        ## length of the exons
		"skippedExons" => 9,       ## skipped exons
		"intronReadCounts" => 10,  ## total reads mapped within the introns
		"intronLengths" => 11,	   ## length of the introns
		"retainedIntrons" => 12,	   ## retained introns
		"genicRegions" => 13,	   ## genic regions affected due to aberrant splicing (5'UTR, CDS, 3'UTR or Non-coding)
		"annotTags" => 14,	   ## annotation tags for the gene (based on the reads mapped to it)
		"multiAnnotFlag" => 15,	   ## multiple annotation flag (0 by default)
		);
	return(%geneFileIndx);
	}  ## function ends
	###########################
	sub expressed_gene_info{
	my %readFileIndx = (
		"transcriptID" => 0,
		"chr" => 1,
		"txSt" => 2,
		"spannedExons" => 3,
		"spannedIntrons" => 4,
		"skippedExons" => 5,
		"retainedIntrons" => 6,
		"genicRegions" => 7,
		"annotTags" => 8,
		"readID" => 9,
		);
	return(%readFileIndx);
	}  ## function ends
	###########################


1;




##match   mis   rep       N      qNumIn   qBsIn   tNumIn   tBsIn  strand    qID    
##0	  1	2	  3	  4	   5	  6	   7	   8	    9	   

##qLen    qSt     qEn     tID      tSize	tSt             tEn             blNum   blLen   qBlSt   tBlEn	qSeq                           tSeq
##10	  11	  12	  13	   14		15		16		 17	18	19	 20	21			   22

#ucID		chr	strand	txSt	txEn	cdsSt	cdsEn	exCount	exSt		exEn		codFlag		geneSymbol
#0		1	2	3	4	5	6	7	8		9		10		11
#uc001aaa.2	chr1	+	1115	4121	1115	1115	3	1115,2475,3083,	2090,2584,4121,	noncoding	BC032353

## 0:transcriptID, 1: exon mapped read count, 2: exons mapped by read, 3: annotation tag, 4: read ID
##NM_001042465	1	15,15,15,	ExonsOnly,	FIC9KKY04JEZDI
