package file_parsers::AllFieldsToDefault;

use shared::CheckStringsNumericals;
use file_parsers::CheckIfGeneTrack;


use strict;
use warnings;


	sub all_fields_to_default{
	my($geneFile,$outFile) = @_;  ## input gene file and output file name
	my $binColumn  = check_if_bin_column_exists($geneFile);
	my $geneNameIndex = get_gene_name_index($geneFile,$binColumn);  ## getting array index of the gene name/symbol (if any)
	#print "\nbinColumn: $binColumn \n";
	#print "geneNameIndex: $geneNameIndex \n";
	open FHAllFieldsToDefaultIn,$geneFile or die "\n can not open file $geneFile \n";  ## input gene track file
	open FHAllFieldsToDefaultOut,">$outFile";  ## out file
	my(@a1,@a2,$str,$i);
		while($str = <FHAllFieldsToDefaultIn>){
		$str =~ s/\n//;
		$str =~ s/\r//;
			if( ($str =~ m/^(\#+)/) || !($str) ){  ## skipping the header blank line
			next;
			}
		@a1 = split(/\t+/,$str);
			if(check_if_gene_track($binColumn,@a1)){  ## if record is gene track
			@a2 = @a1[(0+$binColumn)..(9+$binColumn)];
				if( ($geneNameIndex) && (scalar(@a1) >= ($geneNameIndex+1)) ){
				push(@a2,$a1[$geneNameIndex]);
				}
				elsif(scalar(@a1) > (9+$binColumn)){
				push(@a2,$a1[9+$binColumn+1]);
				}
				else{
				push(@a2,'.');
				}
			print FHAllFieldsToDefaultOut join("\t",@a2),"\n";
			splice(@a2);
			}  ## if(gene_track_ends())
			else{
			die "\n file format not recognized at line $. in gene track file $geneFile \n";
			}
		splice(@a1);
		}  ## while($str = <FHAllFieldsToDefaultIn>)
	close FHAllFieldsToDefaultIn;
	close FHAllFieldsToDefaultOut;
	}  ## function ends
	###########################
	sub check_if_bin_column_exists{
	my($fin) = @_;
	my($str,@a1,$i,$binColumn);
	$binColumn = 0;
	open FHParseFileAllFields,$fin or die "\n can not open file $fin \n";
	$str = <FHParseFileAllFields>;
	$str =~ s/^(\s+)//;
	$str =~ s/(\s+)$//;
		while($str =~ m/^(\#+)/){  ## if header exists
		@a1 = split(/\s+/,$str);
			if($a1[0] =~ m/(bin)/i){  ## 
			$binColumn = 1;
			}
		$str = <FHParseFileAllFields>;
		$str =~ s/^(\s+)//;
		$str =~ s/(\s+)$//;
		}  ## while ends

		if(!($binColumn)){
		$i = 1;
			while($str = <FHParseFileAllFields>){
			$i += 1;
			$str =~ s/^(\s+)//;
			$str =~ s/(\s+)$//;
			@a1 = split(/\s+/,$str);
				if( (check_if_numerical($a1[0])) && (($a1[3] eq '+') || ($a1[3] eq '-')) ){
				$binColumn = 1;
				}
			last if($i >= 10);
			}  ## while($str = <FHParseFileAllFields) ends
		}  ## if(!($binColumn)) ends

	close FHParseFileAllFields;
	return($binColumn);
	}  ## function ends
	###########################
	sub get_gene_name_index{
	my($fin,$binColumn) = @_;
	my($str,@a1,$i);	  
	my $geneNameIndex = 0;

	open FHGetGeneName,$fin or die "\n can not open file $fin \n";
	$str = <FHGetGeneName>;
	$str =~ s/^(\s+)//;
	$str =~ s/(\s+)$//;

		while($str =~ m/^(\#+)/){  ## while heade exists
		@a1 = split(/\s+/,$str);
			if(scalar(@a1) > (10+$binColumn)){  ## if there is any extra information in the array
				for($i = (10+$binColumn);$i<scalar(@a1);++$i){
					if( ($a1[$i] =~ m/(name)/i) || ($a1[$i] =~ m/(symbol)/i) ){
					$geneNameIndex = $i;
					last;
					}
				}  ## ofr($i) ends
				if($geneNameIndex){
				last;
				}
			}  ## if(scala(@a1) > 10+binColumn) ends
		$str = <FHGetGeneName>;
		}  ## while($str = <FHGeneName>) ends
	close FHGetGeneName;
	return($geneNameIndex);
	}  ## function ends
	###########################
1;

	

#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
#	0		1	2	3		4		5		6		7	8																				9
#0	NM_024763	chr1	-	67051159	67163158	67052400	67163102	17	67051159,67060631,67065090,67066082,67071855,67072261,67073896,67075980,67078739,67085754,67100417,67109640,67113051,67129424,67131499,67143471,67162932,	67052451,67060788,67065317,67066181,67071977,67072419,67074048,67076067,67078942,67085949,67100573,67109780,67113208,67129537,67131684,67143646,67163158,	0	WDR78	cmpl	cmpl	0,2,0,0,1,2,0,0,1,1,1,2,1,2,0,2,0,
#0	NM_207014	chr1	-	67075869	67163158	67075923	67163102	10	67075869,67078739,67085754,67100417,67109640,67113051,67129424,67131499,67143471,67162932,	67076067,67078942,67085949,67100573,67109780,67113208,67129537,67131684,67143646,67163158,	0	WDR78	cmpl	cmpl	0,1,1,1,2,1,2,0,2,0,

## if first line #
	## if a1[0] =~ m/bin/i
		## binColumn = TRUE
##else if a1[0] = NUM and strand = 3
	## binColumn = TRUE
##else
## continue assuming default

