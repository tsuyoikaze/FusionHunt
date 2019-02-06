package file_parsers::CheckIfGeneTrack;

require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/check_if_gene_track/;

use strict;
use warnings;

use shared::CheckStringsNumericals;

	sub check_if_gene_track{
	my($binColumn,@a1) = @_;
	my($ifGeneTrack,$i,$str);
	my $numericalValues = 0;

	$ifGeneTrack = 1;  ## true by default

		if(scalar(@a1) < (10 + $binColumn)){  ## array length is not appropriate
		$ifGeneTrack = 0;
		##die "\n array length is not appropriate for $a1[0+$binColumn] \n";
		}
		else{  ## if length of the record is appropriate

			for($i = 3;$i<=7;++$i){

				if(check_if_numerical($a1[$i+$binColumn])){
				$numericalValues += 1;
				}  ## if(numerical) ends
				else{
				##print "\n $i ($a1[$i+$binColumn]) is non numeric in $a1[0+$binColumn] \n";
				}

			}  ## for(i:3-7) ends
			if($numericalValues != 5){
			$ifGeneTrack = 0;
			##die "\n Error: not enough numerical values in $a1[0+$binColumn] \n";

			}
			elsif($numericalValues == 5){
				if( ($a1[3+$binColumn] > $a1[4+$binColumn]) || ($a1[5+$binColumn] > $a1[6+$binColumn]) ||
					($a1[5+$binColumn] < $a1[3+$binColumn]) || ($a1[6+$binColumn] > $a1[4+$binColumn]) ){
				$ifGeneTrack = 0;
				}
			}  ## elsif ends
			for($i = (8+$binColumn);$i<=(9+$binColumn);++$i){
			$a1[$i] =~ s/\D//g;  ## removing non numeric characters
				if(!(check_if_numerical($a1[$i]))){  ## alignment block cooridnates are not numeric
				$ifGeneTrack = 0; 
				}
			}  ## for(i:8-9) ends
			if( ($a1[2+$binColumn] ne '+') && ($a1[2+$binColumn]  ne '-') ){
			$ifGeneTrack = 0;
			}
		}  ## else(array length appropriate) ends
	return($ifGeneTrack);
	}  ## function ends
	###########################
1;

