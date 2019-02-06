package file_parsers::CheckIfBED;

require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/check_if_bed/;

use strict;
use warnings;

use shared::CheckStringsNumericals;


	sub check_if_bed{
	my(@a1) = @_;
	my($ifBED,$numericalValues);
	$ifBED = 1;  ## by default BED
	$numericalValues = 0;  ## minimum number of required numerical values

		if(scalar(@a1) < 12){  ## record is not long enough to be a BED record
		$ifBED = 0;
		##print "\n length not appropriate \n";
		}
		else{  ## if length of the record is >= 12
			for(my $i = 1;$i<=7;++$i){
				if( ($i == 3) || ($i == 4)  || ($i == 5) ){
				next;
				}
				if(check_if_numerical($a1[$i])){
				++$numericalValues;
				}
			}  ## for($i) ends
			if($numericalValues != 4){  ## geneSt, geneEn, cdsSt, cdsEn (all should be numerical)
			$ifBED = 0;
			##print "\n not numerical \n";
			}
			elsif($numericalValues == 4){
				if( ($a1[1] >= $a1[2]) || ($a1[6] > $a1[7]) || ($a1[6] < $a1[1]) ||
					($a1[7] > $a1[2]) ){
				$ifBED = 0;
				##print "\n out of range \n";
				}
			}  ## elsif ends
			if(check_if_numerical($a1[5])){  ## supposed to be DNA strand
			$ifBED = 0;
			##print "\n strand is numeric \n";
			}
			if(!(check_if_numerical($a1[9]))){  ## number of exons
			$ifBED = 0;
			##print "\n exon numbers is not numeric \n";
			}
		$a1[10] =~ s/\D//g;  ## removing non numeric characters
			if(!(check_if_numerical($a1[10]))){  ## exon blocks starts
			$ifBED = 0;
			##print "\n blockLength is not numeric \n";
			}
		$a1[11] =~ s/\D//g;  ## removing non numeric characters
			if(!(check_if_numerical($a1[11]))){  ## exon blocks ends
			$ifBED = 0;
			##print "\n exon block start are not numeric \n";
			}
		}  ## else(arrayLen >= 12) ends
	return($ifBED);
	}  ## function ends
	###########################
1;

