
package file_parsers::VerifyPSLRecords;


use strict;
use warnings;


use shared::CheckStringsNumericals;

	sub verify_psl_records{
	my($inFile,$linesToCheck) = @_;  ## input file and number of the lines in the beginning to check
	open FHVerifyPslRecordsIn,$inFile or die "\n can not open file $inFile \n";

	my($str,@a1,$n,$formatFlagPass);
	$n = 0;
	$formatFlagPass = 1;
	
		while($str = <FHVerifyPslRecordsIn>){
		last if($n > $linesToCheck);  ## format checkinf done
		next if($str =~ m/^(\#)/);  ## skipping the header
		$str =~ s/\n//;
		$str =~ s/\r//;
		@a1 = split(/\t+/,$str);
		++$n;
			if(!(check_if_psl(\@a1))){
			$formatFlagPass = 1;
			last;
			}
		}  ## while($str = <FHVerifyPslRecordsIn>) ends
	close FHVerifyPslRecordsIn;
	return($formatFlagPass);
	}  ## function ends
	###########################
	sub check_if_psl{
	my($a1Ref) = @_;
	my($ifPsl,$i);
	$ifPsl = 1;  ## by default true
	
		if( !((scalar(@{$a1Ref}) == 21) || (scalar(@{$a1Ref}) == 23)) ){  ## if record length is neither psl nor pslx
		$ifPsl = 0;
		}
		else{
			for($i = 0;$i<=17;++$i){
				if( ($i == 8) || ($i == 9) || ($i == 13) ){  ## non-numerical values in the record
				next;
				}
				if(!(check_if_numerical($a1Ref->[$i]))){
				$ifPsl = 0;
				last;
				}
			}  ## for($i) ends

			for($i = 18;$i<=20;++$i){  ## alignment block lengths, block starts for query and target
			$a1Ref->[$i] =~ s/\D//g;  ## removing non numeric characters
				if(!(check_if_numerical($a1Ref->[$i]))){  ## exon blocks starts
				$ifPsl = 0;
				}
			}  ## for($i) ends
		
			if( !(($a1Ref->[8] eq '+') || ($a1Ref->[8] eq '-')) ){  ## supposed to be alignment direction (or strand may be)
			$ifPsl = 0;
			}
		}  ## else ends
	return($ifPsl);
	}  ## function ends
	###########################
1;


##0	1	2	3	4	5	6	7	8	9					10	11	12	13	14		15		16		17	18	19	20		21										22
##71	2	0	0	0	0	0	0	+	HWI-EAS229_75_30DY0AAXX:7:1:0:1547/1	75	1	74	chr14	106368585	105244564	105244637	1	73,	1,	105244564,	gtccacctccgccatgacaacagacaaattgacatgggtggggttacccgccaagcggtcgatggtcttctgt,	gtccacctccgccatgacaacagacacattgacatgggtgggtttacccgccaagcggtcgatggtcttctgt,