package file_parsers::CheckIfGTF;

use strict;
use warnings;

use shared::CheckStringsNumericals;

	sub check_if_GTF{
	my($a1Ref) = @_;
	my $ifGTF = 1;  ## default is true

		if(scalar(@{$a1Ref}) < 9){  ## array is length is not appropriate
		$ifGTF = 0;
		}
		else{
			if( !(check_if_numerical($a1Ref->[3])) || !(check_if_numerical($a1Ref->[4])) ){  ## start and end should be numerical
			$ifGTF = 0;
			}
			else{
				if($a1Ref->[3] > $a1Ref->[4]){  ## start is greater than end of the feature
				$ifGTF = 0;
				}
			}
			if(check_if_numerical($a1Ref->[2])){  ## feature is not a string
			$ifGTF = 0;
			}
			
			##if( !($a1Ref->[6] eq '+') && !($a1Ref->[6] eq '-') ){  ## strand should be eigher '+' or '-' for genes
			##$ifGTF = 0;
			##}
			
			if($a1Ref->[8] !~ m/(transcript_id)/){  ## GTF group should have this attribute
			$ifGTF = 0;
			}
		}  ## else ends
	return($ifGTF);
	}  ## function ends
	###########################
1;
	
	## rules: 
		## length >= 9
		## 3,4,5: numeric
		## 3 < 4
		## 6: + or -
		## group: should have transcriptID

		


#0	1		2	3	4	5		6	7	8
#seq	source		feature	start	end	score		strand	frame	group
#chr10	hg18_knownGene	stop	82997	82999	0.000000	-	.	gene_id "uc001ifi.1"; transcript_id "uc001ifi.1"; 