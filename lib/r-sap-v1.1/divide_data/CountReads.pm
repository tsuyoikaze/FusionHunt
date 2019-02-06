package divide_data::CountReads;

use strict;
use warnings;

	sub count_reads{
	my($columnIndx,$fileName) = @_;
	my($str,@a1,@a2);
	my $totalReads = 0;

	open FHCountReads,$fileName or die "\n can not open file $fileName \n";

		#### skipping the top header in the input file #####
		while($str = <FHCountReads>){
			if($str =~ /^(\#)/){ next; }
			else{ last; }  ## else ends
		}  ## while(<FHDivideDataFileIn>) header skipping ends
		######################################################

	$str =~ s/\n//;
	$str =~ s/\r//;
	@a1  = split(/\t+/,$str);

		if(scalar(@a1) < 21){  ## if record length is short
		die "\n file format is not psl or pslx at line $. in file $fileName \n";
		}

	$totalReads += 1;
	
		while($str = <FHCountReads>){
		$str =~ s/\n//;
		$str =~ s/\r//;
		@a2 = split(/\t+/,$str);
	
			if(scalar(@a2) < 21){  ## if record length is short
			die "\n file format is not psl or pslx at line $. in file $fileName \n";
			}

			if($a1[$columnIndx] ne $a2[$columnIndx]){
			++$totalReads;
			@a1 = @a2;
			}
		}  ## while(<FHCountReads> ends
	
		if(scalar(@a2)){
	        	if($a1[$columnIndx] ne $a2[$columnIndx]){
	        	++$totalReads;
	        	}
	        }  ## if(!(scalar(@a2)) ends
	
	close FHCountReads;
	return($totalReads);
	}  ## function ends
	###########################
1;
