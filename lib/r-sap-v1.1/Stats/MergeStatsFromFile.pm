package Stats::MergeStatsFromFile;


use strict;
use warnings;


	sub merge_hashes_from_file{
	my($hashRef,$fin) = @_;
	
	open FHMergeHashes,$fin or die "\n can not open file $fin \n";

	my($str,@a1);
		while($str = <FHMergeHashes>){
		$str =~ s/\n//;
		$str =~ s/\r//;
	
		@a1 = split(/\t+/,$str);
		$a1[0] =~ s/\s//g;
		$a1[1] =~ s/\s//g;
			if(exists $hashRef->{$a1[0]}){  ## if key matches
			$hashRef->{$a1[0]} += $a1[1];
			}

	} ## while(<FHMergeHashes>) ends
	close FHMergeHashes;

	#open FHMergeHashes,">$fin";  ## deleting the original file (vales are now in hash)
	#close FHMergeHashes
	##shared::FileCheckPrint::print_hash_to_file($hashRef,$fin);
	}  ## function ends
1;
