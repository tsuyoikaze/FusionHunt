package divide_data::DivideDataFile2;

use strict;
use warnings;


	sub divide_data_file2{
	my($inFile,$pslxIndxRef,$N,$chunkStRef,$chunkEnRef) = @_;  ## directory, input file, query file hash, number of parts to be created
	
	open FHDivideDataFileIn,$inFile;

	my($endPoint,$chunkSize,$prevPoint);
	my($i,$str1);
	my(@a1,@a2);
	my($fp0,$fp1,$fp2);  ## variables to keep track of the file pointers around the current location
	my(@chunkSt,@chunkEn);

	$fp1 = tell FHDivideDataFileIn;

		#### skipping the top header in the input file #####
		while($str1 = <FHDivideDataFileIn>){
		$fp2 = tell FHDivideDataFileIn;
			if($str1 =~ /^(\#)/){  $fp1 = $fp2; next; }
			else{ last; }  ## else ends
		}  ## while(<FHDivideDataFileIn>) header skipping ends
		######################################################


	seek(FHDivideDataFileIn,0,2);  ## setting file pointer to the end of the file
	$endPoint = tell FHDivideDataFileIn;

	seek(FHDivideDataFileIn,$fp1,0);  ## going to start of the file


	$chunkSize = int($endPoint/$N);  ## not exact (estimated)

	$prevPoint = $fp1;

		for($i = 0;$i<$N;++$i){
			if(($prevPoint+$chunkSize) >= $endPoint){
			last;
			}

		seek(FHDivideDataFileIn,($prevPoint+$chunkSize),0); ## seetting the pointer to a new ends	

		$str1 = <FHDivideDataFileIn>;  ## reading in the partial line
			if((tell FHDivideDataFileIn) >= $endPoint){
			last;
			}
		$fp0 = tell FHDivideDataFileIn;
		$str1 = <FHDivideDataFileIn>;
		$fp1 = tell FHDivideDataFileIn;
		@a1 = split(/\t+/,$str1);

			while($str1 = <FHDivideDataFileIn>){
			$fp2 = tell FHDivideDataFileIn;
			@a2 = split(/\t+/,$str1);
				if(eof){
				push(@chunkSt,$prevPoint);
				push(@chunkEn,$endPoint);
				last;
				}

				if($a1[$pslxIndxRef->{"qID"}] ne $a2[$pslxIndxRef->{"qID"}]){  ## if no id match
				push(@chunkSt,$prevPoint);
				push(@chunkEn,$fp1);
				$prevPoint = $fp1;
				last;
				}
				else{
				$fp1 = $fp2;
				}
			}  ## while(<FHDivideDataFileIn>) ends
		}  ## for($parts) ends
	
	close FHDivideDataFileIn;

	push(@chunkSt,$prevPoint);
	push(@chunkEn,$endPoint);

	push(@{$chunkStRef},$chunkSt[0]);
	push(@{$chunkEnRef},$chunkEn[0]);

		for($i = 1;$i<scalar(@chunkSt);++$i){
			if( ($chunkSt[$i] ne $chunkStRef->[scalar(@{$chunkStRef})-1]) && ($chunkEn[$i] ne $chunkEnRef->[scalar(@{$chunkEnRef})-1]) ){
			push(@{$chunkStRef},$chunkSt[$i]);
			push(@{$chunkEnRef},$chunkEn[$i]);
			}
		}  ## for(@chunkSt-1) ends
	splice(@chunkSt);
	splice(@chunkEn);
	}  ## function ends
	###########################
1;
