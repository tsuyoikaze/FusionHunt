package file_parsers::BEDToPslx;

use strict;
use warnings;

use file_parsers::CheckIfBED;


	sub bed_to_pslx{
	my($in,$out) = @_;


	open FHBedToPslxIn,$in or die "\n can not open file $in \n";
	open FHBedToPslxOut,">$out";  ## output file

	my($str,@a1,@a2,@temp1,@temp2);
	my($i,$len,$st);

		while($str = <FHBedToPslxIn>){
		$str =~ s/\n//;
		$str =~ s/\r//;
			if( ($str =~ m/^(\#+)/) || !($str) ){  ## skipping the header blank line
			next;
			}

		@a1 = split(/\t+/,$str);

			if(!(check_if_bed(@a1))){  ## if record is not in bed format
			die "\n record at line $. in file $in is not in bed format (terminating) \n";
			}

		$len = get_block_length($a1[10]);  ## function to be defined
		
		push(@a2,$len);  ## match
		push(@a2,(0,0,0,0,0,0,0));
		push(@a2,$a1[5]);  ## strand
		push(@a2,$a1[3]);  ## query name
		push(@a2,$len);  ## (chrEn-chrSt), query length
		push(@a2,(0,$len));  ## query start and end (0 based)
		push(@a2,($a1[0],'.'));  ## chromosome name and size (unknown)
		push(@a2,($a1[1],$a1[2]));  ## chrSt, chrEn
		push(@a2,($a1[9],$a1[10]));  ## block count, block sizes
				
		@temp1 = split(',',$a1[10]);  ## block sizes
	
		#### qSt coordinates ####	
		$st = 0;	
		push(@temp2,$st);
			for($i = 0;$i<scalar(@temp1-1);++$i){
			push(@temp2,($st+$temp1[$i]));
			}
		push(@a2,join(',',@temp2).',');
		splice(@temp2);
		##########################
		
		#### chrSt coordinates ####
		@temp1 = split(',',$a1[11]);  ## blSt
			for($i = 0;$i<scalar(@temp1);++$i){
			push(@temp2,($a1[1]+$temp1[$i]));
			}
		push(@a2,join(',',@temp2).',');
		splice(@temp2);
		###########################
		push(@a2,('.','.'));  ## qSeq and chrSeq
		print FHBedToPslxOut join("\t",@a2),"\n";
	
		splice(@a1);
		splice(@a2);
		}  ## while(<FHBedToPslxIn>) ends

	close FHBedToPslxIn;
	close FHBedToPslxOut;
	}  ## function ends
	###########################
	sub get_block_length{
	my($blSizes) = @_;
	my(@temp,$i,$l);
	@temp = split(',',$blSizes);
	$l = 0;
		for($i = 0;$i<scalar(@temp);++$i){
		$l += $temp[$i];
		}
	return($l);
	}  ## function ends
	###########################
1;	
	

##0	1	2	3			4			5	6	7	8	9	10	11
#chr	chrSt	chrEn	qId			score			strand	thickSt	thickEn	rgb	count	blSize	blSt
##chr10	3134164	3134909	chr10:3134164-3134909	2.6660591931269817	-	3134164	3134909	0,0,0	2	48,53,	0,692,


##248	0	0	0	0	0	0	0	+	##E3VCB0K02F0007	248	0	248	chr10	135374737	99467462	99467710	1	248,	0,	99467462,	gcattccagagctcactgcccttctagatgtgccttcccgcttggcttccagcggcttgtgctcactctgtctgccaggtatgagaagaacacgtaagaccgccaccacactcaccctccctcaaggccctgtgccataggggtggccacccgacctgcccccagaacttttggatactggaggcagttgcataggtctccctctctgggcaccaggactcagtccagcccaagactactctgggcag,	gcattccagagctcactgcccttctagatgtgccttcccgcttggcttccagcggcttgtgctcactctgtctgccaggtatgagaagaacacgtaagaccgccaccacactcaccctccctcaaggccctgtgccataggggtggccacccgacctgcccccagaacttttggatactggaggcagttgcataggtctccctctctgggcaccaggactcagtccagcccaagactactctgggcag,
##793	0	0	0	0	0	0	0	-	chr10:5825494-5833759	793	0	793	chr10	135374737	5825494		5833759		4	292,145,257,99,	0,292,,145,,257,	5825494,5827490,5832228,5833660,

