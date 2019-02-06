#### perl script to conver Cufflinks generated GTF file for assembled and mapped transcript to format similar to pslx that can be used by the pipelie ####
#### Note: Cufflink does not report '0' based start coordinates.
#### Script will look only for exon features in the GTF records sinice Cufflinks feature list contains only transcript and exon.


package file_parsers::GTFToPslx;


use strict;
use warnings;

use file_parsers::CheckIfGTF;
use file_parsers::GTFToDefault;


	sub gtf_to_pslx{
	my($inFile,$outFile) = @_;

	open FHGTFToPslxIn,$inFile or die "\n can not open file $inFile \n";
	open FHGTFToPslxOut,">$outFile";  ## output file


	my($str,$l,$i,$j,$n);
	my($pt1,$pt2);

	my(@a1,@a2,@st,@en);  
	## GTF record line, pslx record, alignment blocks starts, alignment block ends

	my($chr,$strand,$tId);  ## chromosome, strand and transcript ID for the GTF records
	my($tempId,@temp,$tempChr); ## Id of current transcript , temporary array


	$n = 0;
	$pt1 = tell FHGTFToPslxIn;
	
		while($str=<FHGTFToPslxIn>){
		$pt2 = tell FHGTFToPslxIn;
		$str =~ s/\n//;
		$str =~ s/\r//;

			if( ($str =~ m/^(\#+)/) || !($str) ){  ## skipping the header and blank line
			next;
			$pt1 = $pt2
			}

		@a1 = split(/\t+/,$str);

			if(!(file_parsers::CheckIfGTF::check_if_GTF(\@a1))){
			die "\n record at line $. is not in GTF format in file $inFile (terminating) \n";
			}

			if($a1[2] =~ m/(exon)/){  ## if feature is transcript
				if(!($n)){  ## if first record
				$chr = $a1[0];
				$strand = $a1[6];

				@temp = file_parsers::GTFToDefault::separate_goup_elements($a1[8]);
				$tId = $temp[0];
				push(@st,($a1[3]-1));  ## coordinates are 0 based now
				push(@en,$a1[4]);
				$pt1 = $pt2;
				++$n;
					if(eof){  ## if file ends (only one transcript feature)
					@a2 = generate_pslx($tId,$chr,$strand,\@st,\@en);

						if(scalar(@a2)){
						print FHGTFToPslxOut join("\t",@a2),"\n";
						}
					splice(@a2); splice(@st); splice(@en);
					}  ## if(eof) ends
				}  ## if($n == 0) ends

				else{  ## if not first gtf entry
				@temp = file_parsers::GTFToDefault::separate_goup_elements($a1[8]);
				$tempId = $temp[0];
				$tempChr = $a1[0];

					if($tId eq $tempId){  ## for the same transcript Id
					push(@st,($a1[3]-1));  ## coordinates are NOT 0 based currently
					push(@en,$a1[4]);
					$pt1 = $pt2;
					}

					if( (eof) || ($tId ne $tempId) || ($chr ne $tempChr) ){  ## change of transcriptID
					@a2 = generate_pslx($tId,$chr,$strand,\@st,\@en);

						if(scalar(@a2)){
						print FHGTFToPslxOut join("\t",@a2),"\n";
						}
					splice(@a2); splice(@st); splice(@en);
					$n = 0;
					seek(FHGTFToPslxIn,$pt1,0);  ## resetting the file pointer
					}  ## if(record change ends)
				}  ## else(if($n)) ends
			}  ## if(exon) ends

			else{
			next;
			}
		}  ## while(<FHGTFToPslxIn>) ends

	close FHGTFToPslxIn;
	close FHGTFToPslxOut;
	}  ## function ends
##############################
	sub generate_pslx{
	my($tId,$chr,$strand,$stRef,$enRef) = @_;
	my($l,@pslx);
	splice(@pslx);
	
		$l = get_transcript_length($stRef,$enRef);
		push(@pslx,($l,0,0,0,0,0));  ## matches, mismatches, rm, nc,qInserts,qInsertBases

		push(@pslx,(scalar(@{$stRef})-1));  ## tInserts = number of exons - 1

		push(@pslx,($enRef->[scalar(@{$enRef})-1] - $stRef->[0] - $l));  ## intron length

		push(@pslx,$strand);
		push(@pslx,($tId,$l));  ## query name, query length
		push(@pslx,(0,$l));  ## qSt, qEn
		push(@pslx,($chr,'.'));  ## chromosome name, chromosome size (unknown)
		push(@pslx,(($stRef->[0]),$enRef->[scalar(@{$enRef})-1]));  ## chrSt, chrEn
		push(@pslx,scalar(@{$stRef}));  ## number of aligned blocks = number of exons

		push(@pslx,join(',',get_block_lengths($stRef,$enRef)).",");

		push(@pslx,join(',',get_query_starts($strand,$stRef,$enRef)).",");
		##$stRef->[0] -= 1; ## (converging to 0 base)
		push(@pslx,join(',',@{$stRef}).",");
		##push(@pslx,('.','.'));  ## qSeq, tSeq (to make pslx record)

	return(@pslx);
	}  ## function ends
	###########################
	sub get_transcript_length{
	my($stRef,$enRef) = @_;
	my($i,$len);
	$len = 0;
		for($i = 0;$i<scalar(@{$stRef});++$i){
		$len += ($enRef->[$i]-$stRef->[$i]);
		}
	return($len);
	}  # function ends
	###########################
	sub get_block_lengths{
	my($stRef,$enRef) = @_;
	my($i,@len);
		for($i = 0;$i<scalar(@{$stRef});++$i){
		push(@len,($enRef->[$i]-$stRef->[$i]));
		}
	return(@len);
	}  # function ends
	###########################
	sub get_query_starts{
	my($strand,$stRef,$enRef) = @_;
	my($i,@ens,@sts,$len);
	$len = 0;
	
		for($i = 0;$i<scalar(@{$stRef});++$i){
		push(@sts,$len);
		$len += ($enRef->[$i]-$stRef->[$i]);
		push(@ens,$len);
		}
		if($strand eq '-'){
		splice(@sts);
			for($i = (scalar(@ens)-1);$i>=0;--$i){
			push(@sts,($len-$ens[$i]));
			}
		}  ## if(strand '-') ends
	return(@sts);
	}  ## function ends
	###########################

1;	

