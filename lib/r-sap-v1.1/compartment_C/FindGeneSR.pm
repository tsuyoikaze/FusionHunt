#################### Function file FindGeneSR.pm containing supporting sub routines for ##################
#################### FindGene.pl. ########################################################################
package compartment_C::FindGeneSR;


use strict;
use warnings;


	sub make_gene_track_index{
	my($chrNum,$fin) = @_;

	my %ref = ();
	my ($str,$index1,$index2);
	my(@a1,@a2);
	
	open FH21,$fin or die "\n can not open file input file $fin \n";  ## gene tracsk downloaded from UCSC
	$index1 = tell FH21;
	$str = <FH21>;
	
		if($str =~ m/\#/){
		$index1 = tell FH21;
		$str = <FH21>;
		}
	
	$str =~ s/\n//;
	$str =~ s/\r//;
	@a1 = split(/\t+/,$str);
	$index2 = tell FH21;
	
		while($str = <FH21>){
		$str =~ s/\n//;
		$str =~ s/\r//;
		@a2 = split(/\t+/,$str);
			if($a1[$chrNum] ne $a2[$chrNum]){
			$ref{$a1[$chrNum]} = $index1;
			@a1 = @a2;
			splice(@a2);
			$index1 = $index2;
			}
		$index2 = tell FH21;
		}  ## while(<FH21>) ends 
		
	$ref{$a1[$chrNum]} = $index1;
	close FH21;
	return(%ref);
	}  ## function ends 2
	###########################
	sub annotate_out_files{
	my($fileTag)= @_;
	my %returnHash = (
		"ExonOnly" => $fileTag."ExonOnly.out",
		"ExonDeletion" => $fileTag."ExonDeletion.out",
		"ExonExtension" => $fileTag."ExonExtension.out",
		"IntronOnly" => $fileTag."IntronOnly.out",
		"GeneExpansion" => $fileTag."GeneExpansion.out",
		"NeighboringExon" => $fileTag."NeighboringExon.out",
		"MultipleAnnotations" => $fileTag."MultipleAnnotations.out",
		"GeneDesert" => $fileTag."GeneDesert.out",
		);
	create_file_from_hash(\%returnHash);
	return(%returnHash);
	}  ## function ends
	###########################
	sub create_file_from_hash{
	my($hashRef) = @_;
	my $key1;
		for $key1(keys %{$hashRef}){
		open FHCreateOutFile,">".$hashRef->{$key1};
		close FHCreateOutFile;
		}  ## foreach(keys(%{$hashRef})) ends
	}  ## function ends
	###########################
	sub find_gene_overlap{
	my($st1Ref,$en1Ref,$st2,$en2) = @_;
	my($overlap,$i);
	$overlap = 0;
		for($i = 0;$i<scalar(@{$st1Ref});++$i){  ## for each aligned block in the pslx record
		$overlap += GeneralPack::intersect_track_alignment($st1Ref->[$i],$en1Ref->[$i],$st2,$en2);
		}
	return($overlap);
	}  ## function ends
	###########################
	sub find_exon_deletions{
	my($st1Ref,$en1Ref,$st2Ref,$en2Ref,$deletionCutoff) = @_;
	my ($i,$j,$alignmentDeletion,$fullExonSkipping);
	my $exonDeletion = 0;
		for($i = 0;$i<(scalar(@{$st1Ref})-1);++$i){
			for($j = 0;$j<scalar(@{$st2Ref});++$j){
			$alignmentDeletion = GeneralPack::intersect_track_alignment($en1Ref->[$i],$st1Ref->[$i+1],$st2Ref->[$j],$en2Ref->[$j]);
			$fullExonSkipping = $alignmentDeletion/($en2Ref->[$j]-$st2Ref->[$j]);
				if( ($alignmentDeletion > $deletionCutoff) || ($fullExonSkipping == 1) ){
				$exonDeletion += $fullExonSkipping;
				}
			}  ## for($j) ends
		}  ## for($i) ends
	return($exonDeletion);
	}  ## function ends
	###########################
	sub find_gene_expansion_flags{
	my($st1,$en1,$txSt,$txEn,$strand) = @_;
	my($geneExpansionFlag,$end5ExpansionFlag,$end3ExpansionFlag);
	my($expansion);
	$geneExpansionFlag = 0;
	$end5ExpansionFlag = 0;
	$end3ExpansionFlag = 0;
	$expansion = 0;
		if(GeneralPack::intersect_track_alignment($st1,$en1,$txSt,$txEn)){
			if($strand eq '+'){
				if($st1 < $txSt){
				$geneExpansionFlag += 1;
				$expansion += ($txSt-$st1);	
				$end5ExpansionFlag = 1;
				}
				if($en1 > $txEn){
				$geneExpansionFlag += 1;
				$expansion += ($en1-$txEn);
				$end3ExpansionFlag = 1;
				}
			}  ## if(strand eq '+') ends
			elsif($strand eq '-'){
				if($en1 > $txEn){
				$geneExpansionFlag += 1;
				$expansion += ($en1-$txEn);
				$end5ExpansionFlag = 1;
				}
				if($st1 < $txSt){
				$geneExpansionFlag += 1;
				$expansion += ($txSt-$st1);
				$end3ExpansionFlag = 1;
				}
			}  ## elsif($strand eq '-') ends
		}  ## if there is any intersect
	return($geneExpansionFlag,$expansion);
	}  ## function ends
	###########################
	sub define_char_array_field_hash{
	my($hashRef) = @_;
	$hashRef->{"clustIndex"} = 0; 
	$hashRef->{"transcriptID"} = 1; 
	$hashRef->{"exonOverlap"} = 2; 
	$hashRef->{"geneOverlap"} = 3; 
	$hashRef->{"exonDeletions"} = 4; 
	$hashRef->{"intronOverlap"} = 5; 
	$hashRef->{"geneExpansionFlag"} = 6;
	$hashRef->{"geneDistance"} = 7;  ## for neighboring genes (within the permissible distance)
	$hashRef->{"codingFlag"} = 8;
	}  ## function ends
	###########################
	sub sort_char_array{  ## sorting intersected transcripts on different characterstics
	my($hashRef,$clustRef) = @_;
	my($i,$j,@temp);	
	## note: all of the transcripts have same exon overlap

		for($i = 0;$i<(scalar(@{$clustRef})-1);++$i){
			for($j = 0;$j<(scalar(@{$clustRef})-1-$i);++$j){
				if($clustRef->[$j+1][$hashRef->{"exonOverlap"}] > $clustRef->[$j][$hashRef->{"exonOverlap"}]){  ## comparing exon overlap
				GeneralPack::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
				}  ## if(ends)
				elsif($clustRef->[$j+1][$hashRef->{"exonOverlap"}] == $clustRef->[$j][$hashRef->{"exonOverlap"}]){  ## if overlaps with the gene boundaries are equal
					if($clustRef->[$j+1][$hashRef->{"geneOverlap"}] > $clustRef->[$j][$hashRef->{"geneOverlap"}]){
					GeneralPack::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
					}
					elsif($clustRef->[$j+1][$hashRef->{"geneOverlap"}] == $clustRef->[$j][$hashRef->{"geneOverlap"}]){  ## if overlaps with the gene boundaries are equal
						if($clustRef->[$j+1][$hashRef->{"exonDeletions"}] < $clustRef->[$j][$hashRef->{"exonDeletions"}]){
						GeneralPack::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
						}
						elsif($clustRef->[$j+1][$hashRef->{"exonDeletions"}] == $clustRef->[$j][$hashRef->{"exonDeletions"}]){  ## if overlaps with the gene boundaries are equal
							if($clustRef->[$j+1][$hashRef->{"intronOverlap"}] < $clustRef->[$j][$hashRef->{"intronOverlap"}]){
							GeneralPack::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
							}
							elsif($clustRef->[$j+1][$hashRef->{"intronOverlap"}] == $clustRef->[$j][$hashRef->{"intronOverlap"}]){  ## if overlaps with the gene boundaries are equal
								if($clustRef->[$j+1][$hashRef->{"geneExpansionFlag"}] < $clustRef->[$j][$hashRef->{"geneExpansionFlag"}]){
								GeneralPack::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
								}
								elsif($clustRef->[$j+1][$hashRef->{"geneExpansionFlag"}] == $clustRef->[$j][$hashRef->{"geneExpansionFlag"}]){  ## if overlaps with the gene boundaries are equal
									if($clustRef->[$j+1][$hashRef->{"geneDistance"}] < $clustRef->[$j+1][$hashRef->{"geneDistance"}]){
									GeneralPack::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
									}
									elsif($clustRef->[$j+1][$hashRef->{"geneDistance"}] == $clustRef->[$j+1][$hashRef->{"geneDistance"}]){  ## if distance from gene is equal
										if($clustRef->[$j+1][$hashRef->{"codingFlag"}] > $clustRef->[$j][$hashRef->{"codingFlag"}]){
										GeneralPack::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
										}
									}  ## elsif() ends
								}  ## elsif() ends
							}  ## elsif() ends					
						}  ## elsif() ends
					}  ## elsif() ends
				}  ## elsif() ends
			}  ## for($j) ends
		}  ## for($i) ends
	}  ## function ends
	###########################
	sub find_top_replicates{  ## function to find out all the transcripts with same alignment characteristics as the top transcript
	my($hashRef,$clustRef) = @_;
	my($i,$j,$matchFlag,$topReplicates);
	$topReplicates = 0;
		for($i = 0;$i<scalar(@{$clustRef});++$i){
		$matchFlag = 0;
			for($j = 2;$j<scalar(keys(%{$hashRef}));++$j){
				if($clustRef->[$i][$j] == $clustRef->[0][$j]){
				$matchFlag += 1;
				}
			}  ## for(scalar(keys(%{$hashRef}))) ends
			if($matchFlag == (keys(%{$hashRef})-2)){
			$topReplicates += 1;
			}
		}  ## for(scalar(@{$clustRef})) ends
	
		if(!($topReplicates)){
		die "\n no top replicate found for $clustRef->[0][1]\n";
		}
	return($topReplicates);
	}  ## function ends
	###########################
	sub sort_top_gene_ids{
	my($n,$hashRef,$clustRef) = @_;
	my($sortIndex,$i,$j);
	$sortIndex = $hashRef->{"transcriptID"};
		for($i = 0;$i<($n-1);++$i){
			for($j = 0;$j<($n-1-$i);++$j){
				if($clustRef->[$j+1][$sortIndex] lt $clustRef->[$j][$sortIndex]){  ## comparing exon overlap
				GeneralPack::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
				}  ## if(ends)
			}  ## for($j) ends
		}  ## for($i) ends
	}  ## function ends
	###########################
	sub decide_annotation_tag{
	my($hashRef,$clustRef,$currTrackIndx,$annotTagsRef) = @_;
	
		if($clustRef->[$currTrackIndx][$hashRef->{"exonOverlap"}]){  ## if there is any exon overlap
			if( !($clustRef->[$currTrackIndx][$hashRef->{"intronOverlap"}]) && !($clustRef->[$currTrackIndx][$hashRef->{"geneExpansionFlag"}]) ){  ## if neither intron overlap nor gene boundart expansion
				if($clustRef->[$currTrackIndx][$hashRef->{"exonDeletions"}]){
				push(@{$annotTagsRef},"ExonDeletion");
				}
				else{
				push(@{$annotTagsRef},"ExonOnly");
				}
			}
			if($clustRef->[$currTrackIndx][$hashRef->{"intronOverlap"}]){
			push(@{$annotTagsRef},"ExonExtension");
			}
			if($clustRef->[$currTrackIndx][$hashRef->{"geneExpansionFlag"}]){
			push(@{$annotTagsRef},"GeneExpansion");
			}
		}  ## if(exonOverlap) ends
		
		if(!($clustRef->[$currTrackIndx][$hashRef->{"exonOverlap"}])){
			if($clustRef->[$currTrackIndx][$hashRef->{"intronOverlap"}]){
			push(@{$annotTagsRef},"IntronOnly");
			}
			if($clustRef->[$currTrackIndx][$hashRef->{"geneExpansionFlag"}]){
			push(@{$annotTagsRef},"GeneExpansion");
			}
		}  ## if(!ExonOverlap) ends
		if(!(scalar(@{$annotTagsRef}))){
		push(@{$annotTagsRef},"NeighboringExon");
		}
	}  ## function ends
	###########################
	sub check_if_exon_skipping{
	my ($st1Ref,$en1Ref,$st2Ref,$en2Ref) = @_; ## obtaining references to the arrays in the main function
	my($skippingFlag,$skippedExons,$skippedExonsCount);
	my(%exonNums,$exonFlag);
	my($i,$j,@temp);

		for($i = 0;$i<scalar(@{$st1Ref});++$i){  ## for each aligned block in the pslx record
		$exonFlag = 0;
			for($j = 0;$j<scalar(@{$st2Ref});++$j){  ## for each exon
				if(GeneralPack::intersect_track_alignment($st1Ref->[$i],$en1Ref->[$i],
					$st2Ref->[$j],$en2Ref->[$j])){
				$exonNums{$i} = $j;
				$exonFlag = 1;
				}
			}  ## for(exons) ends
		}  ## for(alignment blocks) ends
	$skippingFlag = 0;
	$skippedExons = 0;

		if($exonFlag){
			for($i = 1;$i<scalar(@{$st1Ref});++$i){
				if(scalar(@{$st1Ref}) == keys(%exonNums)){
					if( ($exonNums{$i} - $exonNums{$i-1}) > 1 ){
					$skippingFlag = 1;
						for($j = ($exonNums{$i-1}+1);$j<=($exonNums{$i}-1);++$j){
						$skippedExons .= $j.",";  ## add in skipped exon numbers
						}
					} ## if ends
				}  ## if($nr == keys(%exonNums)) ends
			}  ## for($i) ends
		}  ## if($exonFlag) ends
		if($skippedExons){
		@temp = split(',',$skippedExons);
		$skippedExonsCount = scalar(@temp);
		}
		else{
		$skippedExonsCount = 0;
		}
	return($skippedExonsCount);  ## retutning number of exons skipped
	}  ## function ends
	###########################
	sub initialize_array{
	my($arrRef,$hashRef) = @_;
	my $i;
		for($i = 0;$i<scalar(keys %{$hashRef});++$i){
		$arrRef->[$i] = 0;
		}  ## for($i) ends
	}  ## function ends
	###########################


1;


