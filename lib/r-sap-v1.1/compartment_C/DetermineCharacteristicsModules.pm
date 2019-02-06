
package compartment_C::DetermineCharacteristicsModules;


use strict;
use warnings;

require compartment_C::Pslx2RecordProcessing;

	sub check_if_large_intron{
	my($a1Ref,$pslx2IndxRef,$cutOff) = @_;
	my($i,@exSt,@exEn);
	my $largeIntronFlag = 0;  ## by default false (alignment is OK)
	compartment_C::Pslx2RecordProcessing::get_pslx2_algn_blocks($pslx2IndxRef,$a1Ref,\@exSt,\@exEn);
		for($i = 0;$i<(scalar(@exSt)-1);++$i){
			if( ($exSt[$i+1]-$exEn[$i]) > $cutOff){
			++$largeIntronFlag;
			}
		}  ## for(scalar(@exSt)) ends
	return($largeIntronFlag);
	}  ## function ends
	###########################
	sub find_gene_overlap{
	my($st1Ref,$en1Ref,$st2,$en2) = @_;
	my($overlap,$i);
	$overlap = 0;
		for($i = 0;$i<scalar(@{$st1Ref});++$i){  ## for each aligned block in the pslx record
		$overlap += shared::CoordinateComparisons::intersect_track_alignment($st1Ref->[$i],$en1Ref->[$i],$st2,$en2);
		}
	return($overlap);
	}  ## function ends
	###########################
	sub find_exon_deletions{
	my($st1Ref,$en1Ref,$st2Ref,$en2Ref,$deletionCutoff) = @_;
	my($i,$j,$alignmentDeletion,$fullExonSkipping);
	my $exonDeletion = 0;
		for($i = 0;$i<(scalar(@{$st1Ref})-1);++$i){
			for($j = 0;$j<scalar(@{$st2Ref});++$j){
			$alignmentDeletion = shared::CoordinateComparisons::intersect_track_alignment($en1Ref->[$i],$st1Ref->[$i+1],$st2Ref->[$j],$en2Ref->[$j]);
			$fullExonSkipping = $alignmentDeletion/($en2Ref->[$j]-$st2Ref->[$j]);
				if( ($alignmentDeletion > $deletionCutoff) || ($fullExonSkipping == 1) ){
				$exonDeletion += $fullExonSkipping;
				}
			}  ## for($j) ends
		}  ## for($i) ends
	return($exonDeletion);
	}  ## function ends
	###########################
	sub find_intron_overlaps{
	my($st1Ref,$en1Ref,$st2Ref,$en2Ref,$cutOff) = @_;
	my($i,$j,$overlap,$totalIntronOverlap,$intronOverlapFlag);
	$totalIntronOverlap = 0;
	$intronOverlapFlag = 0;

		for($i = 0;$i<scalar(@{$st1Ref});++$i){
			for($j = 0;$j<scalar(@{$st2Ref});++$j){
			$overlap = shared::CoordinateComparisons::intersect_track_alignment($st1Ref->[$i],$en1Ref->[$i],$st2Ref->[$j],$en2Ref->[$j]);
			$totalIntronOverlap += $overlap;
				if($overlap > $cutOff){
				$intronOverlapFlag += 1;
				}
			}  ## for($j) ends
		}  ## for($i) ends
	return($totalIntronOverlap,$intronOverlapFlag);
	}  ## function ends
	###########################
	sub find_gene_expansion_flags{
	my($st1,$en1,$txSt,$txEn,$cutOff) = @_;
	my($trueGeneExpansionFlag,$totalGeneExpansion,$expansion);
	
	$trueGeneExpansionFlag = 0;
	$totalGeneExpansion = 0;
	$expansion = 0;

		if(shared::CoordinateComparisons::intersect_track_alignment($st1,$en1,$txSt,$txEn)){
			if($st1 < $txSt){
			$expansion = ($txSt-$st1);
			$totalGeneExpansion += $expansion;
				if($expansion > $cutOff){
				$trueGeneExpansionFlag += 1;
				}
			}  ## if($st < $txSt)
			if($en1 > $txEn){
			$expansion = ($en1-$txEn);
			$totalGeneExpansion += $expansion;
				if($expansion > $cutOff){
				$trueGeneExpansionFlag += 1;
				}
			}  ## if($en1 > $txEn)
		}  ## if there is any intersect
	return($totalGeneExpansion,$trueGeneExpansionFlag);
	}  ## function ends
	###########################
	sub find_intron_chain_match{
	my($st1Ref,$en1Ref,$st2Ref,$en2Ref) = @_;
	my($intronChainMatch,$i,$j);
	$intronChainMatch = 0;

		for($i = 0;$i<scalar(@{$st1Ref});++$i){  ## for each aligned block in the pslx record
			for($j = 0;$j<scalar(@{$st2Ref});++$j){
				if($st1Ref->[$i] == $st2Ref->[$j]){
				$intronChainMatch += 1;
				}
				if($en1Ref->[$i] == $en2Ref->[$j]){
				$intronChainMatch += 1;
				}
			}  ## for($j) ends
		}  ## for($i) ends
	return($intronChainMatch);
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
				if(shared::CoordinateComparisons::intersect_track_alignment($st1Ref->[$i],$en1Ref->[$i],
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
1;