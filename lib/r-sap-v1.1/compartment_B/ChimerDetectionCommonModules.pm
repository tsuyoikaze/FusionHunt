package compartment_B::ChimerDetectionCommonModules;

use strict;
use warnings;

	sub check_top_hit_multi{  ## Common
	my($repeatIdentityCutoff,$topMultiFlagRef,$clustRef,$pslx2IndxRef) = @_;
	my($i,$j);
	${$topMultiFlagRef} = 0;
	
		####################### Repeat detection ######################
		for($i = 1;$i<scalar(@{$clustRef});++$i){  ## for all the hits in the cluster
			if( (shared::CoordinateComparisons::intersect_track_alignment($clustRef->[0][$pslx2IndxRef->{"qSt"}],$clustRef->[0][$pslx2IndxRef->{"qEn"}],
						$clustRef->[$i][$pslx2IndxRef->{"qSt"}],$clustRef->[$i][$pslx2IndxRef->{"qEn"}])/($clustRef->[0][$pslx2IndxRef->{"qEn"}]-$clustRef->[0][$pslx2IndxRef->{"qSt"}])) >= 0.9){
				if($clustRef->[$i][$pslx2IndxRef->{"identity"}] >= $repeatIdentityCutoff){
				++${$topMultiFlagRef};
				}  ## if(1) ends
			}  ## if(2) ends
		}  ## for($i) ends
		######################### Repeat detection ####################

		if(${$topMultiFlagRef}){  ## if there is any similar multiple hits
		++${$topMultiFlagRef};  ## increment for the top one
		} ## if ends
	return(${$topMultiFlagRef});
	}  ## function ends
	###########################
	sub screen_alignments{
	my($identityCutoff,$repeatIdentityCutoff,$pslx2IndxRef,$clustRef,$covClustRef,$qualifiedAlignmentsRef) = @_;
	my(@a1,$alignmentScreenTag);

	@a1 = @{$clustRef->[0]}[0..(scalar(keys(%{$pslx2IndxRef}))-1)];  ## top hit
	
	my ($i,$j,$c,$n,$id,$same,$same_index);
	my($sameFlag,@temp);
	$sameFlag = 0;

		for($i = 0;$i<scalar(@{$covClustRef});++$i){    ## sorting array_cover based on coverage 
			for($j = 0;$j<(scalar(@{$covClustRef})-1-$i);++$j){
				if($covClustRef->[$j][(scalar(keys(%{$pslx2IndxRef}))-1)+1] < $covClustRef->[$j+1][(scalar(keys(%{$pslx2IndxRef}))-1)+1]){ ## sorting on coverage
				shared::ArrayOperations::swap_elements($j,($j+1),(scalar(keys(%{$pslx2IndxRef}))+1),$covClustRef);
				}
				elsif($covClustRef->[$j][(scalar(keys(%{$pslx2IndxRef}))-1)+1] == $covClustRef->[$j+1][(scalar(keys(%{$pslx2IndxRef}))-1)+1]){  ## if coverage is same
					if($covClustRef->[$j][$pslx2IndxRef->{"identity"}] < $covClustRef->[$j+1][$pslx2IndxRef->{"identity"}]){  ## soring on identity
					shared::ArrayOperations::swap_elements($j,($j+1),(scalar(keys(%{$pslx2IndxRef}))+1),$covClustRef);
					}
				}  ## elsif ends
			}  ## for(j) ends
		  }  ## for(i) ends, sorting ends
						
		$c = $covClustRef->[0][(scalar(keys(%{$pslx2IndxRef}))-1)+1];  ## highest coverage
		$id = $covClustRef->[0][$pslx2IndxRef->{"identity"}];  ## identity of that read
		$n = 0;
		$same = 0;
	
			#### multi-hits and same chromosome (intra-chromosomal)hits ####
			for($i = 0;$i<scalar(@{$covClustRef});++$i){
				if( ($covClustRef->[$i][(scalar(keys(%{$pslx2IndxRef}))-1)+1] >= (0.9*$c)) && ($covClustRef->[$i][$pslx2IndxRef->{"identity"}] >= $repeatIdentityCutoff) ){
				++$n;  ## incrementing multi-hits
					if($a1[$pslx2IndxRef->{"tID"}] eq $covClustRef->[$i][$pslx2IndxRef->{"tID"}]){
					++$same;
					}  ## if($chr1 eq chr2) ends
				}  ## if(1) ends
			}  ## for loop end


		#### outfile print conditions ####
		if($same == 0){  ## if no intra-chromosomal (intra-chromosomal gets preferecnce over inter-chromosomal)
			if($n > 1){
				for($i = 0;$i<scalar(@{$covClustRef});++$i){
					if(($covClustRef->[$i][(scalar(keys%{$pslx2IndxRef})-1)+1] >= 0.9*$c) && 
						($covClustRef->[$i][$pslx2IndxRef->{"identity"}] >= $repeatIdentityCutoff)){
					$alignmentScreenTag = "MultiHits";
					push(@{$qualifiedAlignmentsRef},$i);
					}
				}  ## for loop end
			}  ## if($n >1)
			
			elsif($n == 1){  ## if only one hit
				if($covClustRef->[0][$pslx2IndxRef->{"identity"}] >= $identityCutoff){
				$alignmentScreenTag = "Recombination";
				push(@{$qualifiedAlignmentsRef},0);
				}  ## if($a1[0][identity] >= $identityCutoff) ends     

				elsif($covClustRef->[0][$pslx2IndxRef->{"identity"}] < $identityCutoff){
				$alignmentScreenTag = "NoRecombination";
				}  ## elsif($a2[0][identity] < $identityCutoff) ends
			}  ## if($n == 1) ends
		}  ## if($same == 0) ends
		
		elsif($same == 1){
			for($i = 0;$i<scalar(@{$covClustRef});++$i){
				if(($covClustRef->[$i][(scalar(keys%{$pslx2IndxRef})-1)+1] >= (0.9*$c)) &&
					($covClustRef->[$i][$pslx2IndxRef->{"identity"}] >= $identityCutoff) && ($a1[$pslx2IndxRef->{"tID"}] eq $covClustRef->[$i][$pslx2IndxRef->{"tID"}])){

				$sameFlag = 1;
				$alignmentScreenTag = "SameRecombination";
				push(@{$qualifiedAlignmentsRef},$i);
				last;
				}  ## if ends
			}  ## for loop ends
		
			if(!($sameFlag)){  ## if same chr hit didn't fulfill the criteria
				for($i = 0;$i<scalar(@{$covClustRef});++$i){
					if(($covClustRef->[$i][(scalar(keys%{$pslx2IndxRef})-1)+1] >= (0.9*$c)) && 
						($covClustRef->[$i][$pslx2IndxRef->{"identity"}] >= $repeatIdentityCutoff) && ($a1[$pslx2IndxRef->{"tID"}] eq $covClustRef->[$i][$pslx2IndxRef->{"tID"}])){

					$alignmentScreenTag = "RejectedAlignments";
					push(@{$qualifiedAlignmentsRef},$i);
					}  ## if ends
				}  ## for loop ends
			}  ## if(!($sameFlag)) ends
		}  ## elsif($same == 1) ends

		elsif($same > 1){
			for($i = 0;$i<scalar(@{$covClustRef});++$i){
				if(($covClustRef->[$i][(scalar(keys%{$pslx2IndxRef})-1)+1] >= (0.9*$c)) && 
					($covClustRef->[$i][$pslx2IndxRef->{"identity"}] >= $repeatIdentityCutoff) && ($a1[$pslx2IndxRef->{"tID"}] eq $covClustRef->[$i][$pslx2IndxRef->{"tID"}])){
				$alignmentScreenTag = "MultiHits";
				push(@{$qualifiedAlignmentsRef},$i);
				}				
			}  ## for loop ends
		}  ## elsif($same == 1) ends
	
		# temporary warning  ## to be removed later  ##
		if(!($alignmentScreenTag)){  ## if no alignment tag was defined
		die "\n no alignment tag was defined for ",$clustRef->[0][$pslx2IndxRef->{"qID"}]," at line $. \n";
		}
		#################################################

	
	
	return($alignmentScreenTag);
	}  ## function ends
	###########################
1;




