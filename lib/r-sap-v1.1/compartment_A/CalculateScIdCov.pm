############## Package to calculate score, %coverage and %identity for a given BLAT pslx output file. ################
############## %Identity and original psl record. pslx file must be verified before using this package ###############
############## make sure that only numerical values are being used and file header is removed. Three #################
############## new columns score, %coverage and %identity respectively will be added to the input array. ##############

package compartment_A::CalculateScIdCov;

use strict;
use warnings;

require shared::CheckStringsNumericals;

	sub main_module_CalculateScIdCov{
	my($pslxIndxRef,$clustRef)  = @_;
	my($i,@pslx2Clust,@temp);	
	my ($id,$score,$coverage);  ## temporary
		for($i = 0;$i<scalar(@{$clustRef});++$i){
		$id = 100.00 - (calculate_Id($pslxIndxRef,$i,$clustRef))*0.1;
		$score = calculate_score($pslxIndxRef,$i,$clustRef);
		$coverage = calculate_coverage($pslxIndxRef,$i,$clustRef);
		push(@temp,($score,sprintf("%.2f",$coverage),sprintf("%.2f",$id)));
		push(@temp,@{$clustRef->[$i]}[0..(scalar(keys(%{$pslxIndxRef}))-1)]);
		push(@pslx2Clust,[@temp]);
		splice(@temp);
		}  ## for(scalar(@{$clustRef})) ends
	return(@pslx2Clust);
	}  ## function ends
	###########################

	############## double calculate_ID(@a1) ##################
	sub calculate_Id{
	my($pslxIndxRef,$currTrack,$clustRef)= @_;
	my($qst,$qen,$tst,$ten,$qnum,$tnum,$match,$mismatch,$rep_match);
	my($sizeMul,$qAliSize,$tAliSize,$aliSize,$milliBad,$sizeDif,$insertFactor,$total,$isMrna);
	my($return_value);
	$match = $clustRef->[$currTrack][$pslxIndxRef->{"matches"}];
	$mismatch = $clustRef->[$currTrack][$pslxIndxRef->{"misMatch"}];
	$rep_match = $clustRef->[$currTrack][$pslxIndxRef->{"repMatch"}];
	$qnum = $clustRef->[$currTrack][$pslxIndxRef->{"qNumIn"}];
	$tnum = $clustRef->[$currTrack][$pslxIndxRef->{"tNumIn"}];
	$qst = $clustRef->[$currTrack][$pslxIndxRef->{"qSt"}];
	$qen = $clustRef->[$currTrack][$pslxIndxRef->{"qEn"}];
	$tst = $clustRef->[$currTrack][$pslxIndxRef->{"tSt"}];
	$ten = $clustRef->[$currTrack][$pslxIndxRef->{"tEn"}];
	$sizeMul = 1;
	$milliBad = 0;
	$isMrna = 1;
	$qAliSize = $sizeMul * ($qen - $qst);
	$tAliSize = $ten - $tst;
	$aliSize = shared::ArrayOperations::find_min($qAliSize, $tAliSize);    #### function to be written
		if($aliSize <= 0){
		$return_value = 0;
		}
		else{
		$sizeDif = $qAliSize - $tAliSize;    
			if($sizeDif < 0){
			    if($isMrna){
			    $sizeDif = 0;
			    }
			    else{
			    $sizeDif = -sizeDif;
			    }
			}
		$insertFactor = $qnum;
			if(!($isMrna)){
			$insertFactor += $tnum;
			}
		$total = ($sizeMul * ($match + $rep_match + $mismatch));
			if($total != 0){
			$milliBad = (1000 * ($mismatch*$sizeMul + $insertFactor + shared::CheckStringsNumericals::round((3*log(1+$sizeDif))))) / $total;
			}
		$return_value =  $milliBad;
		}  ## else ends
	 return($return_value);
	} ## function ends
	###########################

	############# double calculate_score(@a1) ###############
	sub calculate_score{
	my($pslxIndxRef,$currTrack,$clustRef)= @_;
	my($sizeMul,$match,$mismatch,$rep_match,$qnum,$tnum,$return_value);
	$match = $clustRef->[$currTrack][$pslxIndxRef->{"matches"}];
	$mismatch = $clustRef->[$currTrack][$pslxIndxRef->{"misMatch"}];
	$rep_match = $clustRef->[$currTrack][$pslxIndxRef->{"repMatch"}];
	$qnum = $clustRef->[$currTrack][$pslxIndxRef->{"qNumIn"}];
	$tnum = $clustRef->[$currTrack][$pslxIndxRef->{"tNumIn"}];
	$sizeMul = 1;
	$return_value = ($sizeMul * ($match + $rep_match) - $sizeMul * $mismatch - $qnum - $tnum);
	return($return_value);
	}  ## function ends
	###########################

	############ double calculate_coverage(@a1) #############
	sub calculate_coverage{
	my($pslxIndxRef,$currTrack,$clustRef)= @_;
	return (($clustRef->[$currTrack][$pslxIndxRef->{"qEn"}] - $clustRef->[$currTrack][$pslxIndxRef->{"qSt"}])*100.0/($clustRef->[$currTrack][$pslxIndxRef->{"qLen"}]));
	}  ## function ends
	###########################

1;

