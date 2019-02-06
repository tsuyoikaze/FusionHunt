############ Module file for chimer detection at left. ##############
package compartment_B::DetectRightChimerasSR;

use strict;
use warnings;

	sub define_potential_right_chimers{  ## function to detect all the possible alignment fillers. These alignments will be screened again later
	my($gapCutoff,$repeatIdentityCutoff,$scoreCutoff,$qLenCutoff,$pslx2IndxRef,$clustRef) = @_;

	my ($gpSt,$gpEn,$gapLen,@covClust);
	my ($c,$OvrLap,$i,@temp);

	$c = 0;
	$gpSt = $clustRef->[0][$pslx2IndxRef->{"qEn"}];  ## top hit end (start of the gap)
	$gpEn = $clustRef->[0][$pslx2IndxRef->{"qLen"}]; ## query length (gap end)
	$gapLen = ($gpEn-$gpSt+1);  ## gap length
		
		if($gapLen >= $gapCutoff){  ## if the gap, to be covered, is long enough
		
			for($i = 1;$i<scalar(@{$clustRef});++$i){  ## for all hits other than the top hit
			@temp = @{$clustRef->[$i]}[0..(scalar(keys(%{$pslx2IndxRef}))-1)];
				
			 	if( ($temp[$pslx2IndxRef->{"identity"}] >= $repeatIdentityCutoff) && ($temp[$pslx2IndxRef->{"score"}] >= $scoreCutoff) && ($temp[$pslx2IndxRef->{"qLen"}] >= $qLenCutoff) ){
					if( ($temp[$pslx2IndxRef->{"qSt"}] >= $gpSt) && (($temp[$pslx2IndxRef->{"qSt"}] - $gpSt) < $gapCutoff) ){  ## first case: when no overlap
			  		$c = $temp[$pslx2IndxRef->{"qEn"}] - ($temp[$pslx2IndxRef->{"qSt"}]+1) + 1;
			  		$OvrLap = 0;
						if($c >= $gapCutoff){
						push(@temp,$c);
			  			push(@covClust,[@temp]);
			  			}
			  		}  ## if  ends
	
			  		if( ($temp[$pslx2IndxRef->{"qEn"}] > $gpSt) && ($temp[$pslx2IndxRef->{"qSt"}] < $gpSt) ){  ## second case: overlap
			  		$c = $temp[$pslx2IndxRef->{"qEn"}] - $gpSt;
			  		$OvrLap = $gpSt-$temp[$pslx2IndxRef->{"qSt"}];
			  		 	if($c >= $gapCutoff && ($OvrLap < ($temp[$pslx2IndxRef->{"qLen"}]/3)) && ($c > $OvrLap) ){
			  			push(@temp,$c);
			  			push(@covClust,[@temp]); 
			         		}  ## if(coverage condition) ends
			         	}  ## if ends	
				}  ## if((id>=$repeatIdentiyCutoff) && ($score>=$scoreCutoff)) ends
			splice(@temp);
			}  ###for(all other hits) ends
		}  ## if($lt >= $gapCutoff) ends
	return(@covClust);
	}  ## function ends
	###########################
	sub print_right_chimeras_to_file{
	my($pslx2IndxRef,$clustRef,$covClustRef,$qualifiedAlignmentsRef)= @_;
	my($i,@outArray);
	my $lt = $clustRef->[0][$pslx2IndxRef->{"qLen"}] - $clustRef->[0][$pslx2IndxRef->{"qEn"}];   ## length of the gap
		
		if(scalar(@{$qualifiedAlignmentsRef})){
			for($i = 0;$i<scalar(@{$qualifiedAlignmentsRef});++$i){
			push(@outArray,($clustRef->[0][$pslx2IndxRef->{"qID"}],$clustRef->[0][$pslx2IndxRef->{"qLen"}]));
			push(@outArray,"(".$clustRef->[0][$pslx2IndxRef->{"tID"}]."-".$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"tID"}].")");
			push(@outArray,"(".$clustRef->[0][$pslx2IndxRef->{"identity"}]."-".$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"identity"}].")");
			push(@outArray,"(".$lt.",".($lt-$covClustRef->[$qualifiedAlignmentsRef->[$i]][(scalar(keys%{$pslx2IndxRef})-1)+1]).")");
			push(@outArray,$covClustRef->[$qualifiedAlignmentsRef->[$i]][(scalar(keys%{$pslx2IndxRef})-1)+1]);
			push(@outArray,($clustRef->[0][$pslx2IndxRef->{"qEn"}]-$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"qSt"}]));  ############# Overlap #########
			push(@outArray,"(".$qualifiedAlignmentsRef->[$i].")");
			push(@outArray,"(".$clustRef->[0][$pslx2IndxRef->{"qSt"}]."-".$clustRef->[0][$pslx2IndxRef->{"qEn"}]);
			push(@outArray,$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"qSt"}]."-".$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"qEn"}].")");
			push(@outArray,("(".$clustRef->[0][$pslx2IndxRef->{"score"}],$clustRef->[0][$pslx2IndxRef->{"coverage"}],$clustRef->[0][$pslx2IndxRef->{"strand"}],
			$clustRef->[0][$pslx2IndxRef->{"tLen"}],$clustRef->[0][$pslx2IndxRef->{"tSt"}],$clustRef->[0][$pslx2IndxRef->{"tEn"}],
			$clustRef->[0][$pslx2IndxRef->{"blNum"}],$clustRef->[0][$pslx2IndxRef->{"blLen"}],$clustRef->[0][$pslx2IndxRef->{"qBlSt"}],
			$clustRef->[0][$pslx2IndxRef->{"tBlSt"}].")"));
			push(@outArray,("(".$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"score"}],$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"coverage"}],
			$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"strand"}],$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"tLen"}],
			$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"tSt"}],$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"tEn"}],
			$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"blNum"}],$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"blLen"}],
			$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"qBlSt"}],$covClustRef->[$qualifiedAlignmentsRef->[$i]][$pslx2IndxRef->{"tBlSt"}].")"));
			##splice(@outArray);
			}  ## for($i) ends
		}  ## if(qualifiedAlignments) ends
		else{
		}
	return(@outArray);
	}  ## function ends
	###########################
	
##################################
1;

