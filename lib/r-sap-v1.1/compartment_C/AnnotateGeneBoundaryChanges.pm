package compartment_C::AnnotateGeneBoundaryChanges;

use strict;
use warnings;

require compartment_C::FindGenicRegions;

	sub get_gene_expansion_shrinkage_info{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$cutOff) = @_;

	my(@gExSt,@gExEn);

	my($annotationFlag,$annotationTag,@annotationInfo) = (0,"","");
	my($rightExpansionExternal,$rightExpansionInternal,$leftExpansionExternal,$leftExpansionInternal) = (0,0,0,0);
	my($whichIntron,$transcriptLenChangeType,$transcriptLenChange) = ("","",0);


	compartment_C::FindGenicRegions::get_exon_cords($clustRef,$currTrackIndx,$knGnIndxRef,\@gExSt,\@gExEn);

		if( (scalar(@gExSt) > 1) && ($a1Ref->[$pslx2IndxRef->{"tSt"}] < $gExSt[1]) && ($a1Ref->[$pslx2IndxRef->{"tSt"}] > $gExEn[0]) ){
		## exon has been extended in the first intron  (irrespective of strand)
			if( ($gExSt[1] - $a1Ref->[$pslx2IndxRef->{"tSt"}]) > $cutOff){
			$leftExpansionInternal = 1;
			}
		}  ## if(condition1) ends


		if( ($clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}] - $a1Ref->[$pslx2IndxRef->{"tSt"}]) > $cutOff){
		## start of the transcript has been extended in the intergenic region
		$leftExpansionExternal = 1;
		}  ## if(condition2) ends
		
		if( (scalar(@gExSt) > 1) && ($a1Ref->[$pslx2IndxRef->{"tEn"}] > $gExEn[scalar(@gExEn)-2]) && 
			($a1Ref->[$pslx2IndxRef->{"tEn"}] < $gExSt[scalar(@gExSt)-1]) ){
		## exon has been extended in to last intron
			if( ($a1Ref->[$pslx2IndxRef->{"tEn"}] - $gExEn[scalar(@gExEn)-2]) > $cutOff){
			$rightExpansionInternal = 1;
			}
		}  ## if(condition3) ends
		
		if( ($a1Ref->[$pslx2IndxRef->{"tEn"}] - $clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]) > $cutOff){
		## end of the transcript has been extended in to the intergenic region
		$rightExpansionExternal = 1;
		}  ##  if(condition 4) ends

		
		## gatehring annotation info ##
		if($leftExpansionInternal){
		$annotationFlag = 1;
		$whichIntron = compartment_C::FindGenicRegions::which_intron(1,scalar(@gExSt),$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]);
		$transcriptLenChangeType = "contraction";
		$transcriptLenChange = $gExEn[0]-$gExSt[0];  ## possible exclusion of the first exon (irrespective of strand)
		}
		elsif($leftExpansionExternal){
		$annotationFlag = 1;
		$whichIntron = "Intergenic";
		$transcriptLenChangeType = "extension";
		$transcriptLenChange = ($clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}] - $a1Ref->[$pslx2IndxRef->{"tSt"}]);
		}
		elsif($rightExpansionInternal){
		$annotationFlag = 1;
		$whichIntron  = compartment_C::FindGenicRegions::which_intron( (scalar(@gExSt)-1),scalar(@gExSt),$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]);
		$transcriptLenChangeType = "contraction";
		$transcriptLenChange = $gExEn[scalar(@gExEn)-1] - $gExSt[scalar(@gExSt)-1];  ## exclusion of last intron (irrespective of strand);
		}
		elsif($rightExpansionExternal){
		$annotationFlag = 1;
		$whichIntron = "Intergenic";
		$transcriptLenChangeType = "extension";
		$transcriptLenChange = ($a1Ref->[$pslx2IndxRef->{"tEn"}] - $clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]);
		}
		
		if($annotationFlag){
		$annotationTag = polyA_or_tss($rightExpansionExternal,$rightExpansionInternal,$leftExpansionExternal,
					$leftExpansionInternal,$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]);
		}		

	push(@annotationInfo,($whichIntron,$transcriptLenChangeType,$transcriptLenChange));
	return($annotationFlag,"GeneBoundaryChange",$annotationTag,@annotationInfo);
	}  # function ends
	###########################
	sub polyA_or_tss{
	my($rightExpansion1,$rightExpansion2,$leftExpansion1,$leftExpansion2,$strand) = @_;
	my $annotationTag = "";
		if( ($leftExpansion1) || ($leftExpansion2) ){
			if($strand eq '+'){
			$annotationTag = "AlternativeTSS";
			}
			else{
			$annotationTag = "AlternativePolyA";
			}
		}  ## if ends
		elsif( ($rightExpansion1) || ($rightExpansion2) ){
			if($strand eq '+'){
			$annotationTag = "AlternativePolyA";
			}
			else{
			$annotationTag = "AlternativeTSS";
			}
		}  ## else ends
		
	return($annotationTag);
	}  ## function ends
	###########################
1;
		