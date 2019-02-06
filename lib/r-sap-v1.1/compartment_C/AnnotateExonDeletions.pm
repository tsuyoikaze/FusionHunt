package compartment_C::AnnotateExonDeletions;

use strict;
use warnings;

require compartment_C::Pslx2RecordProcessing;
require compartment_C::FindGenicRegions;
require compartment_C::IntersectWithGenicRegions;


	sub get_exon_deletion_annot_info{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$deletionCutoff) = @_;
	my(@exSt,@exEn,@gExSt,@gExEn);
	my($i,$j,@annotInfo);
	my($st,$en,@intSt,@intEn);  ## intersect coordinates
	my($alignmentDeletion,$fullExonSkipping);  ## temporary variables
	my($totalDeletedBases,@exonDeletions,@exonsDeleted);
	my($skippedExonsNum,@skippedExons);
	my(@regions,$geneLen);
	$totalDeletedBases = 0;
	$skippedExonsNum = 0;

	compartment_C::Pslx2RecordProcessing::get_pslx2_algn_blocks($pslx2IndxRef,$a1Ref,\@exSt,\@exEn);
	compartment_C::FindGenicRegions::get_exon_cords($clustRef,$currTrackIndx,$knGnIndxRef,\@gExSt,\@gExEn);

		for($i = 0;$i<(scalar(@exSt)-1);++$i){
			for($j = 0;$j<scalar(@gExSt);++$j){

			($alignmentDeletion,$st,$en) = 	shared::CoordinateComparisons::get_intersecting_coordinates($exEn[$i],$exSt[$i+1],$gExSt[$j],$gExEn[$j]);

			$fullExonSkipping = $alignmentDeletion/($gExEn[$j]-$gExSt[$j]);

				if($alignmentDeletion > $deletionCutoff){
				$totalDeletedBases += $alignmentDeletion;
				}

				if($fullExonSkipping == 1){
				$skippedExonsNum += 1;
				push(@skippedExons,compartment_C::FindGenicRegions::which_exon($j,scalar(@gExSt),$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]));
				}

				if( ($alignmentDeletion > $deletionCutoff) || ($fullExonSkipping == 1) ){
				push(@intSt,$st);  ## pushing intersect coordinates
				push(@intEn,$en);
				push(@exonDeletions,($en-$st));  ## deleted bases from each deleted exon
				push(@exonsDeleted,compartment_C::FindGenicRegions::which_exon($j,scalar(@gExSt),$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]));  ## which exons have deletions
				}
			}  ## for($j) ends
		}  ## for($i) ends

		if(scalar(@intSt)){  ## if there is deletion
		@regions = compartment_C::IntersectWithGenicRegions::find_intersected_genic_regions(\@intSt,\@intEn,$knGnIndxRef,$clustRef,$currTrackIndx);
		@regions = shared::ArrayOperations::remove_array_redundancy(@regions);
		$geneLen = $clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]-$clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}];

		
		## deleted exons and number of bases deleted
		push(@annotInfo,join(',',@exonsDeleted),(join',',@exonDeletions));
		push(@annotInfo,(join(',',@regions),$skippedExonsNum));
			if($skippedExonsNum){
			push(@annotInfo,join(',',@skippedExons));
			}
			else{
			push(@annotInfo,'NA');
			}
		}
		return(@annotInfo);
	}  ## function ends
	###########################
1;
