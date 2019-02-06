package compartment_C::AnnotateExonExtensions;

use strict;
use warnings;

require compartment_C::Pslx2RecordProcessing;
require compartment_C::FindGenicRegions;
require compartment_C::CharacterizeExonsIntrons;


	sub get_exon_extension_annot_info{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$cutOff) = @_;
	my(@exSt,@exEn,@gExSt,@gExEn,@intronSt,@intronEn);
	my($exonNums,$exons);
	my($extension,$st,$en);
	my($totalExtension,@exonExtensions,@introns);
	my($compIntronInclusionNums,@intronsIncluded);
	my($geneLen,$i,$j,,@annotInfo);
	my($annotationFlag,$annotationTag) = (0,"");

	$annotationFlag = 0;
	$annotationTag = "";
	$totalExtension = 0;
	$compIntronInclusionNums = 0;
	splice(@intronsIncluded);  ## completely introns retained
	splice(@exonExtensions);  ## number of bases extended
	splice(@introns);  ## introns extended into


	#### variables to find out which genic region got extended ####
	my($overlap,@exonOverlaps,@exonChars,@regions,@temp);
	my($closestLeftExon,$closestRightExon,$whichExonAnnotate);

	
	compartment_C::Pslx2RecordProcessing::get_pslx2_algn_blocks($pslx2IndxRef,$a1Ref,\@exSt,\@exEn);  ## alignment blocks
	compartment_C::FindGenicRegions::get_exon_cords($clustRef,$currTrackIndx,$knGnIndxRef,\@gExSt,\@gExEn);  ## exon coordinates
	compartment_C::FindGenicRegions::get_intron_regions(\@gExSt,\@gExEn,\@intronSt,\@intronEn);  ## intron coordinates

		## getting overlap of the alignment blocks with the known exons ##
		for($j = 0;$j<scalar(@gExSt);++$j){  ## for all the exons
			for($i = 0;$i<scalar(@exSt);++$i){  ## for all the alignment blocks
			($overlap,$st,$en) = shared::CoordinateComparisons::get_intersecting_coordinates($exSt[$i],$exEn[$i],$gExSt[$j],$gExEn[$j]);
			$exonOverlaps[$j] += $overlap;
			}  ## for($j) ends
		}  ## for($i) ends
		####

	@exonChars  = compartment_C::CharacterizeExonsIntrons::characterize_exons(\@gExSt,\@gExEn,$knGnIndxRef,$clustRef,$currTrackIndx);

		if(scalar(@exonOverlaps) != scalar(@exonChars) ){
		die "\n overlap and char arrays are not equal \n";
		}
		for($i = 0;$i<scalar(@exSt);++$i){  ## for all the alignment blocks
			for($j = 0;$j<scalar(@intronSt);++$j){  ## for all the introns
			($extension,$st,$en) = shared::CoordinateComparisons::get_intersecting_coordinates($exSt[$i],$exEn[$i],$intronSt[$j],$intronEn[$j]);
				if(!($intronEn[$j]-$intronSt[$j])){
				next;
				}
				if($extension > $cutOff){  ## extension to jth intron
				$totalExtension += $extension;
				push(@exonExtensions,$extension);
				push(@introns,compartment_C::FindGenicRegions::which_intron($j,scalar(@gExSt),$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]));

				#### determining the genic region affected ####
				$closestLeftExon = compartment_C::CharacterizeExonsIntrons::get_closest_left_exon(\@exonOverlaps,$j);  ## function to be written
				$closestRightExon = compartment_C::CharacterizeExonsIntrons::get_closest_right_exon(\@exonOverlaps,$j);  ## function to be written
								
				$whichExonAnnotate = compartment_C::CharacterizeExonsIntrons::which_exon_to_annotate($closestLeftExon,$closestRightExon,\@exonOverlaps);  ## function to be written

					if($whichExonAnnotate == 1){  ## left exon
					push(@regions,compartment_C::CharacterizeExonsIntrons::get_left_region(\@exonChars,$closestLeftExon));
					}
					elsif($whichExonAnnotate == 2){  ## right exon
					push(@regions,compartment_C::CharacterizeExonsIntrons::get_right_region(\@exonChars,$closestRightExon));
					}
					else{  ## if both
					push(@regions,compartment_C::CharacterizeExonsIntrons::get_left_region(\@exonChars,$closestLeftExon));
					push(@regions,compartment_C::CharacterizeExonsIntrons::get_right_region(\@exonChars,$closestRightExon));
					}
				}  ## if($extension > $cutOff) ends

				if( ($extension/($intronEn[$j]-$intronSt[$j])) == 1 ){
				$compIntronInclusionNums += 1;
				push(@intronsIncluded,compartment_C::FindGenicRegions::which_intron($j,scalar(@gExSt),$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]));

				$closestLeftExon = compartment_C::CharacterizeExonsIntrons::get_closest_left_exon(\@exonOverlaps,$j);  ## function to be written
				$closestRightExon = compartment_C::CharacterizeExonsIntrons::get_closest_right_exon(\@exonOverlaps,$j);  ## function to be written
				push(@regions,compartment_C::CharacterizeExonsIntrons::get_left_region(\@exonChars,$closestLeftExon));
				push(@regions,compartment_C::CharacterizeExonsIntrons::get_right_region(\@exonChars,$closestRightExon));
				}  ## if(complete intron retention) ends
			}  ## for($j) ends
		}  ## for($i) ends

		if($totalExtension){
		$annotationFlag = 1;
		$annotationTag = "InternalExonExtension";
		}

	@regions = shared::ArrayOperations::remove_array_redundancy(@regions);
	push(@annotInfo,(join(',',@introns),join(',',@exonExtensions)));
	push(@annotInfo,join(',',@regions));
	push(@annotInfo,$compIntronInclusionNums);
		if(scalar(@intronsIncluded)){
		push(@annotInfo,join(',',@intronsIncluded));
		}
		else{
		push(@annotInfo,'NA');
		}
	shared::DeleteDataStructures::undefine_arrays(\@intronSt,\@intronEn,\@regions);
	return($annotationFlag,$annotationTag,@annotInfo);
	}  ## function ends
	###########################
1;
