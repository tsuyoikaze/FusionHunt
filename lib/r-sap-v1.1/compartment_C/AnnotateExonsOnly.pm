package compartment_C::AnnotateExonsOnly;

use strict;
use warnings;

require compartment_C::Pslx2RecordProcessing;
require compartment_C::FindGenicRegions;
require compartment_C::IntersectWithGenicRegions;


	sub get_exons_only_annot_info{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx) = @_;
	my($intersect,$st,$en,@intSt,@intEn);  ## intersect coordinates
	my($i,$j,,@annotInfo);
	my(@exSt,@exEn,@gExSt,@gExEn);
	my($exonNums,$exons,@regions);

	## function needs to be defined ##
	compartment_C::Pslx2RecordProcessing::get_pslx2_algn_blocks($pslx2IndxRef,$a1Ref,\@exSt,\@exEn);
	compartment_C::FindGenicRegions::get_exon_cords($clustRef,$currTrackIndx,$knGnIndxRef,\@gExSt,\@gExEn);
	
	$exonNums = 0;  ## number of exons
	$exons = "";  ## which exons

		for($i = 0;$i<scalar(@exSt);++$i){
			for($j = 0;$j<scalar(@gExSt);++$j){
			## function to be defined
			($intersect,$st,$en) = shared::CoordinateComparisons::get_intersecting_coordinates($exSt[$i],$exEn[$i],$gExSt[$j],$gExEn[$j]);
				if($intersect){  ## if there is any intersect
				push(@intSt,$st);
				push(@intEn,$en);
				++$exonNums;
				#$exons .= ($j+1).',';
				$exons .= compartment_C::FindGenicRegions::which_exon($j,scalar(@gExSt),$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]).',';
				}
			}  ## for($j) ends
		}  ## for($i) ends

	## function to be defined
	@regions = compartment_C::IntersectWithGenicRegions::find_intersected_genic_regions(\@intSt,\@intEn,$knGnIndxRef,$clustRef,$currTrackIndx);
	## function to be defined
	@regions = shared::ArrayOperations::remove_array_redundancy(@regions);
	push(@annotInfo,($exonNums,$exons,join(',',@regions)));
	return(@annotInfo);
	}  ## function ends
	###########################
	sub update_exons_only_hash{
	my($regionsRef,$exonsOnlyHashRef) = @_;
	my($tempStr,@tempArray);

		if(scalar(@{$regionsRef}) > 0){
		$tempStr = join(',',@{$regionsRef});
		@tempArray = split(',',$tempStr);	
		@tempArray = shared::ArrayOperations::remove_array_redundancy(@tempArray);
			if(scalar(@tempArray) > 1){
			$exonsOnlyHashRef->{"MultiMapping"} += 1;
			}
			elsif(scalar(@tempArray) == 1){
			$exonsOnlyHashRef->{$tempArray[0]} += 1;
			}
		}  ## if(scalar(@regionsRef) > 0) ends
	}  ## function ends
	###########################
		
			



1;
