package compartment_C::AnnotateIntronsOnly;

use strict;
use warnings;

require compartment_C::Pslx2RecordProcessing;
require compartment_C::FindGenicRegions;


	sub get_introns_only_annot_info{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx) = @_;
	my(@exSt,@exEn,@gExSt,@gExEn,@intronSt,@intronEn);
	my($i,$j,@annotInfo);
	my($extension,$st,$en);
	my(@introns);

	compartment_C::Pslx2RecordProcessing::get_pslx2_algn_blocks($pslx2IndxRef,$a1Ref,\@exSt,\@exEn);
	compartment_C::FindGenicRegions::get_exon_cords($clustRef,$currTrackIndx,$knGnIndxRef,\@gExSt,\@gExEn);
	compartment_C::FindGenicRegions::get_intron_regions(\@gExSt,\@gExEn,\@intronSt,\@intronEn);

		for($i = 0;$i<scalar(@exSt);++$i){
			for($j = 0;$j<scalar(@intronSt);++$j){
			($extension,$st,$en) = shared::CoordinateComparisons::get_intersecting_coordinates($exSt[$i],$exEn[$i],
							$intronSt[$j],$intronEn[$j]);
				if($extension){
				push(@introns,compartment_C::FindGenicRegions::which_intron($j,scalar(@gExSt),$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]));
				}
			}  ## for($j) ends
		}  ## for($i) ends
	@introns = shared::ArrayOperations::remove_array_redundancy(@introns);
	push(@annotInfo,scalar(@introns),join(',',@introns));
	return(@annotInfo);
	}  ## function ends
	###########################
1;
