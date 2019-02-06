package compartment_C::AnnotateNeighboringGenes;

use strict;
use warnings;


	sub get_neighboring_genes_annot_info{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx) = @_;
	
	my($distance,$expandedEnd);
	my(@annotInfo);
	## neighboring genes: Which gene got expanded, which side, how much is the distance
		##between the closest ends of the read and the expanded gene	

	## function to be written
	$distance = shared::CoordinateComparisons::calculate_end_to_end_distance($a1Ref->[$pslx2IndxRef->{"tSt"}],$a1Ref->[$pslx2IndxRef->{"tEn"}],
			$clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],
			$clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]);
		
		if($a1Ref->[$pslx2IndxRef->{"tEn"}] <= $clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}]){
			if($clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}] eq '+'){
			$expandedEnd = "5End";
			}
			else{
			$expandedEnd = "3End";
			}
		}  ## if(read is on the left side of the gene) ends
		elsif($a1Ref->[$pslx2IndxRef->{"tSt"}] >= $clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]){
			if($clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}] eq '+'){
			$expandedEnd = "3End";
			}
			else{
			$expandedEnd = "5End";
			}
		}  ## else ends
		else{
		$expandedEnd = "overlap";
		}
	push(@annotInfo,($expandedEnd,$distance));
	return(@annotInfo);
	}  ## function ends
	###########################
1;
	
		