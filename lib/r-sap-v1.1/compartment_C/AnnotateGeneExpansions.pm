package compartment_C::AnnotateGeneExpansions;

use strict;
use warnings;



	sub get_gene_expansion_annot_info{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx) = @_;
	my(@expandedBases,@expandedRegions);
	my(@annotInfo);

		if($a1Ref->[$pslx2IndxRef->{"tSt"}] < $clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}] ){  ## start position expansion
		push(@expandedBases,($clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}] - $a1Ref->[$pslx2IndxRef->{"tSt"}]));
			if($clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}] eq '+'){
			push(@expandedRegions,"5End");
			}
			else{
			push(@expandedRegions,"3End");
			}
		}  ## start position expansion ends
		if($a1Ref->[$pslx2IndxRef->{"tEn"}] > $clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]){  ## end position expanded
		push(@expandedBases,($a1Ref->[$pslx2IndxRef->{"tEn"}] - $clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]));
			if($clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}] eq '+'){
			push(@expandedRegions,"3End");
			}
			else{
			push(@expandedRegions,"5End");
			}
		}  ## end position expansion ends
	push(@annotInfo,(join(',',@expandedBases),join(',',@expandedRegions)));
	return(@annotInfo);
	}  ## function ends
	###########################
1;
	
	