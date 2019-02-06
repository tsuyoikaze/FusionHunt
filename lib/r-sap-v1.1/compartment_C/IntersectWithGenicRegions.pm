package compartment_C::IntersectWithGenicRegions;

use strict;
use warnings;

require compartment_C::FindGenicRegions;

	sub find_intersected_genic_regions{
	my($stArrayRef,$enArrayRef,$knGnIndxRef,$clustRef,$currTrackIndx) = @_;
	my(@exSt,@exEn,$cdsSt,$cdsEn,$strand);  ## rerefence info to be passed on to different functions
	my(@regSt,@regEn);  ## coodrinates for region start, end
	my(@regions);

	compartment_C::FindGenicRegions::get_cds_pos($clustRef,$currTrackIndx,$knGnIndxRef,\$cdsSt,\$cdsEn);  ## coding star and coding end
	compartment_C::FindGenicRegions::get_exon_cords($clustRef,$currTrackIndx,$knGnIndxRef,\@exSt,\@exEn); ## exon start and exon end

	$strand = $clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}];  ## gene strand
	
	
		if($cdsEn-$cdsSt){  ## if gene has protein coding region
			if($strand eq '+'){
			compartment_C::FindGenicRegions::get_5UTR_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,$strand,\@regSt,\@regEn);
				if(shared::CoordinateComparisons::find_block_overlaps($stArrayRef,$enArrayRef,\@regSt,\@regEn)){  ## if overlap with the 5'UTR
				push(@regions,"5UTR");
				}
			shared::DeleteDataStructures::undefine_arrays(\@regSt,\@regEn);
			
			compartment_C::FindGenicRegions::get_coding_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,\@regSt,\@regEn);
				if(shared::CoordinateComparisons::find_block_overlaps($stArrayRef,$enArrayRef,\@regSt,\@regEn)){  ## if overlap with the coding regions
				push(@regions,"CDS");
				}
			shared::DeleteDataStructures::undefine_arrays(\@regSt,\@regEn);
			
			compartment_C::FindGenicRegions::get_3UTR_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,$strand,\@regSt,\@regEn);
				if(shared::CoordinateComparisons::find_block_overlaps($stArrayRef,$enArrayRef,\@regSt,\@regEn)){ ## if overlap with the 3'UTR
				push(@regions,"3UTR");
				}
			shared::DeleteDataStructures::undefine_arrays(\@regSt,\@regEn);
			}  ## if(strand eq '+') ends
			else{
			compartment_C::FindGenicRegions::get_3UTR_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,$strand,\@regSt,\@regEn);
				if(shared::CoordinateComparisons::find_block_overlaps($stArrayRef,$enArrayRef,\@regSt,\@regEn)){ ## if overlap with the 3'UTR
				push(@regions,"3UTR");
				}
			shared::DeleteDataStructures::undefine_arrays(\@regSt,\@regEn);

			compartment_C::FindGenicRegions::get_coding_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,\@regSt,\@regEn);
				if(shared::CoordinateComparisons::find_block_overlaps($stArrayRef,$enArrayRef,\@regSt,\@regEn)){  ## if overlap with the coding regions
				push(@regions,"CDS");
				}
			shared::DeleteDataStructures::undefine_arrays(\@regSt,\@regEn);
		
			compartment_C::FindGenicRegions::get_5UTR_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,$strand,\@regSt,\@regEn);
				if(shared::CoordinateComparisons::find_block_overlaps($stArrayRef,$enArrayRef,\@regSt,\@regEn)){  ## if overlap with the 5'UTR
				push(@regions,"5UTR");
				}
			shared::DeleteDataStructures::undefine_arrays(\@regSt,\@regEn);
			}  ## elsif(strand eq '-') ends
		}  ## if(coding) ends			
		else{  ## if gene is non-coding
		compartment_C::FindGenicRegions::get_noncoding_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,\@regSt,\@regEn);  ## if overlap with the non-coding region
			if(shared::CoordinateComparisons::find_block_overlaps($stArrayRef,$enArrayRef,\@regSt,\@regEn)){
			push(@regions,"NonCoding");
			}
		}  ## else(!(coding region)) ends
	return(@regions);
	}  ## function ends
	###########################
1;


