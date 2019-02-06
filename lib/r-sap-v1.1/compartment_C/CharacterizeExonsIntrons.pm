#### Characterize introns ####

package compartment_C::CharacterizeExonsIntrons;

use strict;
use warnings;

require compartment_C::IntersectWithGenicRegions;

	sub characterize_exons{
	my($gExStRef,$gExEnRef,$knGnIndxRef,$clustRef,$currTrackIndx) = @_;
	my(@exonChars,@temp1,@temp2,@regions);
	my($i,$j);
	
		for($i = 0;$i<scalar(@{$gExStRef});++$i){
		push(@temp1,$gExStRef->[$i]);
		push(@temp2,$gExEnRef->[$i]);
		@regions = compartment_C::IntersectWithGenicRegions::find_intersected_genic_regions(\@temp1,\@temp2,$knGnIndxRef,$clustRef,$currTrackIndx);
		push(@exonChars,join(',',@regions));
		shared::DeleteDataStructures::undefine_arrays(\@temp1,\@temp2,\@regions);
		}  ## for($i:all exons) ends
	return(@exonChars);
	}  ## function ends
	###########################
	sub get_closest_left_exon{
	my($exOverlapsRef,$currInIndx) = @_;
	my($i,$closest);
	$closest = "NA";  ## no left closest exon
		for($i = $currInIndx;$i >= 0;--$i){
			if($exOverlapsRef->[$i]){
			$closest = $i;
			last;
			}
		}
	return($closest);
	}  ## function ends
	###########################	
	sub get_closest_right_exon{
	my($exOverlapsRef,$currInIndx) = @_;
	my($i,$closest);
	$closest = "NA";  ## no left closest exon
		for($i = ($currInIndx+1);$i<scalar(@{$exOverlapsRef});++$i){
			if($exOverlapsRef->[$i]){
			$closest = $i;
			last;
			}
		}
	return($closest);
	}  ## function ends
	###########################	
	sub which_exon_to_annotate{
	my($leftIndx,$rightIndx,$exOverlapsRef) = @_;
	my($leftOverlap,$rightOverlap);
	my $whichIndx = "NA";
		if($leftIndx eq "NA"){
		$whichIndx = 2;  ## right exon
		}
		elsif($rightIndx eq "NA"){
		$whichIndx = 1; ## left exon
		}
		elsif( ($leftIndx ne "NA") && ($rightIndx ne "NA") ){
			if($exOverlapsRef->[$leftIndx] > $exOverlapsRef->[$rightIndx]){
			$whichIndx = 1;
			}
			elsif($exOverlapsRef->[$rightIndx] > $exOverlapsRef->[$leftIndx]){
			$whichIndx = 2;
			}
			elsif($exOverlapsRef->[$rightIndx] == $exOverlapsRef->[$leftIndx]){
			$whichIndx = 3;
			}
		}  ## elsif() ends
	return($whichIndx);
	}  ## function ends
	###########################
	sub get_left_region{
	my($exCharsRef,$leftIndx) = @_;
		if($leftIndx ne  "NA"){
		my @temp = split(',',$exCharsRef->[$leftIndx]);
		return($temp[scalar(@temp)-1]);
		}
		else{
		return '';  ## NULL string
		}
	}  ## function end
	###########################
	sub get_right_region{
	my($exCharsRef,$rightIndx) = @_;
		if($rightIndx ne  "NA"){
		my @temp = split(',',$exCharsRef->[$rightIndx]);
		return($temp[0]);
		}
		else{
		return '';
		}
	}  ## function ends
	###########################
	