package compartment_C::PrintUncharacterizedReads;

use strict;
use warnings;

	sub check_if_artifact{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$cutOff,$annotStringency,$annotTagsRef,$outFilesRef) = @_;

	my $uncharacterizedFlag = 0;  ## false by default	
	my @tempClust;
	push(@tempClust,[@{$a1Ref}]);


		if( join('',@{$annotTagsRef}) =~ m/(uncharacterized)/i){  ## if uncharacterized annotation found
			if($annotStringency eq "unique"){  ## if annotation stringency is unique (do nothing for multi)
			shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{"Uncharacterized"},scalar(keys(%{$pslx2IndxRef})),1,\@tempClust);
			}
		$uncharacterizedFlag = 1;
		}  ## if (annotationTags =~ uncharacterized);
		else{
			if( ($clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}] - $a1Ref->[$pslx2IndxRef->{"tSt"}]) > $cutOff){
			shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{"Uncharacterized"},scalar(keys(%{$pslx2IndxRef})),1,\@tempClust);
			$uncharacterizedFlag = 1;
			}
			elsif( ($a1Ref->[$pslx2IndxRef->{"tEn"}] - $clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]) > $cutOff ){
			shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{"Uncharacterized"},scalar(keys(%{$pslx2IndxRef})),1,\@tempClust);
			$uncharacterizedFlag = 1;	
			}
			elsif( ($clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}] > $a1Ref->[$pslx2IndxRef->{"tSt"}]) && 
				($clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}] < $a1Ref->[$pslx2IndxRef->{"tEn"}]) ){
			shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{"Uncharacterized"},scalar(keys(%{$pslx2IndxRef})),1,\@tempClust);
			$uncharacterizedFlag = 1;	
			}

		}  ## else ends
	return($uncharacterizedFlag);  ## returns true if reads annotation type is "Uncharacterized"
	}  ## function ends
	###########################
1;
