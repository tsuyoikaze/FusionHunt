##########################################################################################################################
################################################ Om Ganeshaya Namah ######################################################
##########################################################################################################################

################ This perl script was written to start with the initial processing of pslx files. ########################
################ Score, Covrage and Identity have been calculated and reads have been sorted on ##########################
################ score, identity and query coverage before running this script. ##########################################

package compartment_A::AlignmentBasedClassification;
use strict ;
use warnings;

use compartment_A::ClassifyBasedOnAlignmentScoresSR;

#### main subroutine ####
	sub alignment_based_classification{
	my($cutoffsRef,$pslx2IndxRef,$clustRef,$outFilesRef,$readClassHashRef) = @_;
	my $l = scalar(keys(%{$pslx2IndxRef}));
	my($similarHits,$alignmentAnnotTag);
	
	##shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{"TopHits"},$l,1,$clustRef);  ## printing top hit to a file 

	$similarHits = check_non_uniqueness($pslx2IndxRef,$clustRef);  ## checking for multi-hits

			if($similarHits > 1){                       ## if cluster has non unique mapping to genome
			$alignmentAnnotTag = "MultiHits";
			shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{$alignmentAnnotTag},$l,$similarHits,$clustRef);  
			++$readClassHashRef->{"MultiHits"};  ## increasing stats count
			}  ## if($similarHits>1) ends	

			else{  ## if no multi hits (unique mapping)
			shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{"TopHits"},$l,1,$clustRef);  ## printing top hit to a file 
			$alignmentAnnotTag = check_ID_cov_cutoff($pslx2IndxRef,$cutoffsRef,$clustRef);

				if($alignmentAnnotTag ne "PotentialFusion"){  ## if not potential fusion (annotate or discarded)
					if($alignmentAnnotTag eq "Discarded"){  ## if discarded
					shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{$alignmentAnnotTag},$l,1,$clustRef);  ## printing psl2 record accroding to the alignment tags
					++$readClassHashRef->{"Discarded"};
					}
				}  ## if($CutoffCheck ne "PotentialFusion") ends

				elsif($alignmentAnnotTag eq  "PotentialFusion"){  ## if(identity is OK but coverage is low)
				## change function
				$alignmentAnnotTag = check_if_fusion($pslx2IndxRef,$cutoffsRef,$clustRef);
					if($alignmentAnnotTag eq "Discarded"){  ## if no fusion detected (discard)
					shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{$alignmentAnnotTag},$l,1,$clustRef);  ## printing psl
					++$readClassHashRef->{"Discarded"};
					}
					##else{
					##shared::FileCheckPrint::::print_cluster_to_file($outFilesRef->{$alignmentAnnotTag},$l,scalar(@{$clustRef}),$clustRef);  ## printing psl
					##}
				}  ## elsif($CutoffCheck eq "PotentialFusion") ends
				#### alignment integration  ends 
			}  ## else ends
	return($alignmentAnnotTag);
	}  ## function ends
	###########################
1;
