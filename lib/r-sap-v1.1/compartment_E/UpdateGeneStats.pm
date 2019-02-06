package compartment_E::UpdateGeneStats;

use strict;
use warnings;

	sub update_gene_dist_stats{
	my($key,$a1Ref,$geneFileIndxRef,$geneCharHashRef,$geneDistHashRef) = @_;
	my ($i,@temp);
	
		if($a1Ref->[$geneFileIndxRef->{"annotTags"}] =~ m/[a-z]/i){  ## there is any annotation tag defined (gene is expressed in the sample)
		$geneDistHashRef->{"detectedGenes"} += 1;

			if($a1Ref->[$geneFileIndxRef->{"codingFlag"}] eq "coding"){
			$geneDistHashRef->{"codingDetected"} += 1;
			}
			else{
			$geneDistHashRef->{"non-codingDetected"} += 1;
			}

		$a1Ref->[$geneFileIndxRef->{"annotTags"}] =~ s/\s//g;
		@temp = split(',',$a1Ref->[$geneFileIndxRef->{"annotTags"}]);

			for($i = 0;$i<scalar(@temp);++$i){
			$geneCharHashRef->{$temp[$i]} += 1;  ## updating the stats
			} ## for(scalar(@temp)) ends

			if($a1Ref->[$geneFileIndxRef->{"multiAnnotFlag"}]){  ## if gene is mapped by multi-annotation reads
			$geneCharHashRef->{"MultipleAnnotations"} += 1;
			}


		shared::ArrayOperations::remove_array_element(\@temp,"ExonOnly");
		shared::ArrayOperations::remove_array_element(\@temp,"NeighboringGenes");

			if(!(scalar(@temp))){  ## if no element left
			$geneDistHashRef->{"normallySpliced"} += 1;
			}
			elsif(scalar(@temp) == 1){
			$geneDistHashRef->{"oneAbSplicing"} += 1;
			}

			elsif( (scalar(@temp) > 1) || ($a1Ref->[$geneFileIndxRef->{"multiAnnotFlag"}]) ){
			$geneDistHashRef->{"multiAbSplicing"} += 1;
			}

		}  ## if(annotationTags exist) ends
	splice(@temp);
	}  ## function ends
	###########################
	sub update_exon_skipping_stats{
	my($a1Ref,$geneFileIndxRef,$statsHashRef) = @_;
	my(@temp);
		if($a1Ref->[$geneFileIndxRef->{"skippedExons"}]){
		@temp = split(',',$a1Ref->[$geneFileIndxRef->{"skippedExons"}]);
		$statsHashRef->{"skippedExons"}  += scalar(@temp);
		$statsHashRef->{"skipExonGenes"} += 1;
		}
	}  ## function ends
	###########################
	sub update_retained_intron_stats{
	my($a1Ref,$geneFileIndxRef,$statsHashRef) = @_;
	my(@temp);
		if($a1Ref->[$geneFileIndxRef->{"retainedIntrons"}]){
		@temp = split(',',$a1Ref->[$geneFileIndxRef->{"retainedIntrons"}]);
		$statsHashRef->{"retainedIntrons"} += scalar(@temp);
		$statsHashRef->{"retainIntronGenes"} += 1;
		}
	}  ## function ends
	###########################
	sub update_genic_regions_stats{
	my($a1Ref,$geneFileIndxRef,$statsHashRef) = @_;
	my(@temp);

		if($a1Ref->[$geneFileIndxRef->{"genicRegions"}]){
		@temp = split(',',$a1Ref->[$geneFileIndxRef->{"genicRegions"}]);
			if(scalar(@temp) > 1){
			$statsHashRef->{"MultipleRegions"} += 1;
			}
			else{
			$statsHashRef->{$temp[0]} += 1;
			}
		}  ## if(genic regions exist) ends
	}  ## function ends
	###########################
1;

