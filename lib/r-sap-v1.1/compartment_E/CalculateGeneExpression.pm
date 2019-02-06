
package compartment_E::CalculateGeneExpression;

use strict;
use warnings;

require compartment_C::FindGenicRegions;
require compartment_E::CalculateGeneExpressionSR;
require compartment_E::UpdateGeneStats;

require Begin::Indexes;

	sub calculate_expression{
	my($totalReads,$geneFile,$outFilesRef,$knGnIndxRef,$geneCharHashRef,$geneDistHashRef,
		$exonSkippingStatsRef,$retainedIntronStatsRef,$splicedRegionsStatsRef) = @_;

	my %geneFileIndx = Begin::Indexes::ref_gene_array_info();  ## array structure to store refgene info

	my %readFileIndx = Begin::Indexes::expressed_gene_info();  ## array structure to read in expressed gene info file

	## reading the reference gene file ##
	my %geneInfoHash = compartment_E::CalculateGeneExpressionSR::collect_ref_gene_info($geneFile,$knGnIndxRef);

	my($str,@a1,@temp1,@temp2,$i);
	my($gExp);  ## gene expression value


	#### Consolidating expressed/detected gene info from output file ####
	open FHDGEIn,$outFilesRef->{"ExpressedGeneInfo"} or die "\n can not open file ",$outFilesRef->{"ExpressedGeneInfo"},"\n";  ## gene read characterization file

		while($str = <FHDGEIn>){  ## note: format of the gene stats file is fixed
		$str =~ s/\n//;
		$str =~ s/\r//;
		@a1 = split(/\t+/,$str);  ## read characterization info (explained by %readFileIndx)

		my $key = $a1[$readFileIndx{"transcriptID"}]."____".$a1[$readFileIndx{"chr"}]."____".$a1[$readFileIndx{"txSt"}];  ## unique transcriptID (with chromosome and transcription start position)

			if(exists $geneInfoHash{$key}){  ## if key exists (transcriptID)
			$geneInfoHash{$key}->[$geneFileIndx{"readCount"}] += get_curr_read_count(\@a1,\%readFileIndx);  ## updating read counts (spanned exons)

				## updating the exon mapped read counts
				if($a1[$readFileIndx{"spannedExons"}] !~ m/(NA)/){  ## if exon counts exist
				@temp1 = split(',',$a1[$readFileIndx{"spannedExons"}]);  ## exons of the transcript mapped by sequencing reads
				@temp1 = shared::ArrayOperations::remove_array_redundancy(@temp1);
				@temp2 = split(',',$geneInfoHash{$key}->[$geneFileIndx{"exonReadCounts"}]);  ## read count for each exon
				compartment_E::CalculateGeneExpressionSR::increment_read_counts(\@temp1,\@temp2);
				$geneInfoHash{$key}->[$geneFileIndx{"exonReadCounts"}] = join(',',@temp2);  ## replacing by the updated read counts for each exon
				}
				
				## updating the intron mapped read counts
				if($a1[$readFileIndx{"spannedIntrons"}] !~ m/(NA)/){  ## if intron counts exist
				@temp1 = split(',',$a1[$readFileIndx{"spannedIntrons"}]);  ## exons of the transcript mapped by sequencing reads
				@temp1 = shared::ArrayOperations::remove_array_redundancy(@temp1);
				@temp2 = split(',',$geneInfoHash{$key}->[$geneFileIndx{"intronReadCounts"}]);  ## read count for each exon
				compartment_E::CalculateGeneExpressionSR::increment_read_counts(\@temp1,\@temp2);
				$geneInfoHash{$key}->[$geneFileIndx{"intronReadCounts"}] = join(',',@temp2);  ## replacing by the updated read counts for each exon
				}
				
				## updating the skipped exons
				if($a1[$readFileIndx{"skippedExons"}] !~ m/(NA)/){  ## if there are skipepd exons
				$exonSkippingStatsRef->{"skippedExonReads"} += 1;			
				$geneInfoHash{$key}->[$geneFileIndx{"skippedExons"}] = update_list($a1[$readFileIndx{"skippedExons"}],$geneInfoHash{$key}->[$geneFileIndx{"skippedExons"}]);
				}
				
				if($a1[$readFileIndx{"retainedIntrons"}] !~ m/(NA)/){  ## if there are retained introns
				$retainedIntronStatsRef->{"retainedIntronReads"} += 1;
				$geneInfoHash{$key}->[$geneFileIndx{"retainedIntrons"}] = update_list($a1[$readFileIndx{"retainedIntrons"}],$geneInfoHash{$key}->[$geneFileIndx{"retainedIntrons"}]);
				}
				
				if($a1[$readFileIndx{"genicRegions"}] !~ m/(NA)/){  ## if there are spliced genic regions
				$geneInfoHash{$key}->[$geneFileIndx{"genicRegions"}] = update_list($a1[$readFileIndx{"genicRegions"}],$geneInfoHash{$key}->[$geneFileIndx{"genicRegions"}]);
				}
			

				my $tempMultiFlag = compartment_E::CalculateGeneExpressionSR::check_if_multi_splicing($a1[$readFileIndx{"annotTags"}]);


				if($tempMultiFlag > 0){
				$geneInfoHash{$key}->[$geneFileIndx{"multiAnnotFlag"}] += $tempMultiFlag;
				}
				else{
				$geneInfoHash{$key}->[$geneFileIndx{"annotTags"}] = update_list($a1[$readFileIndx{"annotTags"}],$geneInfoHash{$key}->[$geneFileIndx{"annotTags"}]);
				}

			}  ## if(key) exists ends

			else{
			die "\n no geneInfo in hash for ",$a1[$readFileIndx{"transcriptID"}],"\n";
			}
		}  ## while(<FHDGEIn>) ends
		close FHDGEIn;
	###############################################################



	#### gene expression estimation and stats collection/update stats ####
	open FHDGEOut,">>".$outFilesRef->{"DGEProfile"};  ## output file1 (for gene expression)
	open FHPrintIntrons,">>".$outFilesRef->{"NewIntronicExons"};  ## output file
	splice(@temp1);
	splice(@temp2);
	my($tempExp);

		#### gene/transcript/exon expression calculation ####
		for my $key(keys(%geneInfoHash)){
		$gExp = compartment_E::CalculateGeneExpressionSR::calculate_RPKM($geneInfoHash{$key}->[$geneFileIndx{"transcriptLen"}],$geneInfoHash{$key}->[$geneFileIndx{"readCount"}],$totalReads);

			if($gExp > 0){
			print FHDGEOut "$geneInfoHash{$key}->[0]\t$geneInfoHash{$key}->[1]\t$geneInfoHash{$key}->[2]\t$geneInfoHash{$key}->[3]\t$gExp\n";
			}
		@temp1 = split(',',$geneInfoHash{$key}->[$geneFileIndx{"intronReadCounts"}]);
		@temp2 = split(',',$geneInfoHash{$key}->[$geneFileIndx{"intronLengths"}]);

			for(my $j = 0;$j<scalar(@temp1);++$j){
			$tempExp = 0;
			$tempExp = compartment_E::CalculateGeneExpressionSR::calculate_RPKM($temp2[$j],$temp1[$j],$totalReads);
				if($gExp){
					if( (compartment_E::CalculateGeneExpressionSR::check_if_same_order($gExp,$tempExp)) ){
					print FHPrintIntrons "$geneInfoHash{$key}->[0]\t$geneInfoHash{$key}->[1]\t$geneInfoHash{$key}->[2]\t$geneInfoHash{$key}->[3]\t",($j+1),"\t$tempExp\t$gExp\n";
					}
				}  ## if($gExp) ends
			}  ## for(all the introns) ends

		@temp1 = @{$geneInfoHash{$key}};  ## expressed gene info



		## updating stats hash ##
		compartment_E::UpdateGeneStats::update_gene_dist_stats($key,\@temp1,\%geneFileIndx,$geneCharHashRef,$geneDistHashRef);

		compartment_E::UpdateGeneStats::update_exon_skipping_stats(\@temp1,\%geneFileIndx,$exonSkippingStatsRef);
		compartment_E::UpdateGeneStats::update_retained_intron_stats(\@temp1,\%geneFileIndx,$retainedIntronStatsRef);
		compartment_E::UpdateGeneStats::update_genic_regions_stats(\@temp1,\%geneFileIndx,$splicedRegionsStatsRef);
		####################
		splice(@temp1);
		}  ## for(geneInfoHash) ends

	close FHDGEOut;
	close FHPrintIntrons;
	}  ## function ends
	###########################
	sub update_list{
	my($list1,$list2) = @_;
	my(@temp1,@temp2);	
	@temp1 = split(',',$list1);  ## list is comma separated
	@temp2 = split(',',$list2);  ## list is comma separated
	push(@temp2,@temp1);
	@temp2 = shared::ArrayOperations::remove_array_redundancy(@temp2);
	return(join(',',@temp2));
	}  ## function ends
	###########################
	sub get_curr_read_count{
	my($a1Ref,$readFileIndxRef) = @_;
	my($readCount,@temp);
	$readCount = 0;
		if($a1Ref->[$readFileIndxRef->{"spannedExons"}] !~ m/(NA)/){  ## if exon counts exist
		$readCount = 1;
		}
	return($readCount);
	}  ## function ends
	###########################
	1;




