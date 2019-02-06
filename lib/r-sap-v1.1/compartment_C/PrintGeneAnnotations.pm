package compartment_C::PrintGeneAnnotations;

use strict;
use warnings;

require compartment_C::FindGenicRegions;
require compartment_C::AnnotateExonDeletions;
require compartment_C::AnnotateExonExtensions;
require compartment_C::AnnotateExonsOnly;
require compartment_C::AnnotateIntronsOnly;
require compartment_C::AnnotateNeighboringGenes;
require compartment_C::AnnotateGeneBoundaryChanges;
require compartment_C::AnnotationProcessModules;

require Begin::Indexes;

	sub print_annotations{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$deletionCutoff,
	$exonExtensionCutoff,$annotStringency,$annotTagsRef,$outFilesRef) = @_;

	my($i,@outFileNames,@basicInfo,@annotationInfo,$annotationTagList,@printInfo);

	my($tempAnnotFlag,$tempAnnotTag1,$tempAnnotTag2,@tempArray);
	my(@geneMappingInfo1,@geneMappingInfo2);  ## will print the gene mapping info that will be used for the gene expression calculation
	my %expGeneInfoIndx = Begin::Indexes::expressed_gene_info();  ## hash for the expressed gene info array

	my $tempExonRegions = "";
	$tempAnnotFlag = 0;
	$tempAnnotTag1 = "";  ## key to outfile hash
	$tempAnnotTag2 = "";  ## actual annotation tag for printing in the file
	$annotationTagList = "";

	get_basic_align_annot_info($pslx2IndxRef,$a1Ref,\@basicInfo);
	get_basic_gene_annot_info($knGnIndxRef,$clustRef,$currTrackIndx,\@basicInfo);
	@geneMappingInfo1 = shared::ArrayOperations::initialize_value_array("NA",(scalar(keys(%expGeneInfoIndx))-2));


		if(!(scalar(@basicInfo))){  ## if no basic information was stored
		die "\n no basicInfo was stored for\n",join("\t",@{$a1Ref}),"\n",$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],"\n";
		}				

		for($i = 0;$i<scalar(@{$annotTagsRef});++$i){
		$tempAnnotTag2 = $annotTagsRef->[$i];

			if($annotTagsRef->[$i] eq "ExonOnly"){

			push(@outFileNames,$outFilesRef->{"ExonOnly"});
			@annotationInfo = compartment_C::AnnotateExonsOnly::get_exons_only_annot_info($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx);
			$tempExonRegions = $annotationInfo[scalar(@annotationInfo)-1];  ## exon regions(5'utr, cds, 3'utr or non-coding)
			##print "\ntempExonRegions: $tempExonRegions\n";

			update_gene_map_info1_array(\%expGeneInfoIndx,\@geneMappingInfo1,$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],
			    $clustRef->[$currTrackIndx][$knGnIndxRef->{"chr"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],$annotationInfo[scalar(@annotationInfo)-2],"NA","NA","NA","NA");
			}

			elsif($annotTagsRef->[$i] eq "ExonDeletion"){
			push(@outFileNames,$outFilesRef->{"ExonDeletion"});
			@annotationInfo = compartment_C::AnnotateExonDeletions::get_exon_deletion_annot_info($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$deletionCutoff);
			update_gene_map_info1_array(\%expGeneInfoIndx,\@geneMappingInfo1,$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],
			    $clustRef->[$currTrackIndx][$knGnIndxRef->{"chr"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],"NA","NA",$annotationInfo[scalar(@annotationInfo)-1],"NA",$annotationInfo[scalar(@annotationInfo)-3]);
			}

			elsif($annotTagsRef->[$i] eq "IntronOnly"){
			push(@outFileNames,$outFilesRef->{"IntronOnly"});
			@annotationInfo = compartment_C::AnnotateIntronsOnly::get_introns_only_annot_info($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx);
			update_gene_map_info1_array(\%expGeneInfoIndx,\@geneMappingInfo1,$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],
			    $clustRef->[$currTrackIndx][$knGnIndxRef->{"chr"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],"NA",$annotationInfo[scalar(@annotationInfo)-1],"NA","NA","NA");
			}

			elsif($annotTagsRef->[$i] eq "NeighboringExon"){
			push(@outFileNames,$outFilesRef->{"NeighboringExon"});
			@annotationInfo = compartment_C::AnnotateNeighboringGenes::get_neighboring_genes_annot_info($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx);
			
			## new code ##
			update_gene_map_info1_array(\%expGeneInfoIndx,\@geneMappingInfo1,$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],
			    $clustRef->[$currTrackIndx][$knGnIndxRef->{"chr"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],"NA","NA","NA","NA","NA");
			##

			}  ## for(scalar(@annotationTags) ends

			elsif($annotTagsRef->[$i] eq "GeneExpansion"){
			($tempAnnotFlag,$tempAnnotTag1,$tempAnnotTag2,@annotationInfo) = compartment_C::AnnotateGeneBoundaryChanges::get_gene_expansion_shrinkage_info
				($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$exonExtensionCutoff);
					if( ($tempAnnotFlag) && ($tempAnnotTag2 eq "AlternativeTSS") ){
					update_gene_map_info1_array(\%expGeneInfoIndx,\@geneMappingInfo1,$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],
					    $clustRef->[$currTrackIndx][$knGnIndxRef->{"chr"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],"NA","NA","NA","NA","5UTR");
					}
					elsif( ($tempAnnotFlag) && ($tempAnnotTag2 eq "AlternativePolyA") ){
					update_gene_map_info1_array(\%expGeneInfoIndx,\@geneMappingInfo1,$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],
					    $clustRef->[$currTrackIndx][$knGnIndxRef->{"chr"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],"NA","NA","NA","NA","3UTR");
					}
			push(@outFileNames,$outFilesRef->{$tempAnnotTag1});
			}  ## elsif(geneExpansion) ends
			
			elsif($annotTagsRef->[$i] eq "ExonExtension"){

			($tempAnnotFlag,$tempAnnotTag1,@annotationInfo)  = compartment_C::AnnotateExonExtensions::get_exon_extension_annot_info
						($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$exonExtensionCutoff);

			$tempAnnotTag2 = $tempAnnotTag1;  ## Annotation tag => key for out file hash
			update_gene_map_info1_array(\%expGeneInfoIndx,\@geneMappingInfo1,$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],
			    $clustRef->[$currTrackIndx][$knGnIndxRef->{"chr"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],"NA","NA","NA",$annotationInfo[scalar(@annotationInfo)-1],
			$annotationInfo[scalar(@annotationInfo)-3]);
			push(@outFileNames,$outFilesRef->{$tempAnnotTag1});
			}  ## elsif(InternalExonExtension) ends

			#### testing ####
			else{
			die "\n no annotation function define for $annotTagsRef->[$i] while printing annotation info\n";
			}
			if(!(scalar(@annotationInfo))){  ## if no annotation information was stored
			die "\n no annotationInfo was stored for annotationTag: $annotTagsRef->[$i]\n",join("\t",@{$a1Ref}),"\n",join("\t",@basicInfo),"\n";
			}				
			##################

		$annotationTagList .= $tempAnnotTag2.',';

			#if($tempAnnotTag2 eq "ExonOnly"){
			#push(@geneMappingInfo1,($clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],1,$annotationInfo[scalar(@annotationInfo)-2]));  ## incrementing number of reads mapped to  this gene + exons of the gene mapped to by sequencing reads
			#}

		unshift(@basicInfo,$tempAnnotTag2);
		push(@printInfo,(join("\t",(@basicInfo,@annotationInfo))));
		shift(@basicInfo);
		splice(@annotationInfo);
		}  ## for(@{$annotTagsRef}) ends


		## printing gene mapping stats info ##
		if($geneMappingInfo1[$expGeneInfoIndx{"transcriptID"}] ne "NA"){  ## if the array has been filled
		push(@geneMappingInfo2,($annotationTagList,$a1Ref->[$pslx2IndxRef->{"qID"}]));
		print_full_annot_info(\@geneMappingInfo1,\@geneMappingInfo2,$outFilesRef->{"ExpressedGeneInfo"});
		}
		shared::DeleteDataStructures::undefine_arrays(\@geneMappingInfo1,@geneMappingInfo2);
		#######################################

		@tempArray  = split(',',$annotationTagList);  ## all the collected annotation tags

		if(scalar(@tempArray) == 1){  ## if only one annotation tag
		print_array_to_file(\@printInfo,$outFileNames[0]);
		return($tempArray[0],$tempExonRegions);  ## returning final annotation tag
		}
		elsif(scalar(@tempArray) > 1){  ## if more than one annotation tags
		print_array_to_file(\@printInfo,$outFilesRef->{"MultipleAnnotations"});
		return("MultipleAnnotations",$tempExonRegions);
		}
	shared::DeleteDataStructures::undefine_arrays(\@printInfo,\@basicInfo,\@geneMappingInfo1,\@geneMappingInfo2,\@annotationInfo);
	}  ## function ends
	###########################
	sub get_basic_align_annot_info{
	my($pslx2IndxRef,$a1Ref,$basicInfoRef) = @_;
	push(@{$basicInfoRef},$a1Ref->[$pslx2IndxRef->{"qID"}]);
	push(@{$basicInfoRef},$a1Ref->[$pslx2IndxRef->{"identity"}]);
	push(@{$basicInfoRef},$a1Ref->[$pslx2IndxRef->{"strand"}]);
	push(@{$basicInfoRef},$a1Ref->[$pslx2IndxRef->{"tID"}]);
	push(@{$basicInfoRef},$a1Ref->[$pslx2IndxRef->{"tSt"}]);
	push(@{$basicInfoRef},$a1Ref->[$pslx2IndxRef->{"tEn"}]);
	}  ## function ends
	###########################
	sub get_basic_gene_annot_info{
	my($knGnIndxRef,$clustRef,$currTrackIndx,$basicInfoRef) = @_;
	push(@{$basicInfoRef},$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}]);
	push(@{$basicInfoRef},$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]);
	push(@{$basicInfoRef},compartment_C::FindGenicRegions::check_if_protein_coding($clustRef,$currTrackIndx,$knGnIndxRef));
	}  ## function ends
	###########################
	sub print_full_annot_info{
	my($arr1Ref,$arr2Ref,$outFileName) = @_;
	open FHPrintAnnotation,">>$outFileName";  ## output file
	print FHPrintAnnotation join("\t",@{$arr1Ref}),"\t",join("\t",@{$arr2Ref}),"\n";
	close FHPrintAnnotation;
	}  ## function ends
	###########################
	sub print_array_to_file{
	my($printArrayRef,$outFileName) = @_;
	open FHPrintAnnotation,">>$outFileName";  ## output file
	print FHPrintAnnotation join("\n",@{$printArrayRef}),"\n";
	close FHPrintAnnotation;
	}  ## function ends
	###########################
	sub update_gene_map_info1_array{
	my($expGeneIndxRef,$geneInfoRef,@mapInfo) = @_;
	
	$geneInfoRef->[$expGeneIndxRef->{"transcriptID"}] = $mapInfo[0];
	$geneInfoRef->[$expGeneIndxRef->{"chr"}] = $mapInfo[1];
	$geneInfoRef->[$expGeneIndxRef->{"txSt"}] = $mapInfo[2];



	$geneInfoRef->[$expGeneIndxRef->{"spannedExons"}] = $mapInfo[3];
	$geneInfoRef->[$expGeneIndxRef->{"spannedIntrons"}] = $mapInfo[4];

	
		if($geneInfoRef->[$expGeneIndxRef->{"skippedExons"}] eq "NA"){
		$geneInfoRef->[$expGeneIndxRef->{"skippedExons"}] = $mapInfo[5];
		}
		else{
		$geneInfoRef->[$expGeneIndxRef->{"skippedExons"}] .= ','.$mapInfo[5];
		}


		if($geneInfoRef->[$expGeneIndxRef->{"retainedIntrons"}]  eq "NA"){
		$geneInfoRef->[$expGeneIndxRef->{"retainedIntrons"}]  = $mapInfo[6];
		}
		else{
		$geneInfoRef->[$expGeneIndxRef->{"retainedIntrons"}]  .= ','.$mapInfo[6];
		}

		if($geneInfoRef->[$expGeneIndxRef->{"genicRegions"}] eq "NA"){
		$geneInfoRef->[$expGeneIndxRef->{"genicRegions"}] = $mapInfo[7];
		}
		else{
		$geneInfoRef->[$expGeneIndxRef->{"genicRegions"}] .= ','.$mapInfo[7];
		}
	}  ## function ends
	###########################
1;

	
##Basic read-alingmnent related info: Read ID, length, alignment score, alignment identity, 
	## coverage, chromosome, strand, gene ID, strand

##Annotation specific info:
	## exons only:  number of exons, which exons (number), regions: (coding, non-coding, utrs).
	
	## exon deletions: number of base pairs deleted, which exons, number of complete exon deletetions, 
		## which exons were deleted completely, total number of base pairs deleted, 
		## % the transcript length deleted, region: (coding, non-coding, utrs)
		
		
	## exon extension: number of bases extended, which introns, number of complete introns inclusions,
		## which introns were completely included, total extension, %transcript extended
	
	## introns only: which intron, (region: coding, non-coding, utrs), spliced?
	
	## gene expansion: Which side, how much expansion of the gen boundary, spliced?
	
	## neighboring genes: Which gene got expanded, which side, how much is the distance
		## between the closest ends of the read and the expanded gene
	



			#### paste temporary code here ####
			###################################

