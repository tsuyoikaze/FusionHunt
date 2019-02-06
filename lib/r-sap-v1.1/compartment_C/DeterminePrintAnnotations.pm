package compartment_C::DeterminePrintAnnotations;

use strict;
use warnings;

require compartment_C::Pslx2RecordProcessing;
require compartment_C::FindGenicRegions;
require compartment_C::PrintGeneAnnotations;
require compartment_C::PrintUncharacterizedReads;

require compartment_C::AnnotationProcessModules;
require compartment_C::DetermineCharacteristicsModules;
require compartment_C::SortCharacteristicsMatrix;
require compartment_C::PrintExonsIntronsOnlyGeneInfo;
require compartment_C::AnnotateExonsOnly;

require file_parsers::PSLX2ToBED;

	sub determine_print_annotations{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$cutoffHashRef,$outFilesRef,$bedFilesRef,$annotStatsHashRef,$exonsOnlyHashRef) = @_;

	my($n,$str1,$str2,$i,$j);

	my(@temp,@charArray);
	my(@ExB,@ExE);	####for alignment blocks( on read)
	my(@GExB,@GExE); #### for exons on reference transcript
	my($cdsSt,$cdsEn);  ## for coding region coordinates
	my(@inSt,@inEn);  ## for introns reference transcript
	my($cdsFlag,$geneExpansionFlag,$end5Expansion,$end3Expansion,$expansion);
	my(@annotationTags,$topReplicates);
	my($tempTag,$tempStr,@finalAnnotTags,@exonRegions);  ## values to be obtained from print function
	my(%charArrayFields);
	my(@bed);

	$n = scalar(@{$clustRef});
	

	compartment_C::Pslx2RecordProcessing::get_pslx2_algn_blocks($pslx2IndxRef,$a1Ref,\@ExB,\@ExE);
	compartment_C::AnnotationProcessModules::define_char_array_field_hash(\%charArrayFields);  ## initializing temporary array


		#### determining annotation characteristics ####
		for($i = 0;$i<$n;++$i){  ## for all the intersected genes
		##print "\n Determining characteristics for ",$clustRef->[$i][$knGnIndxRef->{"transcriptID"}],"\n";
		shared::ArrayOperations::initialize_array(\@temp,scalar(keys(%charArrayFields)));  ## initializing a temporary array corresponding to charArrayHash (all values are 0 or NULL)

		$temp[$charArrayFields{"clustIndex"}] = $i;  ## index of the gene track in the cluster (all the possible transcripts/genes for the given read)
		$temp[$charArrayFields{"transcriptID"}] = $clustRef->[$i][$knGnIndxRef->{"transcriptID"}];  ## pushing UCSC ID

		compartment_C::FindGenicRegions::get_exon_cords($clustRef,$i,$knGnIndxRef,\@GExB,\@GExE);
		compartment_C::FindGenicRegions::get_intron_regions(\@GExB,\@GExE,\@inSt,\@inEn);
	
		# new #
		$temp[$charArrayFields{"intronChainMatch"}] = compartment_C::DetermineCharacteristicsModules::find_intron_chain_match(\@ExB,\@ExE,\@GExB,\@GExE);  ## intron chain match
		########

		$temp[$charArrayFields{"exonOverlap"}] = shared::CoordinateComparisons::find_block_overlaps(\@ExB,\@ExE,\@GExB,\@GExE);  ## exon overlap
		$temp[$charArrayFields{"geneOverlap"}] = compartment_C::DetermineCharacteristicsModules::find_gene_overlap(\@ExB,\@ExE,$clustRef->[$i][$knGnIndxRef->{"txSt"}],$clustRef->[$i][$knGnIndxRef->{"txEn"}]);  ## gene overlap
		$temp[$charArrayFields{"exonDeletions"}] = compartment_C::DetermineCharacteristicsModules::find_exon_deletions(\@ExB,\@ExE,\@GExB,\@GExE,$cutoffHashRef->{"DeletionCutoff"});  ## exon deletions

		($temp[$charArrayFields{"intronOverlap"}],$temp[$charArrayFields{"trueIntronOverlapFlag"}]) = 
				compartment_C::DetermineCharacteristicsModules::find_intron_overlaps(\@ExB,\@ExE,\@inSt,\@inEn,$cutoffHashRef->{"ExonExtensionCutoff"});  ## intron overlap

		($temp[$charArrayFields{"geneExpansion"}], $temp[$charArrayFields{"trueGeneExpansionFlag"}]) = 
			compartment_C::DetermineCharacteristicsModules::find_gene_expansion_flags($a1Ref->[$pslx2IndxRef->{"tSt"}],$a1Ref->[$pslx2IndxRef->{"tEn"}],
						$clustRef->[$i][$knGnIndxRef->{"txSt"}],$clustRef->[$i][$knGnIndxRef->{"txEn"}],
						,$cutoffHashRef->{"ExonExtensionCutoff"});
		
		$temp[$charArrayFields{"geneDistance"}] = shared::CoordinateComparisons::calculate_end_to_end_distance(
						$a1Ref->[$pslx2IndxRef->{"tSt"}],$a1Ref->[$pslx2IndxRef->{"tEn"}],
						$clustRef->[$i][$knGnIndxRef->{"txSt"}],$clustRef->[$i][$knGnIndxRef->{"txEn"}]);  ## distance of the sequencing read from the gene (if read lies completely outside the gene)
								
		compartment_C::FindGenicRegions::get_cds_pos($clustRef,$i,$knGnIndxRef,\$cdsSt,\$cdsEn);
		$temp[$charArrayFields{"codingFlag"}] = ($cdsEn-$cdsSt);  ## size of the coding region(0 if non-coding)
		push(@charArray,[@temp]);  ## pushing to collect info on all intersected transcripts

		shared::DeleteDataStructures::undefine_arrays(\@temp,\@GExB,\@GExE,\@inSt,\@inEn);
		}  ## for(all intersected transcripts($n)) ends

		#### sorting characteristics array to find best transcript to for the characterization ####
		if($cutoffHashRef->{"AnnotationMode"} eq "unique"){
		compartment_C::SortCharacteristicsMatrix::sort_char_array(\%charArrayFields,\@charArray);  ## sorting on different characteristics determined above
		$topReplicates = compartment_C::AnnotationProcessModules::find_top_replicates(\%charArrayFields,\@charArray);
		compartment_C::AnnotationProcessModules::sort_top_gene_ids($topReplicates,\%charArrayFields,\@charArray);
		}  ## sorting based on different annotation characterictics
		################################################

		
		############## Substitution code starts here ##############
		if($cutoffHashRef->{"AnnotationMode"} eq "unique"){
		my $intronsOnlyFlag = 0;

			for($i = 0;$i<scalar(@charArray);++$i){  ## for all the known mapped transcripts
			compartment_C::AnnotationProcessModules::decide_annotation_tag(\%charArrayFields,\@charArray,$i,\@annotationTags);


				if(!(scalar(@annotationTags))){  ## if no annotation tag was defined (## change it to uncharacterized).
				##die "\n no annotation tag was defined for pslx record\n", join("\t",@{$a1Ref}),"\n and gene ",$charArray[0][$charArrayFields{"transcriptID"}],"\n";
				next;
				}


				if($i == 0){
					if(compartment_C::PrintUncharacterizedReads::check_if_artifact($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,
						$charArray[$i][$charArrayFields{"clustIndex"}],$cutoffHashRef->{"IntronCutoff"},
		        			$cutoffHashRef->{"AnnotationMode"},\@annotationTags,$outFilesRef)){
					push(@finalAnnotTags,"Uncharacterized");
					## do nothing: Already printed to the file
					next;
					}  ## if(uncharacterized) ends

					else{
					($tempTag,$tempStr) = compartment_C::PrintGeneAnnotations::print_annotations($pslx2IndxRef,$knGnIndxRef,$a1Ref,
						$clustRef,$charArray[$i][$charArrayFields{"clustIndex"}],
						    $cutoffHashRef->{"DeletionCutoff"},$cutoffHashRef->{"ExonExtensionCutoff"},
						    $cutoffHashRef->{"AnnotationMode"},\@annotationTags,$outFilesRef);
					push(@finalAnnotTags,$tempTag);
					push(@exonRegions,$tempStr);

						if($annotationTags[0] eq "ExonOnly"){
						compartment_C::PrintExonsIntronsOnlyGeneInfo::check_print_full_transcript_info($pslx2IndxRef,
						$knGnIndxRef,$a1Ref,$clustRef,$charArray[$i][$charArrayFields{"clustIndex"}],$cutoffHashRef->{"ExonExtensionCutoff"},
						$outFilesRef,$bedFilesRef);					
						}  ## printing if full transcript match						

						if($annotationTags[0] eq "IntronOnly"){
						$intronsOnlyFlag = 1;
						}

					}  ## elsif(!(uncharacterized)) ends
				}  ## if(i == 0) ends


				else{  ## if(i != 0);  ## printing gene expression info for rest of the genes
					if($annotationTags[0] eq "ExonOnly"){
					compartment_C::PrintExonsIntronsOnlyGeneInfo::print_exons_introns_only_gene_info($pslx2IndxRef,
						$knGnIndxRef,$a1Ref,$clustRef,$charArray[$i][$charArrayFields{"clustIndex"}],
						\@annotationTags,$outFilesRef);
					}  ## if(annotationTag == ExonOnly) ends

					elsif( ($annotationTags[0] eq "IntronOnly") && ($intronsOnlyFlag) ){
					compartment_C::PrintExonsIntronsOnlyGeneInfo::print_exons_introns_only_gene_info($pslx2IndxRef,
						$knGnIndxRef,$a1Ref,$clustRef,$charArray[$i][$charArrayFields{"clustIndex"}],
						\@annotationTags,$outFilesRef);
					}  ## elsif($annotationTag == IntronOnly && IntronOnlyFlag == 1) ends
				}  ## elsef(i != 0) ends

			splice(@annotationTags);
			}  ## for(scalar(@charArray)) ends
		}  ## if(annotationStirngency == unique) ends


		if($cutoffHashRef->{"AnnotationMode"} eq "multi"){

			for($i = 0;$i<scalar(@charArray);++$i){
			compartment_C::AnnotationProcessModules::decide_annotation_tag(\%charArrayFields,\@charArray,$i,\@annotationTags);

				if(!(scalar(@annotationTags))){  ## if no annotation tag was defined
				##die "\n no annotation tag was defined for pslx record\n", join("\t",@{$a1Ref}),"\n and gene ",$charArray[0][$charArrayFields{"transcriptID"}],"\n";
				next;
				}

				if(compartment_C::PrintUncharacterizedReads::check_if_artifact($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,
					$charArray[$i][$charArrayFields{"clustIndex"}],$cutoffHashRef->{"IntronCutoff"},
	        			$cutoffHashRef->{"AnnotationMode"},\@annotationTags,$outFilesRef)){
				push(@finalAnnotTags,"Uncharacterized");
				## do nothing: Already printed to the file
				next;
				}  ## if(uncharacterized) ends

				else{
				($tempTag,$tempStr) = compartment_C::PrintGeneAnnotations::print_annotations($pslx2IndxRef,$knGnIndxRef,$a1Ref,
					$clustRef,$charArray[$i][$charArrayFields{"clustIndex"}],
					    $cutoffHashRef->{"DeletionCutoff"},$cutoffHashRef->{"ExonExtensionCutoff"},
					    $cutoffHashRef->{"AnnotationMode"},\@annotationTags,$outFilesRef);
				push(@finalAnnotTags,$tempTag);
				push(@exonRegions,$tempStr);

					if($annotationTags[0] eq "ExonOnly"){
					compartment_C::PrintExonsIntronsOnlyGeneInfo::check_print_full_transcript_info($pslx2IndxRef,
					$knGnIndxRef,$a1Ref,$clustRef,$charArray[$i][$charArrayFields{"clustIndex"}],$cutoffHashRef->{"ExonExtensionCutoff"},
					$outFilesRef,$bedFilesRef);					
					}  ## printing if full transcript match						

				}  ## elsif(!(uncharacterized)) ends
			splice(@annotationTags);
			}  ## for(scalar(@charArray)) ends
		}  ## if(annotationMode == multi) ends
		


		my $bedFileTag;

		if(  (join('',@finalAnnotTags) =~ m/(multi)/i) || (scalar(@finalAnnotTags) > 1) ){  ## if "MultipleAnnotations" tag found
		++$annotStatsHashRef->{"MultipleAnnotations"};
		$bedFileTag = "MultipleAnnotations";
		}
		elsif(scalar(@finalAnnotTags) == 1){  ## only one annotation
		++$annotStatsHashRef->{$finalAnnotTags[0]};
		$bedFileTag = $finalAnnotTags[0];
		}
		elsif(!(scalar(@finalAnnotTags))){  ## if no annotation tag was defined for this read
		my @tempClust;
		push(@tempClust,[@{$a1Ref}]);
		shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{"Uncharacterized"},scalar(keys(%{$pslx2IndxRef})),1,\@tempClust);
		splice(@tempClust);
		++$annotStatsHashRef->{"Uncharacterized"};
		$bedFileTag = "Uncharacterized";
		}

		#### printing to BED file parts ####
		

			if(!(exists($bedFilesRef->{$bedFileTag}))){
			die "\n no bed file was found for annotation tag $bedFileTag \n";
			}

		open FHPrintBED,">>".$bedFilesRef->{$bedFileTag};
		@bed = file_parsers::PSLX2ToBED::pslx2_to_bed($pslx2IndxRef,@{$a1Ref});
		print FHPrintBED join("\t",@bed),"\n";  ## send the full array
		close FHPrintBED;
		#####################################
	

		#### update exons only regions mapping stat hash ####
		compartment_C::AnnotateExonsOnly::update_exons_only_hash(\@exonRegions,$exonsOnlyHashRef);
		############### Substitution code ends here ###############


	  ## paste temporary code here ##
	  ###############################


	shared::DeleteDataStructures::undefine_arrays(\@charArray,\@annotationTags,\@finalAnnotTags,\@GExE,\@GExE,\@ExB,\@ExE);  ## undefining the data types
	shared::DeleteDataStructures::undefine_hashes(\%charArrayFields);
	shared::DeleteDataStructures::undefine_scalars(\$cdsSt,\$cdsEn);
	}  ## function ends
	###########################
1;

=for  ## multiline commented code begins here
	#######
	new char_array with different cheracteristics
	for (each transcript){
	*find_exon_coverage -> @char_array;
	*find_gene_overlap -> @char_array;
	find_exon_skipping -> @char_array; ## included in exon deletions
	*find_internal_exon_deletion -> @char_array;
	*find_intron_overlap -> @char_array;
	*find_gene_expansion -> @char_array;
	*find_coding_noncoding_flag -> @char_array;
	
	sort(
	exon_coverage: maximize
	gene_overlap: maximize
	exon_deletion: minimize
	intron_overlap: minimize
	utr_expansion: minimize
	coding_noncoding: choose coding
	}
	
	if(more than one similar top tags){
	sort_on_transcriptID()
	pick_top_one
	}
		
	print(annotation tags)
	print(annotation details);
	}
		#### to be removed ####
		print "\n\n";
		for($i = 0;$i<$topReplicates;++$i){
			for($j = 0;$j<scalar(keys(%charArrayFields));++$j){
			print "$charArray[$i][$j]\t";
			}
		print "\n";
		}  ## for($i) printing ends
		########################
=cut  ## commented code ends here

		##########################################
		### temporary code to be inseted here ####
