
package Begin::CreateOutputFiles;

use strict;
use warnings;

	sub characterization_out_file{
	my($fileTag) = @_;
	my %returnHash = (
		"MultiHits" => $fileTag."MultiHits.pslx2",
		"Discarded" => $fileTag."DiscardedReads.pslx2",
		"TopHits" => $fileTag."TopHits.pslx2",					

		"ExonOnly" => $fileTag."ExonOnly.out",				## Change it to ExonOnly
		"ExonDeletion" => $fileTag."ExonDeletion.out",
		"IntronOnly" => $fileTag."IntronOnly.out",				## Change it to IntronOnly
		"NeighboringExon" => $fileTag."NeighboringExon.out",  		## change it to NeighboringExons
		"MultipleAnnotations" => $fileTag."MultipleAnnotations.out",
		"GeneDesert" => $fileTag."GeneDesert.out", 				## change it to GeneDesert
		"InternalExonExtension" => $fileTag."InternalExonExtension.out",
		"GeneBoundaryChange" => $fileTag."UTRExpansion.out",		## Change it to GeneBoundaryExpansion
		"Uncharacterized" => $fileTag."UncharacterizedReads.out",

		## new ##
		"FullTranscriptMatch" => $fileTag."FullTranscriptMatch.out",  ## full transcript match file

		
		"BEDOutFile" => $fileTag."High-scoringBEDFile.out",  ## browser uploadable BED out file
		
		#### Merge stats files in to DataSheet ####
		"AlignmentStats" => $fileTag."AlignmentStats.out",  ## outfile for printing score/id/coverage bases classification of alignments
		"GeneAnnotationStats" => $fileTag."GeneAnnotationStats.out",  ## outfile for printing annotation statistics on different classes of reads characterized to


		"ExpressedGeneInfo" => $fileTag."ExpressedGeneInfo.out",  ## outfile for printing gene oriented annotation info. This info will be consolidated later  ## Temporary file (remove it later)
		"DGEProfile" => $fileTag."RPKMExpressionValues.out",  ## digital gene expression out file
		"NewIntronicExons" => $fileTag."NewIntronicExons.out",
		"ChimericTranscriptAnnotation" => $fileTag."ChimericTranscriptAnnotation.out",
		"DataSheet" => $fileTag."DataDistributionReport.out",		
		);
	create_file_from_hash(\%returnHash);  ## creating a new file	
	return(%returnHash);
	}  ## function ends
	###########################
	sub BED_out_file{  ## creating output files for BED file output of characterization files
	my($fileTag) = @_;
	my %returnHash = (
		"ExonOnly" => $fileTag."ExonOnlyBED.out",				## Change it to ExonOnly
		"ExonDeletion" => $fileTag."ExonDeletionBED.out",
		"IntronOnly" => $fileTag."IntronOnlyBED.out",				## Change it to IntronOnly
		"NeighboringExon" => $fileTag."NeighboringExonBED.out",  		## change it to NeighboringExons
		"MultipleAnnotations" => $fileTag."MultipleAnnotationBED.out",
		"GeneDesert" => $fileTag."GeneDesertBED.out", 				## change it to GeneDesert
		"InternalExonExtension" => $fileTag."InternalExonExtensionBED.out",
		"AlternativeTSS" => $fileTag."AlternativeTSSBED.out",
		"AlternativePolyA" => $fileTag."AlternativePolyABED.out",
		"Uncharacterized" => $fileTag."UncharacterizedReadsBED.out",

		## new ##		
		"FullTranscriptMatch" => $fileTag."FullTranscriptMatchBED.out",
		);
	return(%returnHash);
	}  ## function ends
	###########################
	sub pslx2_out_file{
	my($fileTag) = @_;
	my %returnHash = (
		"PSLX2Records" => $fileTag."Pslx2Files.pslx2",
		);
	create_file_from_hash(\%returnHash);  ## creating a new file
	return(%returnHash);
	}  ## function ends
	###########################
	sub pslx2_sorted_file{
	my($fileTag) = @_;
	my %returnHash = (
		"PSLX2SortedRecords" => $fileTag."Pslx2Sorted.pslx2",
		);
	create_file_from_hash(\%returnHash);  ## creating a new file
	return(%returnHash);
	}
	###########################
	sub pslx2_decision_making{
	my($fileTag) = @_;
	my %returnHash = (
		"PSLX2Annotate" => $fileTag."HighScoring.pslx2",
		"MultiHits" => $fileTag."MultiHits.pslx2",
		"Discarded" => $fileTag."DiscardedReads.pslx2",
		"TopHits" => $fileTag."TopHits.pslx2",
		"PotentialFusionAtLeft" => $fileTag."PotentialFusionAtLeft.pslx2",
		"PotentialFusionAtRight" => $fileTag."PotentialFusionAtRight.pslx2",
		"PotentialFusionAtBoth" => $fileTag."PotentialFusionAtBoth.pslx2",
		);

	create_file_from_hash(\%returnHash);
	return(%returnHash);
	}  ## function ends
	###########################
	sub annotate_out_files{
	my($fileTag)= @_;
	my %returnHash = (
		"ExonOnly" => $fileTag."ExonOnly.out",
		"ExonDeletion" => $fileTag."ExonDeletion.out",
		"IntronOnly" => $fileTag."IntronOnly.out",
		"NeighboringExon" => $fileTag."NeighboringExon.out",
		"SingleAnnotation" => $fileTag."SingleAnnotation.out",
		"MultipleAnnotations" => $fileTag."MultipleAnnotations.out",
		"GeneDesert" => $fileTag."GeneDesert.out",
		"InternalExonExtension" => $fileTag."InternalExonExtension.out",
		"GeneBoundaryChange" => $fileTag."GeneBoundaryChange.out",
		"Uncharacterized" => $fileTag."UncharacterizedReads.out",
		"AnnotationStats1" => $fileTag."AnnotationStats1.out",  ## outfile for printing annotation statistics on different classes of reads characterized to
		"AnnotationStats2" => $fileTag."AnnotationStats2.out",  ## outfile for printing gene oriented annotation info. This info will be consolidated later
	);
	create_file_from_hash(\%returnHash);
	return(%returnHash);
	}  ## function ends
	###########################
	sub chimer_out_files{
	my($fileTag) = @_;
	my %returnHash = (
		"Recombination" => $fileTag."Recombinations.out",
		"SameRecombination" => $fileTag."SameRecombinations.out",
		"MultiTopHits" => $fileTag."MultiTopHits.out",
		"MultiHits" => $fileTag."MultiHits.out",
		"NoRecombination" => $fileTag."NoRecombination.out",
		"RejectedAlignments" => $fileTag."RejectedAlignments.out",
		);
	create_file_from_hash(\%returnHash);
	return(%returnHash);
	}  ## function ends
	###########################
	sub chimer_annotation_files{
	my($fileTag) = @_;
	my %returnHash = (
		"GeneDesertChimera" => $fileTag."GeneDesertChimera.out",
		"Intragenic-type-II" => $fileTag."Intragenic-type-II.out",
		"Intragenic-type-I" => $fileTag."Intragenic-type-I.out",
		"Inter-genic" => $fileTag."Inter-genic.out",
		);
	create_file_from_hash(\%returnHash);
	return(%returnHash);
	}  ## function ends
	###########################
	sub pslx_verification_file{  ## outfile containing verified records
	my($fileTag) = @_;
	my %returnHash = (
		"VerifiedPSLRecords" => $fileTag."VerifiedPSLRecords.psl",
		"DiscardedPSLRecords" => $fileTag."DiscardedPSLRecords.psl",
		);
	create_file_from_hash(\%returnHash);
	return(%returnHash);
	}  ## function ends
	###########################
	sub default_align_file{
	my($fileTag) = @_;

	my %returnHash = (
		"alignFile" => $fileTag."DefaultAlignmentFile",
		);
	create_file_from_hash(\%returnHash);
	return(%returnHash);
	}  ## function ends
	###########################
	sub default_gene_file{
	my($fileTag) = @_;
	my %returnHash = (
		"DefaultGeneTrackFile" => $fileTag."DefaultGeneTrackFile.out",
		);
	create_file_from_hash(\%returnHash);
	return(%returnHash);
	}  ## function ends
	###########################
	sub create_file_from_hash{
	my($hashRef) = @_;
	my $key1;
		for $key1(keys %{$hashRef}){
		open FHCreateOutFile,">".$hashRef->{$key1};
		close FHCreateOutFile;
		}  ## foreach(keys(%{$hashRef})) ends
	}  ## function ends
	###########################
1;