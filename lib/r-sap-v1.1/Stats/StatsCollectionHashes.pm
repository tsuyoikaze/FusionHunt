package Stats::StatsCollectionHashes;

use strict;
use warnings;

	sub alignment_class_hash{
	my %returnHash = (
		"HighScoring" => 0,
		"MultiHits" => 0,
		"Discarded" => 0,
		"Chimers" => 0,
		);
	return(%returnHash);
	}  ## function ends
	###########################
	sub annotation_class_hash{
	my %returnHash = (
		"ExonOnly" => 0,
		"ExonDeletion" => 0,
		"IntronOnly" => 0,
		"NeighboringExon" => 0,
		"MultipleAnnotations" => 0,
		"GeneDesert" => 0,
		"InternalExonExtension" => 0,
		"AlternativePolyA" => 0,
		"AlternativeTSS" => 0,
		"MultipleAnnotations" => 0,  ## reads with multiple annotation tags
		"Uncharacterized" => 0,  ## uncharacterized reads (potential artifacts in the data)
		);
	return(%returnHash);
	}  ## function ends
	###########################
	sub chimer_distribution_hash{
	my %returnHash = (
		"GeneDesertChimera" => 0,
		"ExonReorderingChimera" => 0,
		"ExonFlippingChimera" => 0,
		"IntergenicChimera" => 0,
		);
	return(%returnHash);
	}  ## function ends
	###########################
	sub gene_mapping_distribution_hash{
	my %geneDistHash = (
		"detectedGenes" => 0,
		"normallySpliced" => 0,
		"oneAbSplicing" => 0,
		"multiAbSplicing" => 0,
		"codingDetected" => 0,
		"non-codingDetected" => 0,
		);
	return(%geneDistHash);
	}  ## function ends
	###########################
	sub exons_only_read_hash{
	my %returnHash = (
		"5UTR" => 0,
		"CDS" => 0,
		"3UTR" => 0,
		"NonCoding" => 0,
		"MultiMapping" => 0,
		);
	return(%returnHash);
	}  ## function ends
	###########################
	sub exon_skip_hash{
	my %returnHash = (
		"skippedExons" => 0,
		"skipExonGenes" => 0,
		"skippedExonReads" => 0,
		);
	return(%returnHash);
	}  ## function end
	###########################
	sub retained_intron_hash{
	my %returnHash =(
		"retainedIntrons" => 0,
		"retainIntronGenes" => 0,
		"retainedIntronReads" => 0,
		);
	return(%returnHash);
	}  ## function ends
	###########################
	sub spliced_genic_regions_hash{
	my %returnHash = (
		"5UTR" => 0,
		"CDS" => 0,
		"3UTR" => 0,
		"NonCoding" => 0,
		"MultipleRegions" => 0,
		);
	return(%returnHash);
	}  ## function ends
	###########################
1;
