package compartment_D::AnnotateChimers;

use strict;
use warnings;

use compartment_D::ChimerGeneAnnotationSR;
use compartment_D::ChimerGeneFindSR;


	sub annotate_chimeric_transcripts{
	my($chimerIndxRef,$chimerArrayRef,$knGnIndxRef,$geneFileIndxRef,$geneFileRef,$outFilesRef,$chimerReadStatsRef) = @_;

	my %fragGeneCharsIndexes = compartment_D::ChimerGeneAnnotationSR::define_chimer_char_hash();
	
	my(@combChrs,@breakPoints);
	my(@fragment1Genes,@fragment2Genes);  ## matrices
	my($topReplicates);
	my($geneIndex1,$geneIndex2);  ## index for the chosen genes from fragmentGenes array
	my($chimerTag);  ## determines structure of the chimeric transcript
	my($foutName);  ## output file nam (will be decided according to the chimer type)
	my $i;


	## to be removed later (chimer prediction will do it already)
	compartment_D::ChimerGeneAnnotationSR::flip_fusion_record($chimerIndxRef,$chimerArrayRef);  ## converting to right fill if left fill already

	compartment_D::ChimerGeneAnnotationSR::find_connected_chrs($chimerIndxRef,\@combChrs,$chimerArrayRef);  ## getting fused chromosomes

	compartment_D::ChimerGeneAnnotationSR::find_break_points($chimerIndxRef,\@breakPoints,$chimerArrayRef);
	
	compartment_D::ChimerGeneAnnotationSR::fragment_gene_finder($chimerIndxRef,$knGnIndxRef,$combChrs[0],$breakPoints[0],$chimerArrayRef,$geneFileIndxRef,$geneFileRef,\@fragment1Genes);
	compartment_D::ChimerGeneAnnotationSR::fragment_gene_finder($chimerIndxRef,$knGnIndxRef,$combChrs[1],$breakPoints[1],$chimerArrayRef,$geneFileIndxRef,$geneFileRef,\@fragment2Genes);
	
	compartment_D::ChimerGeneAnnotationSR::sort_frag_genes(\%fragGeneCharsIndexes,\@fragment1Genes);  ## sorting on different attributes
	compartment_D::ChimerGeneAnnotationSR::sort_frag_genes(\%fragGeneCharsIndexes,\@fragment2Genes);  ## sorting on different attributes
	
	$topReplicates = compartment_D::ChimerGeneAnnotationSR::find_top_chimer_genes(\%fragGeneCharsIndexes,\@fragment1Genes);

	compartment_D::ChimerGeneAnnotationSR::sort_top_gene_ids($topReplicates,\%fragGeneCharsIndexes,\@fragment1Genes);
	
	$topReplicates = compartment_D::ChimerGeneAnnotationSR::find_top_chimer_genes(\%fragGeneCharsIndexes,\@fragment2Genes);
	compartment_D::ChimerGeneAnnotationSR::sort_top_gene_ids($topReplicates,\%fragGeneCharsIndexes,\@fragment2Genes);

	($geneIndex1,$geneIndex2) = compartment_D::ChimerGeneAnnotationSR::select_appropriate_genes(\%fragGeneCharsIndexes,\@fragment1Genes,\@fragment2Genes,\@combChrs);
	
	$chimerTag = compartment_D::ChimerGeneAnnotationSR::determine_chimer_type($geneIndex1,$geneIndex2,\%fragGeneCharsIndexes,\@fragment1Genes,\@fragment2Genes,\@combChrs,$chimerArrayRef,$chimerIndxRef);
	$chimerReadStatsRef->{$chimerTag} += 1;  ## updating stats hash
		
	$foutName = $outFilesRef->{"ChimericTranscriptAnnotation"};

	compartment_D::ChimerGeneAnnotationSR::print_chimer_annotation_to_file($foutName,$chimerTag,$geneIndex1,$geneIndex2,\%fragGeneCharsIndexes,\@fragment1Genes,\@fragment2Genes,\@combChrs,$chimerArrayRef,$chimerIndxRef);

	compartment_D::ChimerGeneAnnotationSR::print_chimer_gene_stats($outFilesRef->{"ExpressedGeneInfo"},$chimerTag,$geneIndex1,$geneIndex2,\%fragGeneCharsIndexes,\@fragment1Genes,\@fragment2Genes,\@combChrs,$chimerArrayRef,$chimerIndxRef);

	shared::DeleteDataStructures::undefine_arrays($chimerArrayRef,\@fragment1Genes,\@fragment2Genes,\@combChrs,\@breakPoints);
	shared::DeleteDataStructures::undefine_hashes(\%fragGeneCharsIndexes);
	}  ## function ends
	###########################
1;
## function to be ended


## ArgLen = 3; (query file, gene track file and output file tag)
## get pslx2 index and gene index
## read in each chimer record
## flip to right if left already
## get involved chromosomes
## get break points
## get all the possible genes involved
	## check for all the possible regions (coding, non-coding, intron or intergenic).
	## sort on coding, non-coding, intron and then intergenic.
	## start with the one in the codign region
## pick one for each fragment
	## match gene names if same chr
	## detect different chimer structures
	## intrachromosomal
		## co-transcription (if genes are adjacent)
		## loop
		## exon re-ordering
		## normal chimer

	## interchromosomal




