######### Perl script to read read clusters from the file usinf file pointer ##########

package Control::RunCharacterizationModules;

use strict;
use warnings;

require Begin::CreateOutputFiles;

require shared::ArrayOperations;
require shared::CoordinateComparisons;
require shared::DeleteDataStructures;
require shared::FileCheckPrint;

require compartment_A::CalculateScIdCov;
require compartment_A::SortToFindTopHit;
require compartment_A::AlignmentBasedClassification;
require compartment_C::AnnotatePSLX2Reads;
require compartment_B::DetectChimera;
require compartment_D::AnnotateChimers;
require compartment_E::CalculateGeneExpression;
require Stats::StatsCollectionHashes;


##obtain: $inFile, $fpSt, $fpEn, $outfileTag, cutoffs, annotation stringency, pslIndex, pslx2Index,knGnIndex, gene file name, gene file index, 


	sub characterize_alignment_clusters{
	my($inFile,$alignPslFlag,$fpSt,$fpEn,$geneFileRef,$scriptParametersRef,$pslxIndxRef,$pslx2IndxRef,$knGnIndxRef,
		$chimerIndxRef,$geneFileIndxRef,$outFilePrefix) = @_;

	my($str,@a1,@cluster,$n);
	my($pt1,$pt2);
	my($alignmentAnnotTag);
	my($rightChimerFlag,@rightChimerArray,$leftChimerFlag,@leftChimerArray);
	my $chimerDiscardFlag = 0;  ## low coverage but not chimeric transcript
	my $totalReads = 0;  ## total reads in this thread

	
	#### required output files ####
	my %characterizationOutFiles = Begin::CreateOutputFiles::characterization_out_file($outFilePrefix);

	my %bedOutFiles = Begin::CreateOutputFiles::BED_out_file($outFilePrefix);  ## only hash was returned (files have not been created yet)

	Begin::CreateOutputFiles::create_file_from_hash(\%bedOutFiles);  ## creating files


	#### statistics hashes ####
	my %alignmentClassStats = Stats::StatsCollectionHashes::alignment_class_hash();  ## hash to store counts of score baesd alignments classes
	my %readAnnotationStats = Stats::StatsCollectionHashes::annotation_class_hash();  ## hash to store counts of each type of read annotation

	my %geneCharStats = Stats::StatsCollectionHashes::annotation_class_hash();  ## hash to store number of genes annotated to different characterization classes
	my %geneDistHash = Stats::StatsCollectionHashes::gene_mapping_distribution_hash();  ## hash to store number of genes expressed in the given sample
	my %chimerReadStats = Stats::StatsCollectionHashes::chimer_distribution_hash();  ## to store chimeric read counts
	my %chimerGeneStats = Stats::StatsCollectionHashes::chimer_distribution_hash();  ## to store chimeric gene counts
	my %exonsOnlyStats = Stats::StatsCollectionHashes::exons_only_read_hash();  ## to store read mapping regions stats for exons only characteristics
	
	%geneCharStats = (%geneCharStats,%chimerGeneStats);  ## completing gene charasteristics hash (splicing and chimers)

	###########################
	
	open FHCharacterizeIn,$inFile or die "\n can not open file $inFile \n";

	seek(FHCharacterizeIn,$fpSt,0);  ## setting up the starting pointer

	$n = 0;
	$pt1 = tell FHCharacterizeIn;  ## recordin FPs starting location

	##print "\nptSt: $fpSt, ptEn: $fpEn\n";

		while($str = <FHCharacterizeIn>){
		$pt2 = tell FHCharacterizeIn;  ## recording FPs second location
		$str =~ s/\n//;
		$str =~ s/\r//;
		@a1 = split(/\t+/,$str);

			if(scalar(@a1) == 21){
			push(@a1,('.','.'));  ## converting psl to pslx
			}

			if( !($n) || ($cluster[scalar(@cluster)-1][$pslxIndxRef->{"qID"}] eq $a1[$pslxIndxRef->{"qID"}]) ){
			push(@cluster,[@a1]);
			++$n;
			$pt1 = $pt2;  ## incrementing the file pointer by one line
			}

			if( (eof) || ($cluster[$n-1][$pslxIndxRef->{"qID"}] ne $a1[$pslxIndxRef->{"qID"}]) || ((tell FHCharacterizeIn) >= $fpEn) ){
			++$totalReads;	
	
			## PSLX TO PSLX2 ##
			@cluster = compartment_A::CalculateScIdCov::main_module_CalculateScIdCov($pslxIndxRef,\@cluster);
	
			## SORT ALIGNMENTS ##
			compartment_A::SortToFindTopHit::main_module_SortToFindTopHit($pslx2IndxRef,\@cluster);

	
			## DECISION MAKING ##
			$alignmentAnnotTag = compartment_A::AlignmentBasedClassification::alignment_based_classification($scriptParametersRef,$pslx2IndxRef,\@cluster,\%characterizationOutFiles,\%alignmentClassStats);
	
	
			## ANNOTATION OF HIGH SCORING READS ##
		
				if($alignmentAnnotTag eq "HighScoring"){  ## high scoring reads
				compartment_C::AnnotatePSLX2Reads::annotate_pslx2_records($scriptParametersRef,$pslx2IndxRef,
				\@cluster,$knGnIndxRef,$geneFileIndxRef,$geneFileRef,\%characterizationOutFiles,\%bedOutFiles,\%readAnnotationStats,\%exonsOnlyStats);
				++$alignmentClassStats{"HighScoring"};
				}  ## if(alignmentAnnotTag eq HighScoring) ends
	


			## DETECT CHIMERIC TRANSCRIPTS
	
				if($alignmentAnnotTag eq "PotentialFusionAtRight"){
				($rightChimerFlag,@rightChimerArray) = compartment_B::DetectChimera::detect_chimeric_transcripts("Right",$scriptParametersRef,$pslx2IndxRef,\@cluster,\%characterizationOutFiles);


					if($rightChimerFlag){
					compartment_D::AnnotateChimers::annotate_chimeric_transcripts($chimerIndxRef,\@rightChimerArray,$knGnIndxRef,$geneFileIndxRef,$geneFileRef,\%characterizationOutFiles,\%chimerReadStats);
					++$alignmentClassStats{"Chimers"};
					}
					else{
					$chimerDiscardFlag = 1;
					}
				}
				elsif($alignmentAnnotTag eq "PotentialFusionAtLeft"){
				($leftChimerFlag,@leftChimerArray) = compartment_B::DetectChimera::detect_chimeric_transcripts("Left",$scriptParametersRef,$pslx2IndxRef,\@cluster,\%characterizationOutFiles);
	

					if($leftChimerFlag){
					compartment_D::AnnotateChimers::annotate_chimeric_transcripts($chimerIndxRef,\@leftChimerArray,$knGnIndxRef,$geneFileIndxRef,$geneFileRef,\%characterizationOutFiles,\%chimerReadStats);
					++$alignmentClassStats{"Chimers"};
					}
					else{
					$chimerDiscardFlag = 1;
					}
				}
				elsif($alignmentAnnotTag eq "PotentialFusionAtBoth"){
				($rightChimerFlag,@rightChimerArray) = compartment_B::DetectChimera::detect_chimeric_transcripts("Right",$scriptParametersRef,$pslx2IndxRef,\@cluster,\%characterizationOutFiles);

				($leftChimerFlag,@leftChimerArray) = compartment_B::DetectChimera::detect_chimeric_transcripts("Left",$scriptParametersRef,$pslx2IndxRef,\@cluster,\%characterizationOutFiles);


					if( (($rightChimerFlag) && ($leftChimerFlag)) || (!($rightChimerFlag) && !($leftChimerFlag)) ){
					$chimerDiscardFlag = 1;
					}
					elsif( ($rightChimerFlag) && !($leftChimerFlag) ){
					compartment_D::AnnotateChimers::annotate_chimeric_transcripts($chimerIndxRef,\@rightChimerArray,$knGnIndxRef,$geneFileIndxRef,$geneFileRef,\%characterizationOutFiles,\%chimerReadStats);
					++$alignmentClassStats{"Chimers"};
					}
					elsif( !($rightChimerFlag) && ($leftChimerFlag) ){
					compartment_D::AnnotateChimers::annotate_chimeric_transcripts($chimerIndxRef,\@leftChimerArray,$knGnIndxRef,$geneFileIndxRef,$geneFileRef,\%characterizationOutFiles,\%chimerReadStats);
					++$alignmentClassStats{"Chimers"};
					}
				}  ## elsif(PotentialFusionAtBoth) ends
	
				if($chimerDiscardFlag){  ## if not chimer (then discard)
				shared::FileCheckPrint::print_cluster_to_file($characterizationOutFiles{"Discarded"},scalar(keys(%{$pslx2IndxRef})),1,\@cluster);
				++$alignmentClassStats{"Discarded"};
				}


	

		
			## DESTRYOING PREVIOUS CLUSTER TO START OVER ##
			$n = 0;
			seek(FHCharacterizeIn,$pt1,0);  ## resetting the file pointer
			$chimerDiscardFlag = 0;
			shared::DeleteDataStructures::undefine_arrays(\@cluster,\@rightChimerArray,\@leftChimerArray);

				if( (tell FHCharacterizeIn) >= $fpEn){  ## end of input data chunk
				last;
				}

			}  ## if($id1 ne $id2) ends
		}  ## while(<FHFHCharacterizeIn>) ends
	close FHCharacterizeIn;

	## printing stats ##
	shared::FileCheckPrint::print_hash_to_file(\%readAnnotationStats,$characterizationOutFiles{"GeneAnnotationStats"});
	shared::FileCheckPrint::print_hash_to_file(\%alignmentClassStats,$characterizationOutFiles{"AlignmentStats"});
	shared::FileCheckPrint::print_hash_to_file(\%chimerReadStats,$characterizationOutFiles{"GeneAnnotationStats"});
	shared::FileCheckPrint::print_hash_to_file(\%exonsOnlyStats,$characterizationOutFiles{"GeneAnnotationStats"});
		####################
	##print "\ntotalReads: $totalReads\n";
	return($totalReads);
	}  ## function ends
	###########################
1;	
