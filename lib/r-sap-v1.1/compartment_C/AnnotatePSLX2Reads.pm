############# perl scritpt to annotate the sequencing reads to gene tracks from UCSC. Gene tracks should be from the genome ###############
############# assembly same as the one used for the initial alignment of the sequencing reads to the reference genome. ####################

package compartment_C::AnnotatePSLX2Reads;

use strict;
use warnings;

require compartment_C::DeterminePrintAnnotations;
require compartment_C::AnnotationProcessModules;
require compartment_C::DetermineCharacteristicsModules;

require file_parsers::PSLX2ToBED;


	sub annotate_pslx2_records{
	my($cutoffHashRef,$pslx2IndxRef,$clustRef,$knGnIndxRef,$gnFileIndxRef,$geneFileRef,$outFilesRef,$bedFilesRef,$readAnnotStatsRef,$exonsOnlyHashRef) = @_;
	
	my($str1,$str2,@a1,@tempClust);
	my(@knGn,$indx);
	my(@cluster,$n,$i);
	my(@bed);
	$n = 0;	

	@a1 = @{$clustRef->[0]}[0..(scalar(keys(%{$pslx2IndxRef}))-1)];  ## tophit from the cluster

		if(compartment_C::DetermineCharacteristicsModules::check_if_large_intron(\@a1,$pslx2IndxRef,$cutoffHashRef->{"IntronCutoff"})){  ## checking for large introns inserts
		push(@tempClust,[@a1]);
		shared::FileCheckPrint::print_cluster_to_file($outFilesRef->{"Uncharacterized"},scalar(keys(%{$pslx2IndxRef})),1,\@tempClust);
		splice(@tempClust);
		++$readAnnotStatsRef->{"Uncharacterized"};

		open FHPrintBED,">>".$bedFilesRef->{"Uncharacterized"}; #### printing to BED file parts ####
		@bed = file_parsers::PSLX2ToBED::pslx2_to_bed($pslx2IndxRef,@a1);
		print FHPrintBED join("\t",@bed),"\n";
		close FHPrintBED;
		}
		
		else{
			if(exists $gnFileIndxRef->{$a1[$pslx2IndxRef->{"tID"}]}){  ## chromosome name

			$indx = $gnFileIndxRef->{$a1[$pslx2IndxRef->{"tID"}]};  ## getting gene track index using chromosome name
	
			##open FH21,$ARGV[1] or die "\n can not open file $ARGV[0] \n";
			##seek(FH21,$indx,0);  ## setting pointer to the start of the chromosome in the gene track file

			##open FHAnnotatePSLX21,$geneFile or die "\n can not open file $geneFile \n";
			##seek(FHAnnotatePSLX21,$indx,0);  ## setting pointer to the start of the chromosome in the gene track file

				##while($str2 = <FHAnnotatePSLX21>){
				##$str2 =~ s/\n//;
				##$str2 =~ s/\r//;
				##@knGn = split(/\t+/,$str2);

				for($i = $indx;$i<scalar(@{$geneFileRef});++$i){
				@knGn = split(/\t+/,$geneFileRef->[$i]);

					if($a1[$pslx2IndxRef->{"tID"}] ne $knGn[$knGnIndxRef->{"chr"}]){  ## if chr doesn't match
					last;
					}
				
						
					if(shared::CoordinateComparisons::intersect_track_alignment($a1[$pslx2IndxRef->{"tSt"}],$a1[$pslx2IndxRef->{"tEn"}],
					($knGn[$knGnIndxRef->{"txSt"}]-$cutoffHashRef->{"GeneRadius"}),
					($knGn[$knGnIndxRef->{"txEn"}]+$cutoffHashRef->{"GeneRadius"}))){  ## if there is any intersect

						if($cutoffHashRef->{"stMatch"} =~ m/(t)/i){  ## if strand match required
							if($a1[$pslx2IndxRef->{"strand"}] eq $knGn[$knGnIndxRef->{"strand"}]){  ## if strand matches	
							push(@cluster,[@knGn]);  ## collecting gene tracks for read annotation
							++$n;
							}  ## if(strand matches ends
						}  ## if(strand match required) ends
						else{
						push(@cluster,[@knGn]);  ## collecting gene tracks for read annotation
						++$n;
						}
					}  ## if(intersecting) ends
				splice(@knGn);
				}  ## for(@geneFile) ends

				##}  ## while(<FHFHAnnotatePSLX21>) ends

			##close FHAnnotatePSLX21;  ## closing gene track file

			}  ## if(chr index exists) ends

			if($n){  ## if at least one gene intersect was found
			compartment_C::DeterminePrintAnnotations::determine_print_annotations($pslx2IndxRef,$knGnIndxRef,\@a1,\@cluster,$cutoffHashRef,$outFilesRef,$bedFilesRef,$readAnnotStatsRef,$exonsOnlyHashRef);  ## function call 
			}  ## if($n) ends
	
			else{  ## if no gene intersect was found
			open FHAnnotatePSLX22,">>".$outFilesRef->{"GeneDesert"};  ## gene deserts
			print FHAnnotatePSLX22 join("\t",@a1),"\n";
			close FHAnnotatePSLX22;
			++$readAnnotStatsRef->{"GeneDesert"};

			open FHPrintBED,">>".$bedFilesRef->{"GeneDesert"}; #### printing to BED file parts ####
			@bed = file_parsers::PSLX2ToBED::pslx2_to_bed($pslx2IndxRef,@a1);
			print FHPrintBED join("\t",@bed),"\n";
			close FHPrintBED;
			}
		}  ## elsif(!(uncharacterized)) ends
	}  ## function ends
	###########################
1;
