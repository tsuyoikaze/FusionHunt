package compartment_C::PrintExonsIntronsOnlyGeneInfo;

require compartment_C::AnnotateExonsOnly;
require compartment_C::AnnotateIntronsOnly;
require compartment_C::PrintGeneAnnotations;

require compartment_C::FindGenicRegions;

use strict;
use warnings;

	sub print_exons_introns_only_gene_info{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$annotTagsRef,$outFilesRef) = @_;
	
	my $annotationTagList;
	my(@geneMappingInfo1,@geneMappingInfo2);  ## will print the gene mapping info that will be used for the gene expression calculation

	$annotationTagList = "";
	$annotationTagList .= $annotTagsRef->[0].',';

		if($annotTagsRef->[0] eq "ExonOnly"){
		my @annotationInfo = compartment_C::AnnotateExonsOnly::get_exons_only_annot_info($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx);
		push(@geneMappingInfo1,($clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"chr"}],
		    $clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],$annotationInfo[scalar(@annotationInfo)-2],"NA","NA","NA","NA"));  ## incrementing number of reads mapped to  this gene + exons of the gene mapped to by sequencing reads
		}
		elsif($annotTagsRef->[0] eq "IntronOnly"){
		my @annotationInfo = compartment_C::AnnotateIntronsOnly::get_introns_only_annot_info($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx);
		push(@geneMappingInfo1,$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"chr"}],
		    $clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],"NA",$annotationInfo[scalar(@annotationInfo)-1],"NA","NA","NA");
		}

	push(@geneMappingInfo2,($annotationTagList,$a1Ref->[$pslx2IndxRef->{"qID"}]));	
	compartment_C::PrintGeneAnnotations::print_full_annot_info(\@geneMappingInfo1,\@geneMappingInfo2,$outFilesRef->{"ExpressedGeneInfo"});

	splice(@geneMappingInfo1);

	splice(@geneMappingInfo2);
	}  ## function ends	
	###########################
	sub check_print_full_transcript_info{
	my($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$cutoff,$outFilesRef,$bedFilesRef) = @_;

	## assuming that the read is already classified as "exon-only" ##
	my(@outArray,@bed);

		## matching boundary conditions for full transcript match  ##
		if(  (($clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}]-$a1Ref->[$pslx2IndxRef->{"tSt"}]) >= 0) && 
			(($clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}]-$a1Ref->[$pslx2IndxRef->{"tSt"}]) <= $cutoff) && 
			(($a1Ref->[$pslx2IndxRef->{"tEn"}] - $clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]) >= 0)  && 
			(($a1Ref->[$pslx2IndxRef->{"tEn"}] - $clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]) <= $cutoff) ){
					
		push(@outArray,($a1Ref->[$pslx2IndxRef->{"qID"}],$a1Ref->[$pslx2IndxRef->{"identity"}],$a1Ref->[$pslx2IndxRef->{"strand"}]));
		push(@outArray,$a1Ref->[$pslx2IndxRef->{"tID"}],($a1Ref->[$pslx2IndxRef->{"tSt"}],$a1Ref->[$pslx2IndxRef->{"tEn"}]));

		push(@outArray,($clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]));

		push(@outArray,($clustRef->[$currTrackIndx][$knGnIndxRef->{"txSt"}],$clustRef->[$currTrackIndx][$knGnIndxRef->{"txEn"}]));
		push(@outArray,(compartment_C::FindGenicRegions::check_if_protein_coding($clustRef,$currTrackIndx,$knGnIndxRef)));
		push(@outArray,$clustRef->[$currTrackIndx][$knGnIndxRef->{"exCount"}]);
		open FHPrintFullTranscript,">>".$outFilesRef->{"FullTranscriptMatch"};
		print FHPrintFullTranscript join("\t",@outArray),"\n";
		close FHPrintFullTranscript;
		
		open FHPrintBED,">>".$bedFilesRef->{"FullTranscriptMatch"};
		@bed = file_parsers::PSLX2ToBED::pslx2_to_bed($pslx2IndxRef,@{$a1Ref});
		print FHPrintBED join("\t",@bed),"\n";  ## send the full array
		close FHPrintBED;
		}
	}  ## function ends
	###########################


1;				



