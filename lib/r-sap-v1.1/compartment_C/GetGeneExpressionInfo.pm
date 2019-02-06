package compartment_C::GetGeneExpressionInfo;

use strict;
use warnings;

require compartment_C::AnnotateExonsOnly;
require compartment_C::FindGenicRegions;

	sub get_gene_exp_info{
	my ($pslx2IndxRef,$knGnIndxRef,$a1Ref,$clustRef,$currTrackIndx,$outFileName) = @_;
	my(@annotationInfo,@outArray);
	my($cdsSt,$cdsEn);  ## ORF start, ORF end
	@annotationInfo = compartment_C::AnnotateExonsOnly::get_exons_only_annot_info($pslx2IndxRef,$knGnIndxRef,
		$a1Ref,$clustRef,$currTrackIndx);
	
	push(@outArray,$clustRef->[$currTrackIndx][$knGnIndxRef->{"transcriptID"}]);
	push(@outArray,(compartment_C::FindGenicRegions::get_transcript_length($clustRef,
		$currTrackIndx,$knGnIndxRef)));  ## lengt of the mature transcript (sum of exons)

	push(@outArray,$clustRef->[$currTrackIndx][$knGnIndxRef->{"chr"}]);  ## chromosome
	push(@outArray,$clustRef->[$currTrackIndx][$knGnIndxRef->{"strand"}]);  ## strand
	push(@outArray,$annotationInfo[1]);  ## included exons
	compartment_C::FindGenicRegions::get_cds_pos($clustRef,$currTrackIndx,$knGnIndxRef,\$cdsSt,\$cdsEn);
		if($cdsEn-$cdsSt){
		push(@outArray,"protein_coding");
		}
		else{
		push(@outArray,"non_protein_coding");
		}
	open FHGeneExpInfo,">>$outFileName";  ## out file
	print FHGeneExpInfo join("\t",@outArray),"\n";
	close FHGeneExpInfo;
	splice(@annotationInfo);
	}  ## function ends
	###########################

