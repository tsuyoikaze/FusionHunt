package compartment_C::FindGenicRegions;

use strict;
use warnings;

	sub check_if_protein_coding{
	my($geneClustRef,$currTrack,$geneIndxRef) = @_;
	my $cdsFlag = "";
		if($geneClustRef->[$currTrack][$geneIndxRef->{"cdsEn"}] - $geneClustRef->[$currTrack][$geneIndxRef->{"cdsSt"}] ){
		$cdsFlag = "coding";
		}
		else{
		$cdsFlag = "non-coding";
		}
	return($cdsFlag);
	}  ## function ends
	###########################
	sub get_cds_pos{
	my($geneClustRef,$currTrack,$geneIndxRef,$cdsStRef,$cdsEnRef) = @_;
	${$cdsStRef} = $geneClustRef->[$currTrack][$geneIndxRef->{"cdsSt"}];
	${$cdsEnRef} = $geneClustRef->[$currTrack][$geneIndxRef->{"cdsEn"}];
	}  ## function ends
	###########################
	sub get_exon_cords{
	my($geneClustRef,$currTrack,$geneIndxRef,$exStRef,$exEnRef) = @_;
	@{$exStRef} = split(',',$geneClustRef->[$currTrack][$geneIndxRef->{"exSt"}]);  ## exon starts
	@{$exEnRef} = split(',',$geneClustRef->[$currTrack][$geneIndxRef->{"exEn"}]);  ## exon ends
	}  ## function ends
	###########################
	sub get_coding_regions{
	my($exStRef,$exEnRef,$cdsSt,$cdsEn,$codStRef,$codEnRef) = @_;
	my $i;
		for($i = 0;$i<scalar(@{$exStRef});++$i){
			if( ($cdsEn-$cdsSt)){ ## if coding region exists
				if( ($exStRef->[$i] >= $cdsSt) && ($exEnRef->[$i] <= $cdsEn) ){
				push(@{$codStRef},$exStRef->[$i]);
				push(@{$codEnRef},$exEnRef->[$i]);
				}
				elsif( ($cdsSt >= $exStRef->[$i]) && ($cdsSt < $exEnRef->[$i]) ){  ## if CDS start lies in the exon
				push(@{$codStRef},$cdsSt);
				push(@{$codEnRef},$exEnRef->[$i]);
				}
				elsif( ($cdsEn > $exStRef->[$i]) && ($cdsEn <= $exEnRef->[$i]) ){  ## if CDS end lies in the exon
				push(@{$codStRef},$exStRef->[$i]);
				push(@{$codEnRef},$cdsEn);
				}
			}  ## if(CDS) ends
		}  ## for (Exons) ends
	}  ## function ends
	###########################
	sub get_noncoding_regions{
	my($exStRef,$exEnRef,$cdsSt,$cdsEn,$nonCodStRef,$nonCodEnRef) = @_;
		if(!($cdsEn-$cdsSt)){  ## if no coding region is present
		@{$nonCodStRef} = @{$exStRef};
		@{$nonCodEnRef} = @{$exEnRef};
		}
	}  ## function ends
	###########################
	sub get_5UTR_regions{
	my($exStRef,$exEnRef,$cdsSt,$cdsEn,$strand,$utr5StRef,$utr5EnRef) = @_;
		if($strand eq "+"){
		get_left_noncod_cords($exStRef,$exEnRef,$cdsSt,$cdsEn,$strand,$utr5StRef,$utr5EnRef);
		}
		else{
		get_right_noncod_cords($exStRef,$exEnRef,$cdsSt,$cdsEn,$strand,$utr5StRef,$utr5EnRef);
		}
	}  ## function ends
	###########################
	sub get_3UTR_regions{
	my($exStRef,$exEnRef,$cdsSt,$cdsEn,$strand,$utr3StRef,$utr3EnRef) = @_;
		if($strand eq "+"){
		get_right_noncod_cords($exStRef,$exEnRef,$cdsSt,$cdsEn,$strand,$utr3StRef,$utr3EnRef);		
		}
		else{
		get_left_noncod_cords($exStRef,$exEnRef,$cdsSt,$cdsEn,$strand,$utr3StRef,$utr3EnRef);		
		}
	}  ## function ends
	###########################
	sub get_left_noncod_cords{
	my($exStRef,$exEnRef,$cdsSt,$cdsEn,$strand,$stRef,$enRef) = @_;
	my $i;
		if( ($cdsEn-$cdsSt)){  ## if coding region exists
			for($i = 0;$i<scalar(@{$exStRef});++$i){
				if($exEnRef->[$i] < $cdsSt){  ## if the exon doesn't contain cds start
				push(@{$stRef},$exStRef->[$i]);
				push(@{$enRef},$exEnRef->[$i]);
				}
				elsif( ($cdsSt > $exStRef->[$i]) && ($cdsSt < $exEnRef->[$i]) ){  ## if CDS start lies in the exon
				push(@{$stRef},$exStRef->[$i]);
				push(@{$enRef},$cdsSt);
				}
			}  ## for ($i) ends
		}  ## if(cds region exists) ends
	}  ## function ends
	###########################
	sub get_right_noncod_cords{
	my($exStRef,$exEnRef,$cdsSt,$cdsEn,$strand,$stRef,$enRef) = @_;
	my $i;
		if( ($cdsEn-$cdsSt)){  ## if coding region exists
			for($i = 0;$i<scalar(@{$exStRef});++$i){
				if($cdsEn < $exStRef->[$i]){  ## if the exon doesn't contain cds end
				push(@{$stRef},$exStRef->[$i]);
				push(@{$enRef},$exEnRef->[$i]);
				}
				elsif( ($cdsEn > $exStRef->[$i]) && ($cdsEn < $exEnRef->[$i]) ){  ## if exon contains cds end
				push(@{$stRef},$cdsEn);
				push(@{$enRef},$exEnRef->[$i]);
				}
			}  ## for($i) ends
		}  ## if(cds region exists) ends
	}  ## function ends
	###########################
	sub get_intron_regions{
	my($exStRef,$exEnRef,$inStRef,$inEnRef) = @_;
	my $i;
		for($i = 0;$i<(scalar(@${exStRef})-1);++$i){
		push(@{$inStRef},$exEnRef->[$i]);
		push(@{$inEnRef},$exStRef->[$i+1]);
		}  ## for(scalar(@{$exStRef})-1) ends
	}  ## function ends
	###########################
	sub which_intron{
	my($currIntron,$totalExons,$strand) = @_;
	my $whichIntron;		
			if($strand eq '+'){
			$whichIntron = ($currIntron+1);
			}
			else{
			$whichIntron = (($totalExons-1)-$currIntron); ## which intron got extension
			}
	return($whichIntron);
	}  ## function ends
	###########################
	sub which_exon{
	my($currExon,$totalExons,$strand) = @_;
	my $whichExon;		
			if($strand eq '+'){
			$whichExon = ($currExon+1);
			}
			else{
			$whichExon = ($totalExons-$currExon); ## which intron got extension
			}
	return($whichExon);
	}  ## function ends
	###########################
	sub get_transcript_length{
	my($geneClustRef,$currTrack,$geneIndxRef) =  @_;
	my(@exSt,@exEn,$transcriptLen,$i);
	get_exon_cords($geneClustRef,$currTrack,$geneIndxRef,\@exSt,\@exEn);  ## getting exon coordinates;
	$transcriptLen = 0;
		for($i = 0;$i<scalar(@exSt);++$i){
		$transcriptLen += $exEn[$i]-$exSt[$i];
		}
	splice(@exSt);
	splice(@exEn);
	return($transcriptLen);
	}  ## function ends
	###########################
1;