
package compartment_A::ClassifyBasedOnAlignmentScoresSR;

require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/check_non_uniqueness check_ID_cov_cutoff check_if_fusion/;


use strict;
use warnings;

	sub check_non_uniqueness{  ##9
	my($pslx2IndxRef,$clustRef) = @_;
	my($i,$SimilarHits);

	$SimilarHits = 0;
		if(scalar(@{$clustRef}) > 1){  ## if there are more than one records
			for($i = 1;$i<scalar(@{$clustRef});++$i){
				if( ($clustRef->[$i][$pslx2IndxRef->{"score"}] == $clustRef->[0][$pslx2IndxRef->{"score"}]) 
					&& ($clustRef->[$i][$pslx2IndxRef->{"coverage"}] == $clustRef->[0][$pslx2IndxRef->{"coverage"}]) &&
					($clustRef->[$i][$pslx2IndxRef->{"identity"}] == $clustRef->[0][$pslx2IndxRef->{"identity"}]) ){
				++$SimilarHits;
				}  ## conditions end
			}  ## for($i:1-n) ends
		}  ## if(scalar(@{$clustRef}) > 1) ends
	return(($SimilarHits+1));
	}  ## function ends
	###########################
	sub check_ID_cov_cutoff{  ##10
	my($pslx2IndxRef,$cutoffsRef,$clustRef) = @_;

	my $IdFlag = 0;  ## if identity cutoff is paassing
	my $CovFlag = 0; ## if % coverage cutoff passing
	my $CutoffCheck;  ## 0: OK, 3: low identity,4: low coverage, 5: low id and low cov
	my $LeftCutoff = 0;
	my($ltlv, $rtlv);

		if($clustRef->[0][$pslx2IndxRef->{"identity"}] >= $cutoffsRef->{"IdentityCutoff"}){
		$IdFlag = 1;
		}
		if($clustRef->[0][$pslx2IndxRef->{"coverage"}] >= $cutoffsRef->{"CovCutoff"}){
		$CovFlag = 1;
		}
	
	$ltlv = $clustRef->[0][$pslx2IndxRef->{"qSt"}];  ## counting bases left unaligned on the left end of query
	$rtlv = $clustRef->[0][$pslx2IndxRef->{"qLen"}] - $clustRef->[0][$pslx2IndxRef->{"qEn"}];  ## counting bases left unaligned on the right end of query
	
		if( ($ltlv <= $cutoffsRef->{"LeftLeaveCutoff"}) && ($rtlv <= $cutoffsRef->{"RightLeaveCutoff"}) ){  ## cutoff passed
		$LeftCutoff = 1;
		}

	#### Decide ####
		if( ($IdFlag == 1) && ($CovFlag == 1) && ($LeftCutoff == 1) ){  ## OK to be annotated
		$CutoffCheck = "HighScoring";
		}
		elsif( ($IdFlag == 0) && ($CovFlag == 1) && ($LeftCutoff == 1) ){  ## ID cutoff is not passed
		$CutoffCheck = "Discarded";
		}
		elsif( ($IdFlag == 1) && (($CovFlag == 0) ||  ($LeftCutoff == 0)) ){  ## ID cutoff is OK but Cov is not passed
		$CutoffCheck = "PotentialFusion";
		}
		elsif( ($IdFlag == 0) && (($CovFlag == 0) || ($LeftCutoff == 0)) ){  ## Both ID and Cov cutoffs are not passed
		$CutoffCheck = "Discarded";
		}
	return($CutoffCheck);
	}  ## function ends 
	###########################
	sub check_if_fusion{
	my($pslx2IndxRef,$cutoffsRef,$clustRef) = @_;
	my($LeftCovOK,$RightCovOK,$IfLowCovOK);  ## OK for fusion (0 means no fusion on that side)
	
	$LeftCovOK = 0;
	$RightCovOK = 0;
	
	$IfLowCovOK = "Discarded";  ### default (low cov, high identity but no fusion)
	
		if($clustRef->[0][$pslx2IndxRef->{"qSt"}] >= $cutoffsRef->{"LeftLeaveCutoff"} ){
		$LeftCovOK = 1;  ## left side fusion
		}
		if( ($clustRef->[0][$pslx2IndxRef->{"qLen"}] - $clustRef->[0][$pslx2IndxRef->{"qEn"}]) > $cutoffsRef->{"RightLeaveCutoff"} ){
		$RightCovOK = 1;  ## right side fusion
		}
		
		if(scalar(@{$clustRef}) > 1){
			if( ($LeftCovOK == 1) && ($RightCovOK == 0) ){
			$IfLowCovOK = "PotentialFusionAtLeft";
			}
			elsif( ($LeftCovOK == 0) && ($RightCovOK == 1) ){
			$IfLowCovOK = "PotentialFusionAtRight";
			}
			elsif( ($LeftCovOK == 1) && ($RightCovOK == 1) ){
			$IfLowCovOK = "PotentialFusionAtBoth";
			}
		}  ## if($n>1) ends
	return($IfLowCovOK);
	}  ## function ends
	############################
	sub print_log_file{
	my($out,$fIndx,%LogCounts) = @_;
	my(@outFiles,$i);
	@outFiles = split(/\t+/,$out);  ## outfiles

	open FH18,">>$outFiles[$fIndx]";  ## log file to write
	print FH18 "TotalHits\t",$LogCounts{scalar(@outFiles)},"\n";
		for($i = 0;$i<scalar(@outFiles);++$i){
		next if($i == $fIndx);
		$outFiles[$i] =~ s/(\.*)$//;
		print FH18 "$outFiles[$i]\t$LogCounts{$i}\n";
		}
	
	close FH18;
	}  ## function ends
	################################
1;
