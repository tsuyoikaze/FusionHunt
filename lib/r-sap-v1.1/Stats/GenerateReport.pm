package Stats::GenerateReport;

use shared::CheckStringsNumericals;
use shared::DeleteDataStructures;


use strict;
use warnings;

	sub generate_report{
	my $outFileName = shift;	

	my $totalReads = shift;
	my $alignmentStatsHashRef = shift;
	my $readAnnotationStatsRef = shift;
	my $geneCharStatsRef = shift;
	my $chimerReadStatsRef = shift;
	my $geneDistHashRef = shift;
	my $exonsOnlyHashRef = shift;
	my $exonSkippingHashRef = shift;
	my $retainedIntronHashRef = shift;
	my $splicedGenicRegionStats = shift;
	

	
	open FHPrintReport,">>$outFileName";  ## output file
	print FHPrintReport "Total number of mapped reads: $totalReads\n\n";

	my($mainKey,$subKey);  ## hash key to access hash references
	my($tempStr,$tempPer);  ## for storing number percentage string

	#### string lengths for formating ####
	my $maxCharLen = 22;
	my $maxNumPerLen = 21;
	my $maxNumLen = 7;
	my $longTableLen = 72;
	my $shortTableLen = 40;

	#### temporary arrays for storing lengths ####
	my(@tabL1,@tabL2,@tabL3);  ## table lines: = - _ respectively
	my(@charArray,@numPerArray,@numArray);  ## value arrays
	my(@lineArray);  ## to store the arrays above and print to file at onece
	my(@temp1,@temp2);  ## temporary arrays (to store variable epmpty values)
	

	#### Initializing arrays ####
	@tabL1 = shared::ArrayOperations::initialize_value_array("=",$longTableLen);
	@tabL2 = shared::ArrayOperations::initialize_value_array("-",$longTableLen);
	@tabL3 = shared::ArrayOperations::initialize_value_array("_",$longTableLen);
	


	## Table 1 header ##
	print FHPrintReport "Table1: Distribution of mappable reads.\n";
	print FHPrintReport join('',@tabL1),"\n";
	push(@temp1,split('',"Description"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t\t"));
	splice(@temp1);  splice(@temp2);

	push(@temp1,split('',"Reads(%)"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxNumPerLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"Reference transcripts"));
	print FHPrintReport join('',@lineArray),"\n",join('',@tabL2),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	#### header printing ends ####

	###########################
	## printing high-scoring ##
	###########################
	$mainKey = "HighScoring";

	push(@temp1,split('',$mainKey));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t\t\t"));
	splice(@temp1);  splice(@temp2);
	
	$tempPer = shared::CheckStringsNumericals::get_percentage($alignmentStatsHashRef->{$mainKey},$totalReads);  ## also checking for illigal division by zero
	$tempStr = $alignmentStatsHashRef->{$mainKey}." (".sprintf("%.2f",$tempPer)."%)";
	
	push(@temp1,split('',$tempStr));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxNumPerLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2));
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	
		## printing annotation stats ##
		for $subKey(keys(%{$readAnnotationStatsRef})){
		push(@temp1,split('',$subKey));
		@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
		push(@lineArray,("\t",@temp1,@temp2));
		splice(@temp1);  splice(@temp2);
		
		$tempPer = shared::CheckStringsNumericals::get_percentage($readAnnotationStatsRef->{$subKey},$alignmentStatsHashRef->{$mainKey});
		$tempStr = $readAnnotationStatsRef->{$subKey}." (".sprintf("%.2f",$tempPer)."%)";
		push(@temp1,split('',$tempStr));
		@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxNumPerLen-scalar(@temp1)));
		push(@lineArray,("\t\t",@temp1,@temp2));
		splice(@temp1);  splice(@temp2);


		$tempStr = $geneCharStatsRef->{$subKey};
		push(@temp1,split('',$tempStr));
		push(@lineArray,("\t",@temp1));
		print FHPrintReport join('',@lineArray),"\n";
		splice(@temp1); splice(@temp2); splice(@lineArray);
		}  ## for high-scoring annotation hash ends
	print FHPrintReport "\n";

	###################################
	## printing chimeric transcripts ##
	###################################
	$mainKey = "Chimers";
	push(@temp1,split('',$mainKey));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t\t\t"));
	splice(@temp1);  splice(@temp2);
	
	$tempPer = shared::CheckStringsNumericals::get_percentage($alignmentStatsHashRef->{$mainKey},$totalReads);  ## also checking for illigal division by zero
	$tempStr = $alignmentStatsHashRef->{$mainKey}." (".sprintf("%.2f",$tempPer)."%)";
	
	push(@temp1,split('',$tempStr));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxNumPerLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2));
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);

		## printing chimeric transcript annotation stats ##
		for $subKey(keys(%{$chimerReadStatsRef})){
		push(@temp1,split('',$subKey));
		@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
		push(@lineArray,("\t",@temp1,@temp2));
		splice(@temp1);  splice(@temp2);
		

		$tempPer = shared::CheckStringsNumericals::get_percentage($chimerReadStatsRef->{$subKey},$alignmentStatsHashRef->{$mainKey});
		$tempStr = $chimerReadStatsRef->{$subKey}." (".sprintf("%.2f",$tempPer)."%)";

		push(@temp1,split('',$tempStr));
		@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxNumPerLen-scalar(@temp1)));
		push(@lineArray,("\t\t",@temp1,@temp2));
		splice(@temp1);  splice(@temp2);


		$tempStr = $geneCharStatsRef->{$subKey};
		push(@temp1,split('',$tempStr));
		push(@lineArray,("\t",@temp1));
		print FHPrintReport join('',@lineArray),"\n";
		splice(@temp1); splice(@temp2); splice(@lineArray);
		}  ## for chimeric transcript annotation hash ends
	print FHPrintReport "\n";

	###############################
	## printing multi-hits reads ##
	###############################
	$mainKey = "MultiHits";
	push(@temp1,split('',$mainKey));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t\t\t"));
	splice(@temp1);  splice(@temp2);
	
	$tempPer = shared::CheckStringsNumericals::get_percentage($alignmentStatsHashRef->{$mainKey},$totalReads);  ## also checking for illigal division by zero
	$tempStr = $alignmentStatsHashRef->{$mainKey}." (".sprintf("%.2f",$tempPer)."%)";
	
	push(@temp1,split('',$tempStr));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxNumPerLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2));
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);

	##############################
	## printing Discarded reads ##
	##############################
	$mainKey = "Discarded";
	push(@temp1,split('',$mainKey));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t\t\t"));
	splice(@temp1);  splice(@temp2);
	
	$tempPer = shared::CheckStringsNumericals::get_percentage($alignmentStatsHashRef->{$mainKey},$totalReads);  ## also checking for illigal division by zero
	$tempStr = $alignmentStatsHashRef->{$mainKey}." (".sprintf("%.2f",$tempPer)."%)";
	
	push(@temp1,split('',$tempStr));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxNumPerLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2));
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	
	print FHPrintReport join('',@tabL3),"\n\n\n\n\n";
	#### Table 1 printing ends ####



	######################
	## printing Table 2 ##
	######################

	## Table 2 header ##
	@tabL1 = shared::ArrayOperations::initialize_value_array("=",$shortTableLen-$maxNumLen);
	@tabL2 = shared::ArrayOperations::initialize_value_array("-",$shortTableLen-$maxNumLen);
	@tabL3 = shared::ArrayOperations::initialize_value_array("_",$shortTableLen-$maxNumLen);


	print FHPrintReport "Table2: Distribution of Exon-only\nreads over different genic regions.\n",join('',@tabL1),"\n";
	push(@temp1,split('',"Description"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t"));
	splice(@temp1);  splice(@temp2);

	push(@lineArray,("#Reads"));
	print FHPrintReport join('',@lineArray),"\n",join('',@tabL2),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	#### header printing ends ####	


		for $subKey(sort keys(%{$exonsOnlyHashRef})){
		push(@temp1,split('',$subKey));
		@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
		push(@lineArray,(@temp1,@temp2,"\t"));
		splice(@temp1);  splice(@temp2);

		$tempStr = $exonsOnlyHashRef->{$subKey};
		push(@lineArray,$tempStr);
		print FHPrintReport join('',@lineArray),"\n";
		splice(@lineArray);
		}  ## for (exon-only) distribution ends
	
	print FHPrintReport join('',@tabL3),"\n\n\n\n\n";
	#### Table 2 printing ends ####
	
	
	
	######################
	## printing Table 3 ##
	######################	
	@tabL1 = shared::ArrayOperations::initialize_value_array("=",$shortTableLen+$maxNumLen);
	@tabL2 = shared::ArrayOperations::initialize_value_array("-",$shortTableLen+$maxNumLen);
	@tabL3 = shared::ArrayOperations::initialize_value_array("_",$shortTableLen+$maxNumLen);

	print FHPrintReport "Table3: Distribution of detected reference transcripts.\n";
	
	print FHPrintReport join('',@tabL1),"\n";
	push(@temp1,split('',"Description"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t"));
	splice(@temp1);  splice(@temp2);
	push(@lineArray,"#Reference transcripts");
	print FHPrintReport join('',@lineArray),"\n";
	print FHPrintReport join('',@tabL2),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	#### header printing ends ###

	
	push(@temp1,split('',"Total detected"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t\t"));
	push(@lineArray,$geneDistHashRef->{"detectedGenes"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);

	push(@temp1,split('',"protein-coding"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,("\t",@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$geneDistHashRef->{"codingDetected"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);

	push(@temp1,split('',"non-protein-coding"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,("\t",@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$geneDistHashRef->{"non-codingDetected"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	@temp2 = shared::ArrayOperations::initialize_value_array("_",($maxCharLen+8));
	print FHPrintReport "\t",join('',@temp2),"\n\n";
	

	push(@temp1,split('',"normally-spliced"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,("\t",@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$geneDistHashRef->{"normallySpliced"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	
	push(@temp1,split('',"Single novel event"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,("\t",@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$geneDistHashRef->{"oneAbSplicing"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);

	push(@temp1,split('',"multiple-novel-event"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,("\t",@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$geneDistHashRef->{"multiAbSplicing"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);

	print FHPrintReport join('',@tabL3),"\n\n\n\n\n";
	#### printing Table 3 ends ####

	

	######################
	## printing Table 4 ##
	######################
	@tabL1 = shared::ArrayOperations::initialize_value_array("=",$shortTableLen);
	@tabL2 = shared::ArrayOperations::initialize_value_array("-",$shortTableLen);
	@tabL3 = shared::ArrayOperations::initialize_value_array("_",$shortTableLen);

	print FHPrintReport "Table 4: Exon-skipping reads, reference\ntrasncripts and exons.\n";
	print FHPrintReport join('',@tabL1),"\n";


	push(@temp1,split('',"Exon-skipping reads"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen+$maxNumLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$exonSkippingHashRef->{"skippedExonReads"});

	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	
	push(@temp1,split('',"Reference transcripts"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen+$maxNumLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$exonSkippingHashRef->{"skipExonGenes"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	
	
	push(@temp1,split('',"Skipped exons"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen+$maxNumLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$exonSkippingHashRef->{"skippedExons"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	
	print FHPrintReport join('',@tabL3),"\n\n\n\n\n";
	#### printing Table 4 ends #####
	
	
	######################	
	## printing Table 5 ##
	######################
	print FHPrintReport "Table 5: Intron-retention reads, reference\ntranscripts and introns.\n";
	print FHPrintReport join('',@tabL1),"\n";
	
	push(@temp1,split('',"Intron-retention reads"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen+$maxNumLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$retainedIntronHashRef->{"retainedIntronReads"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);

	push(@temp1,split('',"Reference transcripts"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen+$maxNumLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$retainedIntronHashRef->{"retainIntronGenes"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	
	
	push(@temp1,split('',"Retained-introns"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen+$maxNumLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\t"));
	splice(@temp1); splice(@temp2);
	push(@lineArray,$retainedIntronHashRef->{"retainedIntrons"});
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	
	print FHPrintReport join('',@tabL3),"\n\n\n\n\n";
	#### printing Table 5 ends ###
	
	
	######################
	## printing Table 6 ##
	######################
	print FHPrintReport "Table 6: Distribution of reference transcripts\nwith novel", 
		"transcriptional events in different\ntranscript regions.\n";
	print FHPrintReport join('',@tabL1),"\n";
	
	push(@temp1,split('',"Description"));
	@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
	push(@lineArray,(@temp1,@temp2,"\tReference transcripts"));
	print FHPrintReport join('',@lineArray),"\n";
	splice(@temp1); splice(@temp2); splice(@lineArray);
	print FHPrintReport join('',@tabL2),"\n";
	
		for $subKey (sort(keys(%{$splicedGenicRegionStats}))){
			if($subKey eq "NonCoding"){
			push(@temp1,split('',"Non-protein-coding"));
			}
			elsif($subKey eq "MultipleRegions"){
			push(@temp1,split('',"Multiple-regions"));
			}
			else{
			push(@temp1,split('',$subKey));
			}
		@temp2 = shared::ArrayOperations::initialize_value_array(" ",($maxCharLen-scalar(@temp1)));
		push(@lineArray,(@temp1,@temp2,"\t"));
		splice(@temp1); splice(@temp2);
		push(@lineArray,$splicedGenicRegionStats->{$subKey});
		print FHPrintReport join('',@lineArray),"\n";
		splice(@temp1); splice(@temp2); splice(@lineArray);
		}  ## for (spliced region hash) ends
	print FHPrintReport join('',@tabL3),"\n\n\n";
	### printing Table 6 ends ####
	}  ## function ends
	
1;
