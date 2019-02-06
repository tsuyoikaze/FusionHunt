package compartment_D::ChimerGeneAnnotationSR;


use strict;
use warnings;

use compartment_D::ChimerGeneFindSR;

	sub flip_fusion_record{
	my($chimerIndxRef,$a1Ref) = @_;
	my($fillSide,@temp1,@temp2);

		$a1Ref->[$chimerIndxRef->{"qCov1"}] =~ s/\(//;
		$a1Ref->[$chimerIndxRef->{"qCov2"}] =~ s/\)//;

		##print "\n",$a1Ref->[$chimerIndxRef->{"qCov1"}],"\n";
		##print $a1Ref->[$chimerIndxRef->{"qCov2"}],"\n\n";

		@temp1 = split('-',$a1Ref->[$chimerIndxRef->{"qCov1"}]);  ## query coordinates part1
		@temp2 = split('-',$a1Ref->[$chimerIndxRef->{"qCov2"}]);  ## query coordinates part2
		
			if($temp1[0] > $temp2[0]){
			$fillSide = "left";
			}
			elsif($temp1[0] < $temp2[0]){
			$fillSide = "right";
			}

			if($fillSide eq "left"){
			flip_side_fusion_record($chimerIndxRef,$a1Ref);  ## changing record format
			}

	}  ## funciton ends
	###########################
	sub flip_side_fusion_record{
	my($chimerIndxRef,$a1Ref) = @_;	
	my(@a2,@temp,$str,$i);	
	push(@a2,@{$a1Ref}[0..1]);  ## pushing first 2 elements	
	find_connected_chrs($chimerIndxRef,\@temp,$a1Ref);  ## getting fused chromosomes

	$str = "(".$temp[1]."-".$temp[0].")";	
	push(@a2,$str);  ## pushing chromosomes	
	splice(@temp);	
		
	$a1Ref->[3] =~ s/\(//;	
	$a1Ref->[3] =~ s/\)//;	
	@temp = split('-',$a1Ref->[3]);	
	$str = $str = "(".$temp[1]."-".$temp[0].")";	
	push(@a2,$str);  ## pushing identity	
	splice(@temp);	
		
	$a1Ref->[4] =~ s/\(//;	
	$a1Ref->[4] =~ s/\)//;	
	@temp = split(',',$a1Ref->[4]);	
	$str = "(".$temp[1].",".$temp[0].")";	
	push(@a2,$str);  ## pushing indexes	
	splice(@temp);	
		
	push(@a2,@{$a1Ref}[5..7]);  ## pushing first 2 elements	
	
	$a1Ref->[9] =~ s/\)//;	
	$str = "(".$a1Ref->[9];	
	push(@a2,$str);  ## pushing left part coordinates of query	
	
	$a1Ref->[8] =~ s/\(//;	
	$str = $a1Ref->[8].")";	
	push(@a2,$str);  ## pushing right part coordinates of query	
	
	push(@a2,@{$a1Ref}[20..29]);  ## pushing left part details	
	push(@a2,@{$a1Ref}[10..19]);  ## pushing right part details	

		for($i = 0;$i<30;++$i){
		$a1Ref->[$i] = $a2[$i];
		}
	splice(@a2);
	}  ## function ends
	###########################
	sub find_connected_chrs{
	my($chimerIndxRef,$combChrsRef,$a1Ref) = @_;
	my($str,$tstr);
	$str = $a1Ref->[$chimerIndxRef->{"chrs"}];  ## chromosome combination
	$tstr = '\(';
	$str =~ s/$tstr//;
	$tstr = '\)';
	$str =~ s/$tstr//;
	push(@{$combChrsRef},(split("-",$str)));
	}  ## function ends
	###########################
	sub find_break_points{
	my($chimerIndxRef,$breakPointsRef,$a1Ref) = @_;
		if($a1Ref->[$chimerIndxRef->{"strand1"}] eq "+"){
		push(@{$breakPointsRef},$a1Ref->[$chimerIndxRef->{"tEn1"}]);
		}
		else{
		push(@{$breakPointsRef},($a1Ref->[$chimerIndxRef->{"tSt1"}]+1));
		}			

		if($a1Ref->[$chimerIndxRef->{"strand2"}] eq "+"){
		push(@{$breakPointsRef},($a1Ref->[$chimerIndxRef->{"tSt2"}]+1));
		}
		else{
		push(@{$breakPointsRef},$a1Ref->[$chimerIndxRef->{"tEn2"}]);
		}
	}  ## function ends
	###########################	
	sub define_chimer_char_hash{
	my %fragGeneCharIndex = (
		"transcriptID" => 0,
		"chr" => 1,
		"txSt" => 2,
		"genicRegion" => 3,
		);
	return(%fragGeneCharIndex);
	}  ## function ends
	###########################
	sub fragment_gene_finder{
	my($chimerIndxRef,$knGnIndxRef,$chr,$algnPoint,$a1Ref,$gnIndxRef,$finRef,$geneCharsRef) = @_;

	## chimer array index, refgene record index, chromosome, chimer genomic coordinate, chimer array, reference transcript file index on chromosome, 
	## reference transcript file, chimer annotation info array

	my(@cluster,$i,$genicRegion);
		
	compartment_D::ChimerGeneFindSR::get_all_intersected_tracks($chr,$algnPoint,$knGnIndxRef,$gnIndxRef,$finRef,\@cluster);  ## all possible gene tracks

		if(scalar(@cluster)){  ## @cluster has all the possible genes
			for($i = 0;$i<scalar(@cluster);++$i){  ## for each gene encompassing the fusion point
			$genicRegion = "";
			$genicRegion = compartment_D::ChimerGeneFindSR::find_alignment_region($algnPoint,\@cluster,$i,$knGnIndxRef);
				if(!($genicRegion)){
				die "\n no genic region found for ",join("\t",@{$a1Ref}),"\n in ",$cluster[$i][$knGnIndxRef->{"transcriptID"}]," \n";
				}
			push(@{$geneCharsRef},[($cluster[$i][$knGnIndxRef->{"transcriptID"}],$cluster[$i][$knGnIndxRef->{"chr"}],$cluster[$i][$knGnIndxRef->{"txSt"}],$genicRegion)]);
			undef $genicRegion
			}  ## for(scalar(@cluster)) ends
		}  ## if(scalar(@cluster)) ends
		else{
		push(@{$geneCharsRef},[("none","none","none","intergenic")]);  ## if no gene found  (make a different function)
		}
	}  ## function ends
	###########################
	sub sort_frag_genes{  ## sorts array elements on priority of the region
	my($geneCharsIndxRef,$fragmentGenesRef) = @_;
	my($i,$j);	
	my %genicRegionIndex = compartment_D::ChimerGeneFindSR::generate_genic_region_indexes();

		for($i = 0;$i<(scalar(@{$fragmentGenesRef})-1);++$i){
			for($j = 0;$j<(scalar(@{$fragmentGenesRef})-1-$i);++$j){
				if($genicRegionIndex{$fragmentGenesRef->[$j+1][$geneCharsIndxRef->{"genicRegion"}]} >
					$genicRegionIndex{$fragmentGenesRef->[$j][$geneCharsIndxRef->{"genicRegion"}]}){
				shared::ArrayOperations::swap_elements($j,$j+1,scalar(keys(%{$geneCharsIndxRef})),$fragmentGenesRef);
				}  ## if(priotity2 > priority1) ends
			}  ## for($j) ends
		}  ## for($i) ends
	}  ## function ends
	###########################
	sub find_top_chimer_genes{
	my($geneCharsIndxRef,$fragmentGenesRef) = @_;
	my($i,$topReplicates);
	my %genicRegionIndex = compartment_D::ChimerGeneFindSR::generate_genic_region_indexes();
	$topReplicates = 0;
		for($i = 0;$i<scalar(@{$fragmentGenesRef});++$i){
			if($genicRegionIndex{$fragmentGenesRef->[$i][$geneCharsIndxRef->{"genicRegion"}]} == 
				$genicRegionIndex{$fragmentGenesRef->[0][$geneCharsIndxRef->{"genicRegion"}]}){
			++$topReplicates;
			}
			else{
			last;
			}
		}  ## for(scalar(@{$fragmentGenesRef})) ends
	return($topReplicates);
	}  ## function ends
	###########################
	sub sort_top_gene_ids{  ## function to be shared with other module (pslx2 gene annotation)
	my($n,$hashRef,$clustRef) = @_;
	my($sortIndex,$i,$j);
	$sortIndex = $hashRef->{"transcriptID"};
		for($i = 0;$i<($n-1);++$i){
			for($j = 0;$j<($n-1-$i);++$j){
				if($clustRef->[$j+1][$sortIndex] lt $clustRef->[$j][$sortIndex]){  ## comparing exon overlap
				shared::ArrayOperations::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
				}  ## if(ends)
			}  ## for($j) ends
		}  ## for($i) ends
	}  ## function ends
	###########################
	sub select_appropriate_genes{
	my($geneCharsIndxRef,$frag1GenesRef,$frag2GenesRef,$chrCombRef) = @_;
	my($firstArrayRef,$secondArrayRef,$i);
	my($flippingFlag,$index1,$index2);
	my %genicRegionIndex = compartment_D::ChimerGeneFindSR::generate_genic_region_indexes();
	$flippingFlag = 0;
		if($genicRegionIndex{$frag1GenesRef->[0][$geneCharsIndxRef->{"genicRegion"}]} >= 
			$genicRegionIndex{$frag2GenesRef->[0][$geneCharsIndxRef->{"genicRegion"}]}){
		$firstArrayRef = $frag1GenesRef;
		$secondArrayRef = $frag2GenesRef;
		}
		else{
		$secondArrayRef = $frag1GenesRef;
		$firstArrayRef = $frag2GenesRef;
		$flippingFlag = 1;
		}
	$index1 = 0;
	$index2 = 0;
	for($i = 0;$i<scalar(@{$secondArrayRef});++$i){
		if($firstArrayRef->[0][$geneCharsIndxRef->{"transcriptID"}] eq $secondArrayRef->[$i][$geneCharsIndxRef->{"transcriptID"}]){
		$index2 = $i;
		last
		}  ## if (transcriptID1 eq transcriptID2) ends
	}  ## for(scalar(@{$secondArrayRef})) ends
	
		if(!($flippingFlag)){
		return($index1,$index2);
		}
		else{
		return($index2,$index1);
		}
	}  ## function ends
	###########################
	sub determine_chimer_type{
	my($geneIndex1,$geneIndex2,$geneCharsIndxRef,$frag1GenesRef,$frag2GenesRef,$combChrRef,$a1Ref,$chimerIndxRef) = @_;
	my($chimerTag,%genicRegionIndex);
	%genicRegionIndex = compartment_D::ChimerGeneFindSR::generate_genic_region_indexes();

		if( ($frag1GenesRef->[$geneIndex1][$geneCharsIndxRef->{"genicRegion"}] eq "intergenic") || 
			($frag2GenesRef->[$geneIndex2][$geneCharsIndxRef->{"genicRegion"}] eq "intergenic") ){
			$chimerTag = "GeneDesertChimera";
		}
		else{  ## if not in gene desert
			if($frag1GenesRef->[$geneIndex1][$geneCharsIndxRef->{"transcriptID"}] eq 
			  $frag2GenesRef->[$geneIndex2][$geneCharsIndxRef->{"transcriptID"}]){
	
				if($a1Ref->[$chimerIndxRef->{"strand1"}] eq $a1Ref->[$chimerIndxRef->{"strand2"}]){
				$chimerTag = "ExonReorderingChimera";
				}
				else{
				$chimerTag = "ExonFlippingChimera";
				}
			}  ## if (same transcript IDs) ends
			else{
			$chimerTag = "IntergenicChimera";
			}
		}  ## else(not intergenic) ends
	return($chimerTag);
	}  ## function ends
	###########################
	sub print_chimer_annotation_to_file{
	my($foutName,$chimerTag,$geneIndex1,$geneIndex2,$geneCharsIndxRef,$frag1GenesRef,$frag2GenesRef,$combChrsRef,$a1Ref,$chimerIndxRef) = @_;
	my(@outArray);

	push(@outArray,($chimerTag,$a1Ref->[$chimerIndxRef->{"qID"}]));

	push(@outArray,("(".$combChrsRef->[0],$frag1GenesRef->[$geneIndex1][$geneCharsIndxRef->{"genicRegion"}]));
	push(@outArray,$frag1GenesRef->[$geneIndex1][$geneCharsIndxRef->{"transcriptID"}].")");

	push(@outArray,("(".$combChrsRef->[1],$frag2GenesRef->[$geneIndex2][$geneCharsIndxRef->{"genicRegion"}]));
	push(@outArray,$frag2GenesRef->[$geneIndex2][$geneCharsIndxRef->{"transcriptID"}].")");
	

	open FHPrintChimerAnnotation,">>$foutName";  ## output file
	print FHPrintChimerAnnotation join("\t",@outArray),"\n";
	close FHPrintChimerAnnotation;
	}  ## function ends
	###########################
	sub print_chimer_gene_stats{
	my($foutName,$chimerTag,$geneIndex1,$geneIndex2,$geneCharsIndxRef,$frag1GenesRef,$frag2GenesRef,$combChrsRef,$a1Ref,$chimerIndxRef) = @_;
	my(@annotArray);
	open FHPrintChimerGeneStats,">>$foutName";  ## output file
	
		if($frag1GenesRef->[$geneIndex1][$geneCharsIndxRef->{"transcriptID"}] ne "none"){  ## if not gene desert
			if($frag1GenesRef->[$geneIndex1][$geneCharsIndxRef->{"genicRegion"}] ne "INTRON"){
			push(@annotArray,$frag1GenesRef->[$geneIndex1][$geneCharsIndxRef->{"transcriptID"}]);


			push(@annotArray,$frag1GenesRef->[$geneIndex1][$geneCharsIndxRef->{"chr"}]);
			push(@annotArray,$frag1GenesRef->[$geneIndex1][$geneCharsIndxRef->{"txSt"}]);

			push(@annotArray,("NA","NA","NA","NA"));
			push(@annotArray,$frag1GenesRef->[$geneIndex1][$geneCharsIndxRef->{"genicRegion"}]);
			push(@annotArray,($chimerTag,$a1Ref->[$chimerIndxRef->{"qID"}]));
			print FHPrintChimerGeneStats join("\t",@annotArray),"\n";
			}  ## if(genicRegion ne INTRON) ends
		}  ## if(1 is not geneDesert ends)
			
	splice(@annotArray);
	
		if($frag2GenesRef->[$geneIndex2][$geneCharsIndxRef->{"transcriptID"}] ne "none"){
			if($frag2GenesRef->[$geneIndex2][$geneCharsIndxRef->{"genicRegion"}] ne "INTRON"){
			push(@annotArray,$frag2GenesRef->[$geneIndex2][$geneCharsIndxRef->{"transcriptID"}]);


			push(@annotArray,$frag2GenesRef->[$geneIndex2][$geneCharsIndxRef->{"chr"}]);
			push(@annotArray,$frag2GenesRef->[$geneIndex2][$geneCharsIndxRef->{"txSt"}]);

			push(@annotArray,("NA","NA","NA","NA"));
			push(@annotArray,$frag2GenesRef->[$geneIndex2][$geneCharsIndxRef->{"genicRegion"}]);
			push(@annotArray,($chimerTag,$a1Ref->[$chimerIndxRef->{"qID"}]));
			print FHPrintChimerGeneStats join("\t",@annotArray),"\n";
			}  ## if(genicRegion ne INTRON) ends
		}  ## if(2 is not geneDesert ends)

	close FHPrintChimerGeneStats;
	}  ## functin ends
	###########################
1;

