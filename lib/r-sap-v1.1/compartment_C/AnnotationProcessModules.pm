
package compartment_C::AnnotationProcessModules;


use strict;
use warnings;

	sub define_char_array_field_hash{
	my($hashRef) = @_;
	$hashRef->{"clustIndex"} = 0; 
	$hashRef->{"transcriptID"} = 1; 
	$hashRef->{"exonOverlap"} = 2; 
	$hashRef->{"geneOverlap"} = 3; 
	$hashRef->{"exonDeletions"} = 4; 
	$hashRef->{"intronOverlap"} = 5; 
	$hashRef->{"geneExpansion"} = 6;
	$hashRef->{"geneDistance"} = 7;  ## for neighboring genes (within the permissible distance)
	$hashRef->{"codingFlag"} = 8;
	$hashRef->{"trueIntronOverlapFlag"} = 9;
	$hashRef->{"trueGeneExpansionFlag"} = 10;
	$hashRef->{"intronChainMatch"} = 11;
	}  ## function ends
	###########################
	sub find_top_replicates{  ## function to find out all the transcripts with same alignment characteristics as the top transcript
	my($hashRef,$clustRef) = @_;
	my($i,$j,$matchFlag,$topReplicates);
	$topReplicates = 0;
		for($i = 0;$i<scalar(@{$clustRef});++$i){
		$matchFlag = 0;
			for($j = 2;$j<scalar(keys(%{$hashRef}));++$j){
				if($clustRef->[$i][$j] == $clustRef->[0][$j]){
				$matchFlag += 1;
				}
			}  ## for(scalar(keys(%{$hashRef}))) ends
			if($matchFlag == (keys(%{$hashRef})-2)){
			$topReplicates += 1;
			}
		}  ## for(scalar(@{$clustRef})) ends
	
		if(!($topReplicates)){
		die "\n no top replicate found for $clustRef->[0][1]\n";
		}
	return($topReplicates);
	}  ## function ends
	###########################
	sub sort_top_gene_ids{
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
	sub decide_annotation_tag{
	my($hashRef,$clustRef,$currTrackIndx,$annotTagsRef) = @_;
	
		if($clustRef->[$currTrackIndx][$hashRef->{"geneOverlap"}]){  ## if there is any gene overlap
			if($clustRef->[$currTrackIndx][$hashRef->{"exonOverlap"}]){  ## if there is any exon overlap
				if( !($clustRef->[$currTrackIndx][$hashRef->{"trueIntronOverlapFlag"}]) && !($clustRef->[$currTrackIndx][$hashRef->{"trueGeneExpansionFlag"}]) && 
					!($clustRef->[$currTrackIndx][$hashRef->{"exonDeletions"}]) ){
				push(@{$annotTagsRef},"ExonOnly");
				}
				else{
					if($clustRef->[$currTrackIndx][$hashRef->{"exonDeletions"}]){
					push(@{$annotTagsRef},"ExonDeletion");
					}
					if($clustRef->[$currTrackIndx][$hashRef->{"trueIntronOverlapFlag"}]){
					push(@{$annotTagsRef},"ExonExtension");
					}
					if($clustRef->[$currTrackIndx][$hashRef->{"trueGeneExpansionFlag"}]){
					push(@{$annotTagsRef},"GeneExpansion");
					}
				}  ## if(!(exonOnly)) ends
			}  ## if(exonOverlap) ends
	
			if(!($clustRef->[$currTrackIndx][$hashRef->{"exonOverlap"}])){  ## if no overlap with any of the exon
				if($clustRef->[$currTrackIndx][$hashRef->{"trueIntronOverlapFlag"}]){  ## if overlap with introns
					if(!($clustRef->[$currTrackIndx][$hashRef->{"trueGeneExpansionFlag"}])){
					push(@{$annotTagsRef},"IntronOnly");
					}
					elsif($clustRef->[$currTrackIndx][$hashRef->{"trueGeneExpansionFlag"}]){
					##push(@{$annotTagsRef},"IntronOnly");
					push(@{$annotTagsRef},"GeneExpansion");
					}					
				}  ## if(intronsOverlap) ends
				else{  ##( no exon, no intron but gene overlap)
				push(@{$annotTagsRef},"Uncharacterized");
				}
			}  ## if(!(exonOverlap)) ends
			if(!(scalar(@{$annotTagsRef}))){  ## just in case
			push(@{$annotTagsRef},"Uncharacterized");
			}
		}  ## if(geneOverlap) ends
		else{  ## if(!(geneOverlap))
		push(@{$annotTagsRef},"NeighboringExon");
		}
	}  ## function ends
	###########################
	sub initialize_array{
	my($arrRef,$hashRef) = @_;
	my $i;
		for($i = 0;$i<scalar(keys %{$hashRef});++$i){
		$arrRef->[$i] = 0;
		}  ## for($i) ends
	}  ## function ends
	###########################
1;


