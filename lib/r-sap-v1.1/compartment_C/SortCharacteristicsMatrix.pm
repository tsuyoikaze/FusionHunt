package compartment_C::SortCharacteristicsMatrix;

use strict;
use warnings;

	sub sort_char_array{  ## sorting intersected transcripts on different characterstics
	my($hashRef,$clustRef) = @_;
	my($i,$j,@temp);	
	## note: all of the transcripts have same exon overlap

		for($i = 0;$i<(scalar(@{$clustRef})-1);++$i){
			for($j = 0;$j<(scalar(@{$clustRef})-1-$i);++$j){

				if($clustRef->[$j+1][$hashRef->{"intronChainMatch"}] > $clustRef->[$j][$hashRef->{"intronChainMatch"}]){  ## comparing exon overlap				
				shared::ArrayOperations::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
				}
				elsif($clustRef->[$j+1][$hashRef->{"intronChainMatch"}] == $clustRef->[$j][$hashRef->{"intronChainMatch"}]){
					if($clustRef->[$j+1][$hashRef->{"exonOverlap"}] > $clustRef->[$j][$hashRef->{"exonOverlap"}]){  ## comparing exon overlap
					shared::ArrayOperations::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
					}  ## if(ends)
					elsif($clustRef->[$j+1][$hashRef->{"exonOverlap"}] == $clustRef->[$j][$hashRef->{"exonOverlap"}]){  ## if overlaps with the gene boundaries are equal
						if($clustRef->[$j+1][$hashRef->{"geneOverlap"}] > $clustRef->[$j][$hashRef->{"geneOverlap"}]){
						shared::ArrayOperations::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
						}
						elsif($clustRef->[$j+1][$hashRef->{"geneOverlap"}] == $clustRef->[$j][$hashRef->{"geneOverlap"}]){  ## if overlaps with the gene boundaries are equal
							if($clustRef->[$j+1][$hashRef->{"exonDeletions"}] < $clustRef->[$j][$hashRef->{"exonDeletions"}]){
							shared::ArrayOperations::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
							}
							elsif($clustRef->[$j+1][$hashRef->{"exonDeletions"}] == $clustRef->[$j][$hashRef->{"exonDeletions"}]){  ## if overlaps with the gene boundaries are equal
								if($clustRef->[$j+1][$hashRef->{"intronOverlap"}] < $clustRef->[$j][$hashRef->{"intronOverlap"}]){
								shared::ArrayOperations::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
								}
								elsif($clustRef->[$j+1][$hashRef->{"intronOverlap"}] == $clustRef->[$j][$hashRef->{"intronOverlap"}]){  ## if overlaps with the gene boundaries are equal
									if($clustRef->[$j+1][$hashRef->{"geneExpansion"}] < $clustRef->[$j][$hashRef->{"geneExpansion"}]){
									shared::ArrayOperations::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
									}
									elsif($clustRef->[$j+1][$hashRef->{"geneExpansion"}] == $clustRef->[$j][$hashRef->{"geneExpansion"}]){  ## if overlaps with the gene boundaries are equal
										if($clustRef->[$j+1][$hashRef->{"geneDistance"}] < $clustRef->[$j+1][$hashRef->{"geneDistance"}]){
										shared::ArrayOperations::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
										}
										elsif($clustRef->[$j+1][$hashRef->{"geneDistance"}] == $clustRef->[$j+1][$hashRef->{"geneDistance"}]){  ## if distance from gene is equal
											if($clustRef->[$j+1][$hashRef->{"codingFlag"}] > $clustRef->[$j][$hashRef->{"codingFlag"}]){
											shared::ArrayOperations::swap_elements($j,$j+1,scalar(keys(%{$hashRef})),$clustRef);
											}
										}  ## elsif() ends
									}  ## elsif() ends
								}  ## elsif() ends					
							}  ## elsif() ends
						}  ## elsif() ends
					}  ## elsif() ends
				}  ## elsif() ends
			}  ## for($j) ends
		}  ## for($i) ends
	}  ## function ends
	###########################
1;