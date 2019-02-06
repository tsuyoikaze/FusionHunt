package compartment_E::CalculateGeneExpressionSR;

use strict;
use warnings;

require compartment_C::FindGenicRegions;
require shared::CheckStringsNumericals;

	###########################
	sub collect_ref_gene_info{  ## storing the gene IDs and other info from the provided gene annotation track
	my($geneFile,$knGnIndxRef) = @_;

	my $transcriptKey;  ## transcript_id + chromosome + transcription start

	my($str,@a1,@temp1,@temp2,@tempClust,$i);
	my($tranID,$chr,$txSt,$txEn,$tranLen,$codingTag,$readCount);
	my($exonReadCounts,$exonLengths,$skippedExons);
	my($intronReadCounts,$intronLengths,$retainedIntrons);
	my $multiAnnotFlag = 0;
	
	my($arrayRef,%geneInfoHash);
	
	open FHDGEIn,$geneFile or die "\n can not open file $geneFile \n";  ## gene track file

		while($str = <FHDGEIn>){  ## creating a hash to store gene information
		$str =~ s/\n//;
		$str =~ s/\r//;
		next if($str =~ m/^(\#+)/);  ## skipping the header


		@a1 = split(/\t+/,$str);

		$transcriptKey = $a1[$knGnIndxRef->{"transcriptID"}]."____".$a1[$knGnIndxRef->{"chr"}]."____".$a1[$knGnIndxRef->{"txSt"}];  ## unique ID for each reference transcript

		push(@tempClust,[@a1]);  ## to pass it the array to the following subroutines		

		$tranID = $a1[$knGnIndxRef->{"transcriptID"}];  ## transcript ID
		$chr = $a1[$knGnIndxRef->{"chr"}];
		$txSt = $a1[$knGnIndxRef->{"txSt"}];
		$txEn = $a1[$knGnIndxRef->{"txEn"}];
		$tranLen = compartment_C::FindGenicRegions::get_transcript_length(\@tempClust,0,$knGnIndxRef);  ## import subroutine
		$codingTag = compartment_C::FindGenicRegions::check_if_protein_coding(\@tempClust,0,$knGnIndxRef);  ## import subroutine
		$readCount = 0;

		## exon related info
		$exonReadCounts = initialize_exon_read_counts(\@a1,$knGnIndxRef);
		$exonLengths = get_exon_lengths(\@a1,$knGnIndxRef);  ## import subroutines (add module)
		$skippedExons = "";
		
		## intron related info
		$intronReadCounts = initialize_intron_read_counts(\@a1,$knGnIndxRef);
		$intronLengths = get_intron_lengths(\@a1,$knGnIndxRef);
		$retainedIntrons = "";

		$arrayRef = [];  ## referencing an anonymous array
		push(@{$arrayRef},($tranID,$chr,$txSt,$txEn,$tranLen,$codingTag,$readCount,$exonReadCounts,$exonLengths,
			$skippedExons,$intronReadCounts,$intronLengths,$retainedIntrons,"","",$multiAnnotFlag));  ## anonymous array reference

		$geneInfoHash{$transcriptKey} = $arrayRef;  ## creating a hash for the transcript id and the information above pair
		splice(@tempClust);
		}  ## while(<FHDGEIn>) ends

	close FHDGEIn;
	return(%geneInfoHash);
	}  ## function ends
	###########################
	sub initialize_exon_read_counts{
	my($a1Ref,$knGnIndxRef) = @_;
	my($exonNums,$i);
	$exonNums = "";
		for($i = 1;$i<=$a1Ref->[$knGnIndxRef->{"exCount"}];++$i){  ## number of exons
		$exonNums .= "0".',';
		}
	return($exonNums);
	}  ## function ends
	###########################
	sub get_exon_lengths{
	my($a1Ref,$knGnIndxRef)= @_;
	my(@exSt,@exEn,@exonLengths,$i,@tempClust);
	push(@tempClust,[@{$a1Ref}]);
	compartment_C::FindGenicRegions::get_exon_cords(\@tempClust,0,$knGnIndxRef,\@exSt,\@exEn);  ## getting exon coordinates;
	splice(@exonLengths);
	
		for($i = 0;$i<scalar(@exSt);++$i){
		push(@exonLengths,($exEn[$i]-$exSt[$i]));
		}

		if($a1Ref->[$knGnIndxRef->{"strand"}] eq '-'){  ## if strand is negative
		@exonLengths = reverse(@exonLengths);
		}

	splice(@exSt);
	splice(@exEn);
	return(join(',',@exonLengths));
	}  ## function ends
	###########################
	sub initialize_intron_read_counts{
	my($a1Ref,$knGnIndxRef) = @_;
	my($intronNums,$i);
	$intronNums = "";
		for($i = 1;$i<=($a1Ref->[$knGnIndxRef->{"exCount"}]-1);++$i){  ## for(number of exons-1)
		$intronNums .= "0".',';
		}
	return($intronNums);
	}  ## function ends
	###########################
	sub get_intron_lengths{
	my($a1Ref,$knGnIndxRef)= @_;
	my(@inSt,@inEn,@intronLengths,$i,@tempClust);
	my(@exSt,@exEn);

	push(@tempClust,[@{$a1Ref}]);

	compartment_C::FindGenicRegions::get_exon_cords(\@tempClust,0,$knGnIndxRef,\@exSt,\@exEn);  ## getting exon coordinates;
	compartment_C::FindGenicRegions::get_intron_regions(\@exSt,\@exEn,\@inSt,\@inEn);  ## getting exon coordinates;

	splice(@intronLengths);
		for($i = 0;$i<scalar(@inSt);++$i){
		push(@intronLengths,($inEn[$i]-$inSt[$i]));
		}
		if($a1Ref->[$knGnIndxRef->{"strand"}] eq '-'){  ## if strand is negative
		@intronLengths = reverse(@intronLengths);
		}
	return(join(',',@intronLengths));
	}  ## function ends
	###########################
	sub increment_read_counts{
	my($a1Ref,$a2Ref) = @_;
	my($i);
		for($i = 0;$i<scalar(@{$a1Ref});++$i){
		++$a2Ref->[($a1Ref->[$i]-1)];
		}
		
	}  ## function ends
	###########################
	sub calculate_RPKM{
	my($l,$c,$n) = @_;
	my $expValue;
		if(($n) && ($l)){
		$expValue = ((10**9)*($c))/($n*$l);
		}
		else{
		$expValue = 0;
		}
	return($expValue);
	}  ## function ends
	###########################
	sub check_if_same_order{
	my($num1,$num2) = @_;
	my $flag = 0;

		if( ($num1 > 0) && ($num2 > 0) ){
			if( (($num1 >= 1) && ($num2 >= 1)) || (($num1 < 1) && ($num2 < 1)) ){
				if(shared::CheckStringsNumericals::order_of_magnitude($num1) == shared::CheckStringsNumericals::order_of_magnitude($num2)){
				$flag = 1;
				}
			}  ## if() ends
		}  ## if($num1 && $num2) ends
	return($flag);
	}  ## function ends
	###########################	
	sub check_if_multi_splicing{
	my($tagList) = @_;
	my(@temp,$multiAnnotFlag);
	@temp = split(',',$tagList);
	$multiAnnotFlag = 0;
	shared::ArrayOperations::remove_array_element(\@temp,"ExonOnly");
	shared::ArrayOperations::remove_array_element(\@temp,"NeighboringGenes");
	@temp = shared::ArrayOperations::remove_array_redundancy(@temp);

		if(scalar(@temp) > 1){
		$multiAnnotFlag = 1;
		}
		
	return($multiAnnotFlag);
	}  ## function ends
	###########################		
1;
