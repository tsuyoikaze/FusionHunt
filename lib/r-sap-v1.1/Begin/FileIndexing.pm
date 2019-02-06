use strict;
use warnings;

package Begin::FileIndexing;
	
        sub read_in_gene_file{
        my($geneFile) = @_;
	my($str,$fp,@geneFileArray);
	
	open FHReadGeneFile,$geneFile or die "\n can not open file $geneFile \n";
	
	$fp = tell FHReadGeneFile;
		while($str = <FHReadGeneFile>){
			if($str =~ m/^(\#+)/){  ## if header found
			$fp = tell FHReadGeneFile
			}
			else{
			last;
			}
		}  ## while(<FHReadGeneFile>) ends

	seek(FHReadGeneFile,$fp,0);
	@geneFileArray = <FHReadGeneFile>;
	close FHReadGeneFile;
	return(@geneFileArray);
	}  ## function ends
	###########################
	sub make_gene_track_index{
	my($chrNum,$fin) = @_;

	my %ref = ();
	my ($str,$index1,$index2);
	my(@a1,@a2);
	
	open FH21,$fin or die "\n can not open file input file $fin \n";  ## gene tracsk downloaded from UCSC
	$index1 = tell FH21;
	$str = <FH21>;
	
		if($str =~ m/\#/){
		$index1 = tell FH21;
		$str = <FH21>;
		}
	
	$str =~ s/\n//;
	$str =~ s/\r//;
	@a1 = split(/\t+/,$str);
	$index2 = tell FH21;
	
		while($str = <FH21>){
		$str =~ s/\n//;
		$str =~ s/\r//;
		@a2 = split(/\t+/,$str);
			if($a1[$chrNum] ne $a2[$chrNum]){
			$ref{$a1[$chrNum]} = $index1;
			@a1 = @a2;
			splice(@a2);
			$index1 = $index2;
			}
		$index2 = tell FH21;
		}  ## while(<FH21>) ends 
		
	$ref{$a1[$chrNum]} = $index1;
	close FH21;
	return(%ref);
	}  ## function ends
	###########################
	sub make_gene_track_index2{
	my($chrNum,$inRef) = @_;

	my %ref = ();
	my ($str,$index1,$index2);

	my(@a1,@a2,$i,$indxStart);
	
        $indxStart = 0;

                for($i = 1;$i<scalar(@{$inRef});++$i){
		@a1 = split("\t",$inRef->[$i-1]);
		@a2 = split("\t",$inRef->[$i]);

                        if($a1[$chrNum] ne $a2[$chrNum]){  ## if chromosome changes
                        $ref{$a1[$chrNum]} = $indxStart;
                        $indxStart = $i;
                        }
                }  ## for(@inRef) ends

	@a1 = split("\t",$inRef->[$i-1]);
       $ref{$a1[$chrNum]} = $indxStart;
        return(%ref);
        }  ## function ends
        ###########################
1;

