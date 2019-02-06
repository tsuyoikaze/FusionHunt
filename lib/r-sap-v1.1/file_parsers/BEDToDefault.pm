package file_parsers::BEDToDefault;

use strict;
use warnings;

use file_parsers::CheckIfBED;
use file_parsers::FormatConversionRules;

	sub bed_to_default{
	my($geneFile,$outFile) = @_;
	open FHBedToDefaultIn,$geneFile or die "\n can not open file $geneFile \n";  ## BED format input file
	open FHBedToDefaultOut,">$outFile";  ## output file (default format)

	my($str,@a1,@a2);
	my($key1,$i,@temp1,@temp2);
	
	my %bedToDefaultConversion = bed_to_default_conversion_rules();
	
		while($str = <FHBedToDefaultIn>){
		$str =~ s/\n//;
		$str =~ s/\r//;
			if( ($str =~ m/^(\#+)/) || !($str) ){  ## skipping the header blank line
			next;
			}
		@a1 = split(/\t+/,$str);
			if(check_if_bed(@a1)){  ## if record is in bed format
			@a2 = change_bed_to_default_array(\@a1,\%bedToDefaultConversion);
			print FHBedToDefaultOut join("\t",@a2),"\n";
			}
			else{
			die "\n file format not recognized as BED at line $. in gene track file $geneFile \n";
			}
		splice(@a1);
		splice(@a2);
		}  ## while(<FHBedToDefaultIn>) ends
	
	close FHBedToDefaultIn;
	close FHBedToDefaultOut;
	}  ##  function ends
	###########################
	sub change_bed_to_default_array{
	my($a1Ref,$conversionRef) = @_;  ## array in BED format and conversion rule
	my(@a2,$i);
	my(@temp1,@temp2);
		for my $key1 (sort(keys(%{$conversionRef}))){
		$a2[$key1] = $a1Ref->[$conversionRef->{$key1}];
		}
	@temp1 = split(',',$a2[8]); ## exon starts (relative to transcription start position)
	@temp2 = split(',',$a2[9]);  ## exon lengths
		for($i = 0;$i<scalar(@temp1);++$i){
		$temp1[$i] += $a2[3];
		$temp2[$i] += $temp1[$i];
		}
	$a2[8] = join(',',@temp1);  ## exon start positions
	$a2[9] = join(',',@temp2);  ## exon end positions
	push(@a2,'.');  ## blank feature for other gene feature
	splice(@temp1);
	splice(@temp2);
	return(@a2);
	}  ## function ends
	###########################
1;	
