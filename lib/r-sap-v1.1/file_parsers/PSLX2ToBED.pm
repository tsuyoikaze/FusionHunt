package file_parsers::PSLX2ToBED;


use strict;
use warnings;

	sub pslx2_to_bed{
	my($pslx2IndxRef,@pslx2) = @_;

	my(@bed,@temp1,@temp2);
	my($i);

	push(@bed,$pslx2[$pslx2IndxRef->{"tID"}]);  ## chromosome name
	push(@bed,$pslx2[$pslx2IndxRef->{"tSt"}]);  ## chrSt
	push(@bed,$pslx2[$pslx2IndxRef->{"tEn"}]);  ## chrEn
	push(@bed,$pslx2[$pslx2IndxRef->{"qID"}]);  ## query name
	push(@bed,'2');  ## score(random)
	push(@bed,$pslx2[$pslx2IndxRef->{"strand"}]);  ## strand

	push(@bed,$pslx2[$pslx2IndxRef->{"tSt"}]);  ## chrSt (thick start)
	push(@bed,$pslx2[$pslx2IndxRef->{"tSt"}]);  ## chrSt (thick end)
	
	push(@bed,'0,0,0');  ## rgb (import from the color hash)
	


	push(@bed,$pslx2[$pslx2IndxRef->{"blNum"}]);  ## block count

	push(@bed,$pslx2[$pslx2IndxRef->{"blLen"}]);  ## block length
	
		
	@temp1 = split(',',$pslx2[$pslx2IndxRef->{"tBlSt"}]);  ## tblock starts

		for($i = 0;$i<scalar(@temp1);++$i){
		push(@temp2,($temp1[$i]-$pslx2[$pslx2IndxRef->{"tSt"}]));  ## relative to the txSt
		}

	push(@bed,join(',',@temp2).',');  ## block starts
	return(@bed);
	}  ## function ends
	###########################
1;

