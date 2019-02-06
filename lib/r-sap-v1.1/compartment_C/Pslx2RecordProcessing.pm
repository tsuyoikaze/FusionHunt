package compartment_C::Pslx2RecordProcessing;

use strict;
use warnings;


	
	sub get_pslx2_algn_blocks{
	my($pslx2IndxRef,$a1Ref,$exStRef,$exEnRef) = @_;
	my($i);
	@{$exStRef} = split(',',$a1Ref->[$pslx2IndxRef->{"tBlSt"}]);  ## starting positions of blocks on  chromosome
	@{$exEnRef} = split(',',$a1Ref->[$pslx2IndxRef->{"blLen"}]);  ## length of the alignment blocks

		for($i = 0;$i<scalar(@{$exStRef});++$i){
		$exEnRef->[$i] += $exStRef->[$i]  ## calculating block end positions on chromosome
		}
	}  ## function ends
	###########################
1;

