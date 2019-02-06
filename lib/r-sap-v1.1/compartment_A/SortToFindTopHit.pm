################ Package to sort pslx2 records based on their BLAT score, identity ####################
################ (if scores are equal) or coverage (if identities are equal). Input ###################
################ is pslx file sorted based on IDs. ####################################################

package compartment_A::SortToFindTopHit;

use strict;
use warnings;

require Begin::Indexes;

	sub main_module_SortToFindTopHit{
	my($pslx2IndxRef,$clustRef) = @_;
	my($i,@cluster);

		for($i = 0;$i<scalar(@{$clustRef});++$i){
		push(@cluster,join("\t",@{$clustRef->[$i]}[0..(scalar(keys(%{$pslx2IndxRef}))-1)]));
		}
	@cluster = merge_sort($pslx2IndxRef,@cluster);
	undef(@{$clustRef});
		for($i = 0;$i<scalar(@cluster);++$i){
		my @temp = split(/\t+/,$cluster[$i]);
		push(@{$clustRef},[@temp]);
		splice(@temp);
		}
	}  ## function ends	
	###########################

	########### merge sort algorithm with equality conditional cases ###########
 	sub merge_sort{
 	my ($pslx2IndxRef,@a1) = @_;  ## passed array containing pslx2 record
 	my(@left,@right,@result);  ## temporary arrays
	my($middle,$i,$j,$l);
	my(@temp_array,$temp_len);
 	$l = scalar(@a1);  ## length of the pslx2 array
 	
 		if($l <= 1){
 		return (@a1);
 		}

 		else{
 		$middle = get_middle_point($l);
			for($i = 0;$i<=$middle;++$i){
			push(@left,$a1[$i]);
			}
			for($i = ($middle+1);$i<$l;++$i){
			push(@right,$a1[$i]);
			}
		@left = merge_sort($pslx2IndxRef,@left);
		@right = merge_sort($pslx2IndxRef,@right);
		@result = merge( ($pslx2IndxRef,join("\n",@left)) ,(join("\n",@right)) );
		return(@result);
		}
	}  ## function ends
	###########################	 		

	########### Inner module 1 for the merge sort algorithm. Prioritize array parts to sorted ###########
	sub merge{
	my($pslx2IndxRef,$str1,$str2) = @_;
	my(@left,@right,@result,@tempLeft,@tempRight);
	my($lenl,$lenr,$pref);
	my($scoreIndx,$coverageIndx,$identityIndx);

	##$scoreIndx = Begin::Indexes::pslx2_indx("score");
	##$coverageIndx = Begin::Indexes::pslx2_indx("coverage");
	##$identityIndx = Begin::Indexes::pslx2_indx("identity");
	
	@left = split(/\n/,$str1);
	@right = split(/\n/,$str2);

	$lenl = @left;
	$lenr = @right;
	
		while( ($lenl) && ($lenr) ){
		@tempLeft = split(/\t+/,$left[0]);
		@tempRight = split(/\t+/,$right[0]);
		$pref = 1;

			if($tempLeft[$pslx2IndxRef->{"score"}]  > $tempRight[$pslx2IndxRef->{"score"}]){
			$pref = 0;
			}
			elsif($tempLeft[$pslx2IndxRef->{"score"}]  == $tempRight[$pslx2IndxRef->{"score"}]){
				if($tempLeft[$pslx2IndxRef->{"identity"}]  > $tempRight[$pslx2IndxRef->{"identity"}]){
				$pref = 0;
				}
				elsif($tempLeft[$pslx2IndxRef->{"identity"}]  == $tempRight[$pslx2IndxRef->{"identity"}]){
					if($tempLeft[$pslx2IndxRef->{"coverage"}]  > $tempRight[$pslx2IndxRef->{"coverage"}]){
					$pref = 0;
					}
				} ## elsif(2) ends
			}  ## elsif(1) ends, all elsif ends
				
			if($pref == 0){
			push(@result,$left[0]);
			shift(@left);
			}
			elsif($pref == 1){
			push(@result,$right[0]);
			shift(@right);
			}

		splice(@tempLeft);
		splice(@tempRight);
		$lenl = @left;
		$lenr = @right;
		}  ## while(lenl and lenr) ends 
		
		$lenl = @left;
		$lenr = @right;
		

		if($lenl > 0){
		push(@result,@left[0..($lenl-1)]);
		}
		if($lenr > 0){
		push(@result,@right[0..($lenr-1)]);
		}
	return(@result);
	
	}  ## function ends 
	###########################

	########### Inner module 2 for the merge sort algorithm. Finds middle based on the number ###########
	########### of elements (pslx2 records)in the array to be sorted. ###################################
	sub get_middle_point{
	my($l) = @_;
	my $middle;
	
		if( ($l/2) > int($l/2) ){
		$middle = int($l/2)+1-1;
		}
		elsif( ($l/2) == int($l/2) ){
		$middle = ($l/2)-1;
		}
	
	return($middle);
	}  ##  function ends
	###########################

1;

