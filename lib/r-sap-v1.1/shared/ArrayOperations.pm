package shared::ArrayOperations;

## shift	sub find_min{
## shift	sub swap_elements{  ## function to swap two rows  (n1 with n2) in matrix (reference: clustRef)
## shift	sub initialize_array{
## shift	sub remove_array_redundancy{  ## returns sorted non-redundant array
## shift	sub check_array_length{  ##3


	########### Finds the minimum value in a given array. All the array elements must be numeric ###########
	sub find_min{
	my(@a1) = @_;
	my @a2 = sort {$a <=> $b} @a1;
	return($a2[0]);
	}  ## function ends 
	###########################
	sub swap_elements{  ## function to swap two rows  (n1 with n2) in matrix (reference: clustRef)
	my($n1,$n2,$l,$clustRef)= @_; 
	my(@temp);
        @temp = @{$clustRef->[$n1]}[0..($l-1)];
        $clustRef->[$n1] = [@{$clustRef->[$n2]}[0..($l-1)]];
        $clustRef->[$n2] = [@temp];
        }  ## function ends
	###########################
	sub initialize_array{  ## rerutn(void)
	my($arrayRef,$n) = @_;
	my $i;
		for($i = 0;$i<$n;++$i){
		$arrayRef->[$i] = 0;
		}
	}  ## function ends		
	############################
	sub initialize_value_array{  ## return(@array)
	my($str,$len) = @_;
	my(@a1,$i);
		for($i = 0;$i<$len;++$i){
		push(@a1,$str);
		}  ## for($i) ends
	return(@a1);
	}  ## function ends
	###########################
	sub remove_array_redundancy{  ## returns sorted non-redundant array
	my(@a1) = @_;
	my($element,%hash);
		for $element (@a1){
		$hash{$element} = 1;
		}
	@a1 = keys(%hash);
	return(@a1);
	}  ## function ends
	###########################
	########### Checks if the length of the passed array contains the desured number of elements ###########
	########### Error will terminate the program ###########

	sub check_array_length{  ##3
	my($ArrayLength,@a1) = @_;
		if(scalar(@a1) != $ArrayLength){
		print "\n length of the psl record is not appropriate for \n\n";
		print join("\t",@a1);
		die "\n\n";
		}
	}  ## function ends
	###########################
	sub remove_array_element{
	my($a1Ref,$element) = @_;
	my($flag,$i);
	$flag = 0;
		for($i = 0;$i<scalar(@{$a1Ref});++$i){
			if($a1Ref->[$i] eq $element){
			splice(@{$a1Ref},$i,1);
			$flag = 1;
			last;
			}
		}  ## for($i) ends
	return($flag);
	}  ## function ends
	###########################
1;
