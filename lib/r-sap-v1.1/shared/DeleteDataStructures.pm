package shared::DeleteDataStructures;


	##sub undefine_arrays{
	##sub undefine_hashes{
	##sub undefine_scalars{
	
	sub undefine_arrays{
		foreach my $ref(@_){
		undef @{$ref};
		}
	}  ## function ends
	###########################
	sub undefine_hashes{
		foreach my $ref(@_){
		undef %{$ref};
		}
	}  ## function ends
	###########################
	sub undefine_scalars{
		foreach my $ref(@_){
		undef ${$ref};
		}
	}  ## function ends
	###########################

1;
