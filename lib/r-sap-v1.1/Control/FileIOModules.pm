package Control::FileIOModules;

use strict;
use warnings;

	sub get_temporary_dir_tag{
	my $localTime = localtime();
	$localTime =~ s/\D+//g;
	return($localTime);
	}  ## function ends
	###########################
	sub get_temporary_input_dir{
	my($tempDirTag) = @_;
	return($tempDirTag."InFileParts");
	}  ## function ends
	###########################
	sub get_temporary_output_dir{
	my($tempDirTag) = @_;
	return($tempDirTag."OutFileParts");
	}  ## function ends
	###########################
1;
