
package file_parsers::FormatConversionRules;

require Exporter;
@ISA = qw/Exporter/;
@EXPORT = qw/bed_to_default_conversion_rules/;

use strict;
use warnings;

	sub bed_to_default_conversion_rules{
	my %returnHash = (
		0 => 3,
		1 => 0,
		2 => 5,	
		3 => 1,	
		4 => 2,
		5 => 6,
		6 => 7,
		7 => 9,
		8 => 11,
		9 => 10,
		);
	return(%returnHash);
	}  ## function ends
	###########################
1;


