package file_parsers::ProcessAlignmentFile;

use strict;
use warnings;

use file_parsers::GTFToPslx;
use file_parsers::BEDToPslx;

	sub process_alignment_file{
	my($alignFile,$outFile,$fileFormat) = @_;

		if($fileFormat =~ m/(gtf)/i){
		file_parsers::GTFToPslx::gtf_to_pslx($alignFile,$outFile);
		}
		elsif($fileFormat =~ m/(bed)/i){
		file_parsers::BEDToPslx::bed_to_pslx($alignFile,$outFile);
		}
		else{
		die "\n alignment file format $fileFormat is not allowed (terminating)\n";
		}
	}  ## function ends
	###########################
1;

