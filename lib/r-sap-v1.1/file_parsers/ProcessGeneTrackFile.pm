package file_parsers::ProcessGeneTrackFile;

use strict;
use warnings;

use file_parsers::AllFieldsToDefault;
use file_parsers::BEDToDefault;
use file_parsers::GTFToDefault;

	sub process_gene_tracks{
	my($geneFile,$outFile,$fileFormat) = @_;

		if($fileFormat =~ m/(default)/i){
		file_parsers::AllFieldsToDefault::all_fields_to_default($geneFile,$outFile);
		}
		elsif($fileFormat =~ m/(bed)/i){
		file_parser::BEDToDefault::bed_to_default($geneFile,$outFile);
		}
		elsif($fileFormat =~ m/(gtf)/i){
		file_parsers::GTFToDefault::gtf_to_default($geneFile,$outFile);
		## do something
		}
		else{
		die "\n specified file format $fileFormat was not recognized \n";
		}
	}  ## function ends
	###########################
1;

