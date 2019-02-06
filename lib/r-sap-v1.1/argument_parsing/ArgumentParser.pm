package argument_parser::ArgumentParser;

use strict;
use warnings;

##require argument_parser::ArgumentErrorReporting;
require shared::CheckStringsNumericals;
require Control::FileIOModules;
require Control::MergeOutFiles;

our $VERSION = "1.1";

	sub parse_arguments{

		if(!(scalar(@ARGV))){
		report_error("NoArgument");
		}
		elsif( ($ARGV[0] eq "--h") || ($ARGV[0] eq "-h") ){  ## user wants help
		report_error("PrintShortHelp");
		}
		elsif( ($ARGV[0] eq "--v") || ($ARGV[0] eq "-v") ){
		report_error("PrintVersion");
		}
		elsif( int(scalar(@ARGV)/2) != scalar(@ARGV)/2){  ## missing value argument
		report_error("NoProperList");
		}
		
		

	my %arguments = (  ## arguments with their values

		"in1" => "",  ## blat output file (psl or pslx format)
		"in2" => "",  ## gene track file 
		"af" => "psl",   ## alignment file format (psl, plsx, gtf (from Cufflinks), bed (from Scripture))
		"rf" => "Default",  ## format of the gene track (standard all fields, BED or GTF)
		"outDir" => "",  ## output directory
		"outTag" => "",  ## prefix for output files
		"tNum" => 1,  ## number of threads to create
		"IdentityCutoff" =>  95,  ## identity cutoff
		"CovCutoff" =>  90,  ## coverage cutoff
		"gapTolerance" =>  20, ## gap tolerance to detect chimeric transcripts (fragmented alignment) 
		"LeftLeaveCutoff" =>  20,  ## gap tolerance for left end
		"RightLeaveCutoff" => 20,  ## gap tolerance for right end
		"DeletionCutoff" => 10,  ## deletion tolerances (to tolerate given number of bases)
		"GeneRadius" => 5000,  ## gene radius (to define gene desert in intergenic regions)
		"AnnotationMode" => "unique",  ## annotation stringency (multi/unique)
		"ExonExtensionCutoff" => 2,  ## tolerable number of bases extended in the intron or out of the gene boundary
		"IntronCutoff"=> 100000,  ## maximum number of base pairs by which alignmnet blokcs are allowd to be separated (maximum possible intron length)
		"stMatch" => 'F',  ## strand match between reference transcript and alignments
		);



		#"in1" =>"string",			## required
		#"in2" => "string",		## required
		#"af" => "psl/pslx/gtf/bed", 	## required (may be optional)
		#"rf" => "Default/bed/gtf",	## required (may be optional)
		#"outDir" => "string",			## required
		#"outTag" => "string",		## optional
		#"tNum" => "int",			## optional
		#"IdentityCutoff" =>  "float",  ## optional
		#"CovCutoff" =>  "float",	##optional
		#"gapTolerance" =>  "int",	## optional
		#"LeftLeaveCutoff" => "int",	## optional
		#"RightLeaveCutoff"=> "int",	## optional
		#"DeletionCutoff" => "int",	## optional
		#"GeneRadius" => "int",	## optional
		#"AnnotationMode" => "multi/unique",	## optional
		#"ExonExtensionCutoff" => "int",		## optional
		#"IntronCutoff"=> "int",			## optional
		#"stMatch" => "str",



	my %argumentTypes = (  ## arguments with their data type (float: floating number, int: integer, str: string, others: fixed values)
		"in1" =>"string",
		"in2" => "string",
		"af" => "psl/pslx/gtf/bed",
		"rf" => "Default/bed/gtf",
		"outDir" => "string",
		"outTag" => "string",
		"tNum" => "int",
		"IdentityCutoff" =>  "float",
		"CovCutoff" =>  "float",
		"gapTolerance" =>  "int",
		"LeftLeaveCutoff" => "int",
		"RightLeaveCutoff"=> "int",
		"DeletionCutoff" => "int",
		"GeneRadius" => "int",
		"AnnotationMode" => "multi/unique",
		"ExonExtensionCutoff" => "int",
		"IntronCutoff"=> "int",
		"stMatch" => "t/T/f/F/True/False",
		);

	my($i,$scriptValue,$listFlag,@temp);
	my $dumpFlag = 0;


		for($i = 0;$i<scalar(@ARGV);++$i){

		$ARGV[$i] =~ s/^(\s+)//;  ## argument name
		$ARGV[$i] =~ s/(\s+)$//;
		$ARGV[$i] =~ s/^(\-+)//;

		$ARGV[$i+1] =~ s/^(\s+)//;  ## argument value
		$ARGV[$i+1] =~ s/(\s+)$//;
		$ARGV[$i+1] =~ s/^(\-+)//;

		push(@temp,$ARGV[$i]);
		push(@temp,$ARGV[$i+1]);
		++$i;

		$listFlag = 0;
		##@temp = split(/\=/,$ARGV[$i]);  ## separating argument name and it's specified value

			if($temp[0] =~ m/((^dump)$)/i){
			$dumpFlag = 1;
			}
			else{

			#($listFlag,$scriptValue) = check_if_in_list($temp[0],\@argLists);  ## matching the argument in the list of acceptable arguments
			($listFlag,$scriptValue) = check_if_in_list($temp[0],\%argumentTypes);  ## matching the argument in the list of acceptable arguments

				if($listFlag){  ## if argument is acceptable
					if(scalar(@temp) != 2){
					report_error("MissingValue",$temp[0]);
					}
	
					if(check_data_type($temp[1],$scriptValue,\%argumentTypes)){  ## checking if type match for argument's value
					$arguments{$scriptValue} =$temp[1];  ## accpeting argument value
					}  ##  if(data type match) ends
					else{
					report_error("ArgTypeMismatch",$temp[0],$argumentTypes{$temp[0]});
					}
				}  ## if(listFlag) ends
				else{
				report_error("InvalidArg",$temp[0]);
				}

			}  ## else (no dump option) ends
		splice(@temp);
		}  ## for(@ARGV) ends


		if( !($arguments{"in1"}) || !($arguments{"in2"}) || !($arguments{"outDir"})){
		report_error("MissingRequired");
		}
		
		#### output directory check ####
		$arguments{"outDir"} =~ s/\s//g;

			if($arguments{"outDir"}){
			$arguments{"outDir"} = polish_directory_path($arguments{"outDir"});  ## checking for the correct path
			}
			else{
			my $tempDirTag = Control::FileIOModules::get_temporary_dir_tag();
			mkdir($tempDirTag);
			$arguments{"outDir"} = $tempDirTag."/";
			}

	return($dumpFlag,%arguments);
	}  ## function ends
	###########################
	sub check_if_in_list{
	my($value,$argsHashRef) = @_;
	my($i,$scriptValue,$matchFlag);
	$matchFlag = 0;
	$scriptValue = "";

		for my $key1(keys(%{$argsHashRef})){
			if($key1 =~ m/$value/i){
			$matchFlag = 1;
			$scriptValue = $key1;
			last;
			}
		}  ## for(@{$argsList}) ends

	return($matchFlag,$scriptValue);
	}  ## function ends
	###########################
	sub check_data_type{
	my($value,$key,$argsHashRef) = @_;
	my(@temp,$i);
	my $dataTypeFlag = 0;
		if($argsHashRef->{$key} eq "int"){
		$dataTypeFlag =  shared::CheckStringsNumericals::check_if_integer($value);
		}
		elsif($argsHashRef->{$key} eq "float"){
		$dataTypeFlag =  shared::CheckStringsNumericals::check_if_float($value);
		}
		elsif($argsHashRef->{$key} eq "string"){
		$dataTypeFlag =  shared::CheckStringsNumericals::check_if_string($value);
		}
		else{  ## predefined values
		@temp = split(/\//,$argsHashRef->{$key});
			for($i = 0;$i<scalar(@temp);++$i){

				##if($temp[$i] eq $value){
				if($temp[$i] =~ m/$value/i){
				$dataTypeFlag = 1;
				last;
				}
			}  ## for(scalar(@temp)) ends
		}  ## else ends
	return($dataTypeFlag);
	}  ## function ends
	###########################
	sub polish_directory_path{
	my($path) = @_;
	$path =~ s/\s//g;

	$path =~ s/(\/)$//;
	$path =~ s/(\\)$//;
	
	my $pathVar = Control::MergeOutFiles::get_path_var();
		if($path){
		$path .= $pathVar;
		}
	return($path);
	}  ## function ends
	###########################
	sub report_error{
	my($signal,$argName,$argType) = @_;
		if($signal eq "NoArgument"){
		print "\nNo argumnet passed";
		}
		if($signal eq "NoProperList"){
		print "\nArgument list is not proper (may be misssing values for some parameters)";
		}
		if($signal eq "MissingValue"){
		print "\nMissing value for the argument $argName";
		}
		if($signal eq "ArgTypeMismatch"){
		print "\nArgument value mismatch for $argName ($argType required)";
		}
		if($signal eq "InvalidArg"){
		print "\nInvalid argument $argName";
		}
		if($signal eq "MissingRequired"){
		print "\nMissing the required argument(s)";
		}
		if($signal eq "PrintShortHelp"){
		##
		}
		if($signal eq "PrintVersion"){
		print "Version: $VERSION\n";
		exit 1;
		}

		print_help_summary();
		exit 1;
	}  ## function ends
	###########################
	sub print_help_summary{
	
	print "\nProgram: R-SAP (RNA-Seq analysis pipeline)\nVersion: $VERSION";
	print "\n\nUsage:  perl R-SAP --in1 <file_name> --in2 <file_name> --outDir <dir_name> <options>
	
	Where:
	======
	in1       Reference genome alignment file of the RNA-Seq reads in psl or pslx format (required).
	in2       Reference transcript annotation file (required).
	outDir    A directory where all the output files will be written (required). 
	

	Options:
	-------
	stMatch   Match reference transcript and alignment strand (T/F, default is F(false)).
	outTag    Prefix for output file names (optional but recommended).
	tNum      Number of parallel threads to be created (optional, default is 1).
	rf	  Reference transcript file format (required only for BED and GTF files).
	af	  Alignment file format (required only for GTF and BED files)
	
	For other optional arguments/parameters and their details, please see the README file in the program directory of R-SAP\n\n";
	}  ## function ends
	###########################

1;

## argument lists from the user ##
##input files : required
	#aligment file: psl/pslx format
	#geneTrack file: 
	#geneTrackFileFormat: default/BED/GTF


##setting : optional
	#identity  %
	#coverage  %
	#tolerateGap int
	#tolerateLeftGap int
	#tolerateRightGap int
	#deletionTolerance int
	#geneExpansion bases


##output option : required/optional
	#outputFilesPrefix string/numeral : required
	#output directory: optional
	#dump/noDump: string : optional

##threads : optional
	#threads: int
## argument parsing procedure
	#if no arguments passed: print_help_summary()
	#if required missing: print_required_needed();
	

