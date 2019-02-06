package Control::MergeOutFiles;

use strict;
use warnings;

	sub stream_out_temp_dir{
	my($outDir,$outFileTag, $tempOutDir, $numFileParts,@outHashRefs) = @_;

	my ($tempOutFileName);
	my ($windowsFlag,$i,$j,$pathVar);
	$pathVar = get_path_var();

	$windowsFlag = check_if_windows();  ## checking if system is windows (returns 0 if not);

		for($i = 0;$i<scalar(@outHashRefs);++$i){  ## for each final output hash reference stored in the array

			for my $key1 (keys(%{$outHashRefs[$i]})){  ## for each element (output file name) in the final hash

				if(!($windowsFlag)){  ## if system is non-windows
				$tempOutFileName = $outHashRefs[$i]->{$key1};  ## output file name (final)
				my $tempTag = $outDir.$outFileTag;
				$tempOutFileName =~ s/^($tempTag)//;  ## removing the output file tag
				$tempOutFileName = $tempOutDir.$pathVar."*".$tempOutFileName;  ## temporary file name
				system("cat $tempOutFileName >> $outHashRefs[$i]->{$key1}");
				}

				else{  ## if system is windows

					for($j = 1;$j<=$numFileParts;++$j){  ## for each thread run
					$tempOutFileName = $outHashRefs[$i]->{$key1};  ## final output file

					my @temp = split(/\\/,$tempOutFileName);
					$tempOutFileName = $temp[scalar(@temp)-1];  ## file name only (no tags)

					$tempOutFileName =~ s/^($outFileTag)//;  ## removing the outFilePrefix
					$tempOutFileName = $tempOutDir.$pathVar.$j.$tempOutFileName;

					##print "\ntempOutFileName: $tempOutFileName\n";					

					append_file($tempOutFileName,$outHashRefs[$i]->{$key1},$outDir);  ## inner dir to final out dir
					}  ## for(each thread) ends
				} ## elsif(windows system) ends

			}  ## for each hash ends
		}  ## for(@outHashRefs) ends
	}  ## function ends
	###########################
	sub append_file{
	my($file1,$file2,$outDir) = @_;  ## append file1 to file2

	open FH1AppendFile,$file1 or die "\n can not open file $file1, check if $outDir still exists \n";
	open FH2AppendFile,">>".$file2;
	my $str;
		while($str = <FH1AppendFile>){
		print FH2AppendFile $str;
		}  ## while($str = <FH1Append>) ends
	close FH1AppendFile;
	close FH2AppendFile;
	}  ## function ends
	###########################


	#### BED file sub-routines ####
	sub merge_create_bed_files{

	my($outDir, $outFileTag, $tempOutDir, $numFileParts,$fileNameHashRef,$bedOutFile) = @_;

	my($tempOutFileName,$windowsFlag);
	my(%colorHash,$trackLine);
	my $pathVar = get_path_var();

	$windowsFlag = check_if_windows();  ## checking if system is windows (returns 0 if not);
	%colorHash = browser_track_color_hash();  ## RGB values for the browser track colors


		for my $key (keys(%{$fileNameHashRef})){  ## for all the final BED file names
		$trackLine = generate_track_name($key,$colorHash{$key});
	
		##### printing track description and track formating ####
		open FHPrintBEDTracks,">>".$bedOutFile;  ## final output file
		print FHPrintBEDTracks "\n$trackLine\n";
		close FHPrintBEDTracks;
		######################

			if(!($windowsFlag)){  ## if system is non-windows
			$tempOutFileName = $fileNameHashRef->{$key};
			my $tempTag = $outDir.$outFileTag;
			$tempOutFileName =~ s/^($tempTag)//;  ## removing the output file tag
			$tempOutFileName = $tempOutDir.$pathVar."*".$tempOutFileName;
			system("cat $tempOutFileName >> $bedOutFile");
			}

			else{  ## if system is windows
				for(my $i = 1;$i<=$numFileParts;++$i){  ## for each thread run
				$tempOutFileName = $tempOutFileName = $fileNameHashRef->{$key};

				my @temp = split(/\\/,$tempOutFileName);
				$tempOutFileName = $temp[scalar(@temp)-1];  ## file name only (no tags)
				$tempOutFileName =~ s/^($outFileTag)//;  ## removing the outFilePrefix
				$tempOutFileName = $tempOutDir.$pathVar.$i.$tempOutFileName;

				##print "\ntempOutFileName: $tempOutFileName\n";

				append_file($tempOutFileName,$bedOutFile,$outDir);  ## inner dir to final out dir
				}  ## for(each thread) ends
			} ## elsif(windows system) ends
			#############################################
		}  ## for(bedOutFilesHashRef) ends
	} ## function ends
	###########################
	sub browser_track_color_hash{
	my %colorHash = (
		"ExonOnly" => "255,153,18",
		"ExonDeletion" => "139,69,0",
		"IntronOnly" => "255,165,79",
		"NeighboringExon" => "205,133,63",
		"MultipleAnnotations" => "139,76,57",
		"GeneDesert" => "132,132,132",
		"InternalExonExtension" => "238,154,73",
		"AlternativeTSS" => "160,82,45",
		"AlternativePolyA" => "205,133,63",
		"Uncharacterized" => "40,40,40",

		## new ##
		"FullTranscriptMatch" => "0,0,204",
		);
	return(%colorHash);
	}  ## function ends
	###########################
	sub generate_track_name{
	my($key,$colorRGB) = @_;
	my $trackLine = "track name=$key description=$key color=$colorRGB visibility=full";
	return($trackLine);
	}  ## function ends
	###########################
	sub get_path_var{

	my $pathVar = "/";  ## default is UNIX style

		if( check_if_windows() ){
		$pathVar ="\\";
		}
	return($pathVar);
	}  ## function ends
	###########################
	sub check_if_windows{
		if($^O =~ m/mswin/i){
		return 1;
		}
		else{
		return 0;
		}
	}  ## function ends
	###########################
	sub get_file_extension{
	my($fileName) = @_;
	my @temp = split(/\./,$fileName);
	my $fileExt = $temp[scalar(@temp)-1];
	return($fileExt);
	}  ## function ends
	###########################

1;



