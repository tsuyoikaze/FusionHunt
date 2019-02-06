package shared::FileCheckPrint;


	##sub check_if_readable{
	##sub print_cluster_to_file{
	##sub print_hash_to_file{

	########### Checks if filenames in the array passed to the fucntion exists ###########
	########### Error will terminate the program ###########
	sub check_if_readable{
	my(@fn) = @_;
	
		for(my $i = 0;$i<scalar(@fn);++$i){

		##print "\n input file name: $fn[$i] \n";

			##if(sysopen(FHCheckFileExist,$fn[$i],O_RDONLY)){
			##print "\n Yes, the file was opened \n";
			##}
			##else{
			##die "\n can not open input file $fn[$i] \n\n";
			##}
	
		
		open FHCheckFileExist,$fn[$i] or die "\n can not open input file $fn[$i] \n\n";
		close FHCheckFileExist;
		}
	splice(@fn);	
	}  ## function ends
	###########################
	sub print_cluster_to_file{
	my($fout,$arrayLen,$n,$clustRef) = @_;
	##print "\n printing cluster to file \n";
	open FHPrintCluster,">>$fout";  ## outfile (opened in append mode)
		for(my $i = 0;$i<$n;++$i){
		print FHPrintCluster join("\t",@{$clustRef->[$i]}[0..($arrayLen-1)]),"\n";
		}
	close FHPrintCluster;
	}  ## function ends
	###########################
	sub print_hash_to_file{
	my($hashRef,$outFileName) = @_;
	open FHPrintHashToFile,">>$outFileName";  ## outfile
	
		for my $key1(sort(keys(%{$hashRef}))){
		print FHPrintHashToFile "$key1\t$hashRef->{$key1}\n";
		}
	close FHPrintHashToFile;
	}  ## function ends
	###########################
	sub check_if_dir_accessible{
	my($dir) = @_;
	$dir =~ s/\s//g;
	opendir MD,$dir or die "\n can not open directory $dir \n";
	closedir MD;
	}  ## function ends
	###########################
1;

