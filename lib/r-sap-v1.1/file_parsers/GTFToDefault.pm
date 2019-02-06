package file_parsers::GTFToDefault;

use strict;
use warnings;

use file_parsers::CheckIfGTF;

	sub gtf_to_default{
	my($geneFile,$outFile) = @_;

	open FHGTFToDefaultIn,$geneFile or die "\n can not open file $geneFile \n";
	open FHGTFToDefaultOut,">$outFile";  ## out file

	my($str,@a1,@temp);
	my(@transcriptInfo,$n,$l,$i);
	my(@exSt,@exEn,@cdsSt,@cdsEn);
	my($orfSt,$orfEn);
	my(@startSt,@startEn);
	my(@stopSt,@stopEn);
	my(@defaultArray);

	my($pt1,$pt2);

	$n = 0;
	$pt1 = tell FHGTFToDefaultIn;

		while($str = <FHGTFToDefaultIn>){
		$pt2 = tell FHGTFToDefaultIn;
		$str =~ s/\n//;
		$str =~ s/\r//;

			if( ($str =~ m/^(\#+)/) || !($str) ){  ## skipping the header blank line
				if(!($n)){
				next;	
				}
				else{
				@a1 = @{$transcriptInfo[$n-1]}[0..($l-1)];
				}
				$pt1 = $pt2;
			}  ## if ends

		@a1 = split(/\t+/,$str);

			if(file_parsers::CheckIfGTF::check_if_GTF(\@a1)){  ## if record is in GTF format
			@temp = separate_goup_elements($a1[8]);  ## returns array containing transcriptID, geneID and geneName etc... in order
			@a1 = @a1[0..7];  ## slice from the original GTF record array
			push(@a1,@temp);
			$l = scalar(@a1);  ## length of the gff recored (tab separated)
			splice(@temp);

				if( !($n) || ($transcriptInfo[$n-1][8] eq $a1[8]) ){  ## if transcript ID is equal
				push(@transcriptInfo,[@a1]);
				++$n;
				$pt1 = $pt2;
				}
			}  ## if(check_if(GTF)) ends
			else{	
			die "\n file format GTF not recognized at line $. in gene track file $geneFile \n";
			}
		
			if( (eof) || ($transcriptInfo[$n-1][8] ne $a1[8]) ){
			sort_gff_records($l,\@transcriptInfo);  ## sorting for features and start position

			get_feature_coordinates(\@transcriptInfo,\@exSt,\@exEn,"exon");  ## exon coordinates
			get_feature_coordinates(\@transcriptInfo,\@cdsSt,\@cdsEn,"cds");  ## CDS coortinates
			get_feature_coordinates(\@transcriptInfo,\@startSt,\@startEn,"start");  ## start codon coordinates
			get_feature_coordinates(\@transcriptInfo,\@stopSt,\@stopEn,"stop");  ## stop codon coordinates

				if($transcriptInfo[0][6] eq '-'){  ## if strand is negative
				interchange_arrays(\@startSt,\@stopSt);
				interchange_arrays(\@startEn,\@stopEn);
				}

				if(scalar(@startSt)){  ## if start codon coordinates found
				$orfSt = $startSt[0];
				}
				elsif(scalar(@cdsSt)){  ## if there is any coding region in the transcript
				$orfSt = $cdsSt[0];
					if($transcriptInfo[0][6] eq '-'){
					$orfSt -= 3;  ## subtracting the 3 bases of the stop codon
					}
				}
				else{
				$orfSt = $exSt[0];
				}

				if(scalar(@stopSt)){
				$orfEn = $stopEn[scalar(@stopEn)-1];
				}
				elsif(scalar(@cdsSt)){
				$orfEn = $cdsEn[scalar(@cdsEn)-1];
					if($transcriptInfo[0][6] eq '+'){
					$orfEn += 3;  ## adding the 3 bases of the stop codon
					}
				}
				else{
				$orfEn = $exSt[0];
				}
			push(@defaultArray,$transcriptInfo[0][8]);  ## transcriptID
			push(@defaultArray,$transcriptInfo[0][0]);  ## chromosome
			push(@defaultArray,$transcriptInfo[0][6]);  ## strand
			push(@defaultArray,($exSt[0],$exEn[scalar(@exEn)-1]));  ## txSt, txEn
			push(@defaultArray,($orfSt,$orfEn));  ## orfSt, orfEn
			push(@defaultArray,scalar(@exSt));  ## number of exons
			push(@defaultArray,join(',',@exSt));  ## exon starts
			push(@defaultArray,join(',',@exEn));  ## exon ends
			push(@defaultArray,$transcriptInfo[0][9]);  ## gene name or gene id  (if any) or '.' otherwise

			print FHGTFToDefaultOut join("\t",@defaultArray),"\n";
		
			shared::DeleteDataStructures::undefine_arrays(\@transcriptInfo,\@exSt,\@exEn,\@cdsSt,\@cdsEn,
				\@startSt,\@startEn,\@stopSt,\@stopEn,\@defaultArray);
			$n = 0;
			seek(FHGTFToDefaultIn,$pt1,0);  ## resetting the file pointer
			}  ## if($id1 ne $id2)
		}  ## while(<FHGTFToDefaultIn>) ends
	close FHGTFToDefaultIn;
	close FHGTFToDefaultOut;
	}  ## function ends
	###########################
	sub separate_goup_elements{
	my($str) = @_;
	my(@temp1,@temp2,@returnArray,%groupHash,$i);
	$str =~ s/^(\s+)//;
	$str =~ s/(\s+)$//;
	@temp1 = split(';',$str);  ## grupp elements are always ';' separated
		for($i = 0;$i<scalar(@temp1);++$i){
		$temp1[$i] =~ s/^(\s+)//;
		$temp1[$i] =~ s/(\s+)$//;
		$temp1[$i] =~ s/\"//g;
		@temp2 = split(/\s+/,$temp1[$i]);
		$groupHash{$temp2[0]} = $temp2[1];
		}
	push(@returnArray,$groupHash{"transcript_id"});
		if(exists $groupHash{"gene_name"}){
		push(@returnArray,$groupHash{"gene_name"});
		}
		elsif(exists $groupHash{"gene_id"}){
		push(@returnArray,$groupHash{"gene_id"});
		}
		else{
		push(@returnArray,".");
		}
	##print "\nreturnArray: ",join("\t",@returnArray),"\n";
	return(@returnArray);  ## returning transcript id and gene name or gene id (if provided in the GTF file)
	}  ## function ends
	###########################
	sub sort_gff_records{
	my($l,$clustRef) = @_;
	my($i,$j,@temp);	

		for($i = 0;$i<(scalar(@{$clustRef})-1);++$i){
			for($j = 0;$j<(scalar(@{$clustRef})-1-$i);++$j){

				if($clustRef->[$j+1][2] lt $clustRef->[$j][2]){  ## comparing features
				shared::ArrayOperations::swap_elements($j,$j+1,$l,$clustRef);
				}  ## if(ends)
				elsif($clustRef->[$j+1][2] eq $clustRef->[$j][2]){  ## comparing chromosomal start positions
					if($clustRef->[$j+1][3] < $clustRef->[$j][3]){
					shared::ArrayOperations::swap_elements($j,$j+1,$l,$clustRef);
					}
				}  ## elsif() ends
			}  ## for($j) ends
		}  ## for($i) ends
	}  ## function ends
	###########################
	sub get_feature_coordinates{
	my($clustRef,$stRef,$enRef,$feature) = @_;
	my($i);
		for($i = 0;$i<scalar(@{$clustRef});++$i){
			if($clustRef->[$i][2] =~ m/$feature/i){  ## if record belongs to the $feature (exon or CDS)
			push(@{$stRef},($clustRef->[$i][3]-1));  ## converting 1 start coordinate to 0 start coordinate
			push(@{$enRef},$clustRef->[$i][4]);
			}
		}  ## for(scalar(@{$clustRef})) ends
	}  ## function ends
	###########################
	sub interchange_arrays{
	my($a1Ref,$a2Ref) = @_;
	my @temp;
	@temp = @{$a1Ref};
	@{$a1Ref} = @{$a2Ref};
	@{$a2Ref} = @temp;
	##my $temp;
	##$temp = $a1Ref;
	##$a1Ref = $a2Ref;
	##$a2Ref = $temp;
	}  ## function ends
	###########################
1;

## sort cluster (first on feature then on start positions);
## print in a temp file to check the sorting
	## get exonSt and exonEn coordinates
	## get cdsSt and cdsEn coordinates
	## define txSt, txEn, exSts, exEns
	## if no CDS,
		##cdsSt: $exSt (subtract -1 for 0 indexing)
		##cdsEn: $exSt
	## if CDS
		## if strand +
		##cdsSt = cdsSt
		##cdsEn = cdsEn+3
		## if strand -
		##cdsSt = cdsSt-3
		##cdsEn = cdsEn
	## subtract -1 from all starts (0 based indexing)
		## Convert the GTF record
## print converted record to out file
## start over for the new  record
## move the file pointer back
		
	


