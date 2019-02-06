package compartment_D::ChimerGeneFindSR;

use strict;
use warnings;

use compartment_D::FindGenicRegions;

	sub generate_genic_region_indexes{
	my %genicRegionIndexes = (  ## genic regions are assigned priorities
		"CDS" => 6, 
		"5UTR" => 5,
		"3UTR" => 4,
		"NonCoding" => 3,
		"INTRON" => 2,
		"intergenic" => 1,
		);
	return(%genicRegionIndexes);
	}  ## function ends
	###########################
	sub get_all_intersected_tracks{
	my($chr,$algnPoint,$knGnIndxRef,$gnIndxRef,$finRef,$clustRef) = @_;
	my ($indx,$str,@geneArray);
	my($i);

		if(exists $gnIndxRef->{$chr}){
		$indx = $gnIndxRef->{$chr};  ## getting gene track index using chromosome name

		##open FHChGeneFind, $fin or die "\n can not open file $fin \n";
		##seek(FHChGeneFind,$indx,0);  ## setting file pointer to the current chromosome

			##while($str = <FHChGeneFind>){
			##$str =~ s/\n//;
			##$str =~ s/\r//;
			##@geneArray = split(/\t+/,$str);

			for($i = $indx;$i<scalar(@{$finRef});++$i){
			@geneArray = split(/\t+/,$finRef->[$i]);

				if($chr ne $geneArray[$knGnIndxRef->{"chr"}]){  ## if chr doesn't match
				last;
				}
				if(shared::CoordinateComparisons::intersect_track_alignment( ($algnPoint-1),$algnPoint,$geneArray[$knGnIndxRef->{"txSt"}],
								$geneArray[$knGnIndxRef->{"txEn"}])){  ## if there is any intersect
				push(@{$clustRef},[@geneArray]);  ## populating possible gene tracks in the matrix
				}			
			splice(@geneArray);
			}  ## for(@fin) ends


			##}  ## while(<FHChGeneFind>) ends
		##close FHChGeneFind;
		}  ## if(chr exists) ends

		#else{
		#die "\n no gene track found for chromosome $chr \n";  ## error message
		#} 	
	}  ## function ends
	###########################	
	sub find_alignment_region{
	my($algnPoint,$clustRef,$currTrack,$knGnIndxRef) = @_;
	my($regionTag,$flag);
	my(@exSt,@exEn,$cdsSt,$cdsEn);
	my(@regionSts,@regionEns);

	compartment_D::FindGenicRegions::get_cds_pos($clustRef,$currTrack,$knGnIndxRef,\$cdsSt,\$cdsEn);  ## getting coding start and end positions
	compartment_D::FindGenicRegions::get_exon_cords($clustRef,$currTrack,$knGnIndxRef,\@exSt,\@exEn);  ## getting exon cordinates

	compartment_D::FindGenicRegions::get_coding_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,\@regionSts,\@regionEns);
		if(check_if_in_the_region($algnPoint,\@regionSts,\@regionEns)){
		$regionTag = "CDS";
		return($regionTag);
		}
		else{
		splice(@regionSts);
		splice(@regionEns);
		compartment_D::FindGenicRegions::get_noncoding_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,\@regionSts,\@regionEns);
			if(check_if_in_the_region($algnPoint,\@regionSts,\@regionEns)){
			$regionTag = "NonCoding";
			return($regionTag);
			}
			else{
			splice(@regionSts);
			splice(@regionEns);
			compartment_D::FindGenicRegions::get_5UTR_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,
					$clustRef->[$currTrack][$knGnIndxRef->{"strand"}],\@regionSts,\@regionEns);
				if(check_if_in_the_region($algnPoint,\@regionSts,\@regionEns)){
				return("5UTR");
				}
				else{
				splice(@regionSts);
				splice(@regionEns);
				compartment_D::FindGenicRegions::get_3UTR_regions(\@exSt,\@exEn,$cdsSt,$cdsEn,
					$clustRef->[$currTrack][$knGnIndxRef->{"strand"}],\@regionSts,\@regionEns);
					if(check_if_in_the_region($algnPoint,\@regionSts,\@regionEns)){
					return("3UTR");
					}
					else{
					splice(@regionSts);
					splice(@regionEns);
					compartment_D::FindGenicRegions::get_intron_regions(\@exSt,\@exEn,\@regionSts,\@regionEns);
						if(check_if_in_the_region($algnPoint,\@regionSts,\@regionEns)){
						return("INTRON");
						}
					}  ## else ends
				}  ## else ends
			}  ## else ends
		}  ## all else blocks end
	}  ## function ends
	###########################
	sub check_if_in_the_region{
	my($point,$stRef,$enRef) = @_;
	my($i,$flag);
	$flag = 0;
		for($i = 0;$i<scalar(@{$stRef});++$i){
			if(shared::CoordinateComparisons::intersect_track_alignment( ($point-1),$point,$stRef->[$i],$enRef->[$i])){
			$flag = 1;
			last;
			}
		}  ## for($i) ends
	return($flag);
	}  ## function ends
	###########################
1;
