package shared::CoordinateComparisons;

	##sub intersect_track_alignment{
	##sub calculate_end_to_end_distance{  ## to calculate the distance if fragments are non-ovelapping (separated from each other)
	##sub get_intersecting_coordinates{  ## intersect two fragments and retruns coordinates of overlapping region
	##sub find_block_overlaps{

	sub intersect_track_alignment{
	my($st1,$en1,$st2,$en2) = @_;  ## note:  0,1 based intersection
	my $intersect = 0;
		if( (($st1+1) >= ($st2+1))  && ($en1 <= $en2) ){  ##1st fragment <= 2nd fragment
		##case 1
		$intersect =  $en1 - ($st1+1) +1;
		}		
		##case2i
		elsif( (($st1+1) >= ($st2+1)) && (($st1+1) <= $en2) && ($en1 > $en2) ){  ## right end is going out
		$intersect = $en2 - ($st1+1) + 1;
		}
		#####case3
		elsif( (($st1+1) < ($st2+1)) && ($en1 >= ($st2+1)) && ($en1 <= $en2) ){  ## left end is coming out
		$intersect = $en1 - ($st2+1) + 1;
		}	
		####case4
		elsif( (($st1+1) < ($st2+1)) && ($en1 > $en2) ){  ## 2nd frag < 1st end
		$intersect = $en2  - ($st2+1) + 1;     
		}	
		if($intersect < 0){  ## just checking for negative overlap
		$intersect = 0;
		}
	return($intersect);
	}  ## function ends
	###########################
	sub calculate_end_to_end_distance{  ## to calculate the distance if fragments are non-ovelapping (separated from each other)
	my($st1,$en1,$st2,$en2) = @_;
		if($en1 < $st2){
		return( ($st2-$en1));
		}
		elsif($st1 > $en2){
		return( ($st1-$en2));
		}
		else{
		return(0);
		}
	}  ## function ends
	###########################
	sub get_intersecting_coordinates{  ## intersect two fragments and retruns coordinates of overlapping region
	my($st1,$en1,$st2,$en2) = @_;
	my($ovSt,$ovEn,$intersect);
	$intersect = 0;
	$ovSt = 0;
	$ovEn = 0;
		if( (($st1+1) >= ($st2+1))  && ($en1 <= $en2) ){  ##1st fragment <= 2nd fragment
		##case 1
		$intersect =  $en1 - ($st1+1) +1;
		$ovSt = $st1;
		$ovEn = $en1;
		}		
		##case2i
		elsif( (($st1+1) >= ($st2+1)) && (($st1+1) <= $en2) && ($en1 > $en2) ){  ## right end is going out
		$intersect = $en2 - ($st1+1) + 1;
		$ovSt = $st1;
		$ovEn = $en2;
		}
		#####case3
		elsif( (($st1+1) < ($st2+1)) && ($en1 >= ($st2+1)) && ($en1 <= $en2) ){  ## left end is coming out
		$intersect = $en1 - ($st2+1) + 1;
		$ovSt = $st2;
		$ovEn = $en1;
		}	
		####case4
		elsif( (($st1+1) < ($st2+1)) && ($en1 > $en2) ){  ## 2nd frag < 1st end
		$intersect = $en2  - ($st2+1) + 1;     
		$ovSt = $st2;
		$ovEn = $en2;
		}	

		if($intersect < 0){  ## just checking for negative overlap
		$intersect = 0;
		$ovSt = 0;
		$ovEn = 0;
		}
	return($intersect,$ovSt,$ovEn);
	}  ## function ends
	###########################
	sub find_block_overlaps{
	my($st1Ref,$en1Ref,$st2Ref,$en2Ref) = @_;
	my($overlap,$i,$j);
	$overlap = 0;
		for($i = 0;$i<scalar(@{$st1Ref});++$i){  ## for each aligned block in the pslx record
			for($j = 0;$j<scalar(@{$st2Ref});++$j){  ## for each exon
			$overlap += intersect_track_alignment($st1Ref->[$i],$en1Ref->[$i],
									$st2Ref->[$j],$en2Ref->[$j]);
			}  ## for($j) ends
		}  ## for($i) ends
		if($overlap < 0){
		$overlap = 0;
		}
	return($overlap);
	}  ## function ends
	###########################

1;
