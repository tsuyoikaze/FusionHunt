package shared::CheckStringsNumericals;

require Exporter;

@ISA = qw/Exporter/;
@EXPORT = qw/check_if_numerical/;

use strict;
use warnings;

	sub check_if_numerical{  ## checks if given string is completely numerical
	my($str1) = @_;
	my($str2,$numFlag);
	$numFlag = 1;  ## by default numeric

	$str1 =~ s/^(\s+)//;
	$str1 =~ s/(\s+)$//;

		if($str1 eq ""){  ## if $str1 is non-blank
		$numFlag = 0;
		}
		else{
		$str1 =~ s/\d//g;
			if($str1 ne ""){  ## if string is not null
			$numFlag = 0;
			}
		}
	return($numFlag);
	}  ## function ends
	###########################
	sub check_if_num{  ##4
	my($str) = @_;
	my $NumFlag  = 1;
	$str =~ s/(\.)//;
	$str =~ s/\d//g;
		if($str){
		$NumFlag = 0;
		}
	return($NumFlag);
	}  ## function ends
	###########################
	sub check_if_integer{
	my($value) = @_;
	my $dataTypeCheck = 1;  ## true by default
	$value =~ s/^(\s+)//;
	$value =~ s/(\s+)$//;
	$value =~ s/\d//g;
		if($value){
		$dataTypeCheck = 0;
		}
	return($dataTypeCheck);
	}  ## function ends
	###########################
	sub check_if_string{
	my($value) = @_;
	my $dataTypeCheck = 1;  ## true by default
	$value =~ s/^(\s+)//;
	$value =~ s/(\s+)$//;
		if(!($value)){
		$dataTypeCheck = 0;
		}
	return($dataTypeCheck);
	}  ## function ends
	###########################
	sub check_if_float{
	my($value) = @_;
	my $dataTypeCheck = 1;  ## true by default
	$value =~ s/^(\s+)//;
	$value =~ s/(\s+)$//;
	$value =~ s/\d//g;

		if( ($value) && ($value !~ '.') ){
		$dataTypeCheck = 0;
		}
	return($dataTypeCheck);
	}  ## function ends
	###########################
	########### Round off a given number ###########
	sub round{
	my($number) = @_;
	return int($number + .5 * ($number <=> 0));
	}  ## function ends
	###########################
	sub order_of_magnitude{
	my($num) = @_;
	my $order = 0;
		if($num > 1){
		$num = int($num);
			while($num){
			$num = int($num/10);
				if($num){
				++$order;
				}
			}  ## while($num) ends
		}  ## if(num > 1) ends
		elsif($num < 1){
			while($num < 1){
			++$order;
			$num = $num*10;
			}   ## while($num < 1) ends
		$order *= (-1);
		}  ## elsif(number < 1) ends
	return($order);
	}  ## function ends
	###########################
	sub get_percentage{
	my($numerator,$denominator) = @_;
		if($denominator){  ## if denominator is non zero
		return(($numerator/$denominator)*100);
		}
		else{
		return(0);
		}
	}  ## function ends
	###########################

1;

