package Control::TimeModule;

use strict;
use warnings;


	
our %timeIndex = (
	"hour" => 0,
	"minute" => 1,
	"second" => 2,
	"month" => 3,
	"day" => 4,
	"year" => 5,
	);

	sub get_current_time{
	my @time = current_time();
	return(@time);
	}  ## function ends
	###########################
	sub current_time{
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$yearOffset += 1900;
	my @currTime = ($hour,$minute,$second,$month,$dayOfMonth,$yearOffset);
	return(@currTime);
	}  ## function ends
	###########################
	sub time_duration{
	my($time1Ref,$time2Ref) = @_;
	my($days1,$days2,$hours1,$hours2,$difference);

	$difference = ( ($time2Ref->[$timeIndex{"year"}]*365 - $time1Ref->[$timeIndex{"year"}]*365) + 
			($time2Ref->[$timeIndex{"month"}]*30 - $time1Ref->[$timeIndex{"month"}]*30) + 
			($time2Ref->[$timeIndex{"day"}] - $time1Ref->[$timeIndex{"day"}]) )*24*60*60; ## in seconds

	$difference += ($time2Ref->[$timeIndex{"hour"}]*60*60 - $time1Ref->[$timeIndex{"hour"}]*60*60) + 
			($time2Ref->[$timeIndex{"minute"}]*60 - $time1Ref->[$timeIndex{"minute"}]*60) + 
			($time2Ref->[$timeIndex{"second"}] - $time1Ref->[$timeIndex{"second"}]);  ## in seconds

	return($difference);  ## returns differences in minutes
	}  ## function ends
	###########################
	
1;