use strict;
use warnings;
package SenFinishReport;

	sub send_mail{
	my $subject = "programLogFile";
	my $home = $ENV{'HOME'};
	$home =~ s/\n//;
	$home =~ s/\r//;
	my $message = $home."/Mail/TextMessage.msg";
	my $email ="vinaykmittal\@gmail.com";
	`mutt -s $subject -a nohup.out $email < $message`;
	}  ## function ends
	###########################
1;
