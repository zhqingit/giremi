#!/usr/bin/perl
#
#
#
#
#
#
#
use Getopt::Long;
use strict;
use warnings;

system("./giremi @ARGV");

my $output = "";
GetOptions(
	'o=s'  => \$output,
);

if ($#ARGV>0){
	my $cmd = "R CMD BATCH --no-save --no-restore '--args finput=\"$output\" fout=\"$output.res\"' giremi.r $output.Rout";
	print $cmd,"\n";
	system($cmd);

	print "Analysis DONE!\n";
}

