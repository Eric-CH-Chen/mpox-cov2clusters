#!/usr/bin/env perl
use strict;
use warnings;

# Renames sequences output from FAVITES_Lite to contain only
# the node (individual) and no other information

# Example: 
#	Input:
#	>170|17876|714.320312
#	Output:
#	>17876
my $file = shift;
my $o_file = shift;

open (FILE, "<", $file) or die;
open (OFILE, ">", $o_file);
while(<FILE>) {
	my $line = $_;
	chomp $line;

	if ($line =~ />/) {
		my @split_line = split(/\|/, $line);

		my $new_header = '>' . $split_line[1];
		print OFILE $new_header . "\n";
		
	} else {
		print OFILE $line . "\n";
		next;
	}
}

close FILE;
close OFILE;
