#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## Usage definition
my $usage = "\nUSAGE = perl rename-cazy-protein.pl [options]\n
EXAMPLE (simple): rename-cazy-protein.pl -i *.fa";
my $hint = "Type rename-cazy-protein.pl -h (--help) for list of options\n";
die "$usage\n$hint\n" unless@ARGV;

## Defining options
my $options = <<'END_OPTIONS';
OPTIONS:

-h (--help)    Display this list of options
-i (--input)    Input files [default: *.fa]


END_OPTIONS

my $help ='';
my @files;
my @list;

GetOptions(
	'h|help' => \$help,
	'i|input=s@{1,}' => \@files
	
);

if ($help){die "$usage\n$options";}


## split and rename
my $file;
while ($file = shift@files){
	open IN, "<$file";
	(my $without_extension = $file) =~ s/\.fa$//;
	open OUT, ">$without_extension.fasta";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^\>(FUN\S+)/){
			my $gene=$1;
			print OUT '>'."$without_extension".'-'."$gene\n";
		}
		else {print OUT "$line\n";}
	}
}