#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Sort::Naturally;

## Usage definition
my $usage = "\nUSAGE = perl list_CAZy.pl [options]\n
EXAMPLE (simple): list_CAZy.pl -i *CAZ.txt";
my $hint = "Type list_CAZy.pl -h (--help) for list of options\n";
die "$usage\n$hint\n" unless@ARGV;

## Defining options
my $options = <<'END_OPTIONS';
OPTIONS:

-h (--help)	Display this list of options
-i (--input)		Input files [default: *.txt]

END_OPTIONS

my $help ='';
my @files;

GetOptions(
	'h|help' => \$help,
	'i|input=s@{1,}' => \@files
);

if ($help){die "$usage\n$options";}

my $file;
my @array;
while ($file = shift@files){
	open IN, "<$file";
	open OUT, ">cazy.list";

	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(FUN\S+)\s+note\s+CAZy\:(\w+)/){
			my $cazy=$2;
			if ($cazy~~ @array){next;}
			else {
				push @array, "$cazy";
			}
		}
	}
}
my@list=nsort @array;
while (my$name = shift@list){
    print OUT "$name\n";
}