#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## Usage definition
my $usage = "\nUSAGE = perl count_CAZy.pl [options]\n
EXAMPLE (simple): count_CAZy.pl -i CAZ.txt(from funannotate)";
my $hint = "Type count_CAZy.pl -h (--help) for list of options\n";
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
open OUT, ">CAZy_result.out";
while ($file = shift@files){
	my $AA=0;
	my $CBM=0;
	my $CE=0;
	my $GH=0;
	my $GT=0;
	my $PL=0;
	open IN, "<$file";
	for (<IN>) {
		$AA++ if /note	CAZy:AA/;
		$CBM++ if /note	CAZy:CBM/;
		$CE++ if /note	CAZy:CE/;
		$GH++ if /note	CAZy:GH/;
		$GT++ if /note	CAZy:GT/;
		$PL++ if /note	CAZy:PL/;
	}
	print OUT "$file\t$AA\t$CBM\t$CE\t$GH\t$GT\t$PL\n";
}