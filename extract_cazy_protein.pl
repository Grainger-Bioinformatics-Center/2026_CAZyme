#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## Usage definition
my $usage = "\nUSAGE = perl extract_cazy_protein.pl [options]\n
EXAMPLE (simple): extract_cazy_protein.pl -i *.txt -l cazy.list";
my $hint = "Type extract_cazy_protein.pl -h (--help) for list of options\n";
die "$usage\n$hint\n" unless@ARGV;

## Defining options
my $options = <<'END_OPTIONS';
OPTIONS:

-h (--help)    Display this list of options
-i (--input)    Input files [default: *.txt]
-l (--list)    exist cayzme list [default: cazy.list]

END_OPTIONS

my $help ='';
my @files;
my @list;

GetOptions(
	'h|help' => \$help,
	'i|input=s@{1,}' => \@files,
	'l|list=s@{1,}' => \@list
	
);

if ($help){die "$usage\n$options";}

## make a hash of all the cazy from the list, and create the folder
my $list,
my %hash;
while ($list = shift@list){
  open IN1, "<$list";
 	while (my $line1 = <IN1>){
		chomp $line1;
	  $hash{$line1}=0;
 }
}
foreach my $key (keys %hash) {system ("mkdir $key");}

## extract proteins to floders
my $file;
while ($file = shift@files){
	open IN, "<$file";
	(my $without_extension = $file) =~ s/\.txt$//;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(FUN\S+)\s+note\s+CAZy\:(\w+)/){
			my $gene=$1;
			my $cazy=$2;
      system ("samtools faidx $without_extension".'.fa '."$gene".' >> '."$cazy".'/'."$cazy".'_'."$without_extension".'.fa');
		}
	}
}