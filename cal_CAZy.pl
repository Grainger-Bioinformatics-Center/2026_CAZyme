#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Sort::Naturally;

## Usage definition
my $usage = "\nUSAGE = perl cal_CAZy.pl [options]\n
EXAMPLE (simple): cal_CAZy.pl -i *CAZ.txt -l cazy.list";
my $hint = "Type cal_CAZy.pl -h (--help) for list of options\n";
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

## make a hash of all the cazy from the list
open OUT, ">>all_cazy.tsv";
my $list,
my %hash;
while ($list = shift@list){
  open IN1, "<$list";
 	while (my $line1 = <IN1>){
		chomp $line1;
	  $hash{$line1}=0;
 }
}

## print out the header
print OUT "\t";
foreach my $key (nsort keys %hash) {print OUT "$key" . "\t";}
print OUT "\n";

my $file;
while ($file = shift@files){
	open IN, "<$file";
  print OUT "$file\t";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(FUN\S+)\s+note\s+CAZy\:(\w+)/){
			my $cazy=$2;
      $hash{$cazy} ++ if exists $hash{$cazy};
		}
	}
  foreach my $key (nsort keys %hash) {print OUT "$hash{$key}" . "\t";}
  print OUT "\n";
  foreach my $key (nsort keys %hash) {$hash{$key}=0;} ##clear the count for each key
}