#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Usage and help message
my $usage = <<'USAGE';
USAGE:
    perl count_CAZy.pl [options]

EXAMPLE:
    perl count_CAZy.pl -i CAZ.txt

DESCRIPTION:
    This script counts the number of CAZy family annotations (AA, CBM, CE, GH, GT, PL)
    in Funannotate output files and writes the results to an output file.
USAGE

my $options = <<'OPTIONS';
OPTIONS:
    -h (--help)    Display this help message
    -i (--input)   Input files (space-separated, supports multiple files)
OPTIONS

# Variables for command-line options
my $help = 0;
my @input_files;

# Parse command-line options
GetOptions(
    'h|help'    => \$help,
    'i|input=s{1,}' => \@input_files
) or die "Error in command-line arguments. Use -h for help.\n";

# Display help and exit if requested
if ($help) {
    print "$usage\n$options";
    exit;
}

# Check that input files are provided
if (!@input_files) {
    die "Error: No input files specified. Use -i option to provide input files.\n";
}

# Open output file
my $output_file = 'CAZy_result.out';
open my $out_fh, '>', $output_file or die "Could not open output file '$output_file': $!\n";

# Print header to output file
print $out_fh "File\tAA\tCBM\tCE\tGH\tGT\tPL\n";

# Process each input file
foreach my $file (@input_files) {
    # Skip if file does not exist
    unless (-e $file) {
        warn "Warning: Input file '$file' does not exist. Skipping...\n";
        next;
    }

    # Initialize counts for CAZy families
    my %counts = (
        'AA'  => 0,
        'CBM' => 0,
        'CE'  => 0,
        'GH'  => 0,
        'GT'  => 0,
        'PL'  => 0
    );

    # Open input file
    open my $in_fh, '<', $file or die "Could not open input file '$file': $!\n";

    # Count occurrences of each CAZy family
    while (my $line = <$in_fh>) {
        $counts{'AA'}++  if $line =~ /note\s+CAZy:AA/;
        $counts{'CBM'}++ if $line =~ /note\s+CAZy:CBM/;
        $counts{'CE'}++  if $line =~ /note\s+CAZy:CE/;
        $counts{'GH'}++  if $line =~ /note\s+CAZy:GH/;
        $counts{'GT'}++  if $line =~ /note\s+CAZy:GT/;
        $counts{'PL'}++  if $line =~ /note\s+CAZy:PL/;
    }

    close $in_fh;

    # Write results to output file
    print $out_fh join("\t", $file, @counts{qw(AA CBM CE GH GT PL)}), "\n";
}

close $out_fh;

# Completion message
print "CAZy family counts have been written to '$output_file'.\n";
