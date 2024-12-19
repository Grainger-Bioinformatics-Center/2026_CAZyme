#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Sort::Naturally;

# Usage definition
my $usage = <<'USAGE';

USAGE:
    perl cal_CAZy.pl [options]#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Sort::Naturally;

# Usage and help messages
my $usage = <<'USAGE';
USAGE:
    perl cal_CAZy.pl [options]

EXAMPLE:
    perl cal_CAZy.pl -i *CAZ.txt -l cazy.list -o all_cazy.tsv

DESCRIPTION:
    This script calculates the frequency of CAZyme annotations from input files
    and writes the results in a tab-delimited output file.
USAGE

my $options = <<'OPTIONS';
OPTIONS:
    -h (--help)    Display this help message
    -i (--input)   Input files (space-separated, supports multiple files)
    -l (--list)    Existing CAZyme list file (default: cazy.list)
    -o (--out)     Output file name (default: all_cazy.tsv)
OPTIONS

# Variables for command-line options
my $help = 0;
my @input_files;
my $list_file = 'cazy.list';
my $out_file = 'all_cazy.tsv';

# Parse command-line options
GetOptions(
    'h|help'    => \$help,
    'i|input=s{1,}' => \@input_files,
    'l|list=s' => \$list_file,
    'o|out=s'  => \$out_file
) or die "Error in command-line arguments. Use -h for help.\n";

# Display help and exit if requested
if ($help) {
    print "$usage\n$options";
    exit;
}

# Check for required input files
if (!@input_files) {
    die "Error: No input files specified. Use -i option to provide input files.\n";
}

# Open output file
open my $out_fh, '>', $out_file or die "Could not open output file '$out_file': $!\n";

# Load CAZyme list into a hash
my %cazy_counts;
open my $list_fh, '<', $list_file or die "Could not open list file '$list_file': $!\n";
while (my $line = <$list_fh>) {
    chomp $line;
    $cazy_counts{$line} = 0;
}
close $list_fh;

# Print header in output file
print $out_fh "\t", join("\t", nsort keys %cazy_counts), "\n";

# Process each input file
foreach my $file (@input_files) {
    # Skip if file does not exist
    unless (-e $file) {
        warn "Warning: Input file '$file' does not exist. Skipping...\n";
        next;
    }

    open my $in_fh, '<', $file or die "Could not open input file '$file': $!\n";

    # Print file name as the first column
    print $out_fh "$file\t";

    # Reset counts for each file
    my %file_counts = %cazy_counts;

    while (my $line = <$in_fh>) {
        chomp $line;

        # Match lines containing FUNxxxx and CAZy:xxxx pattern
        if ($line =~ /^(FUN\S+)\s+note\s+CAZy:(\w+)/) {
            my $cazy = $2;
            $file_counts{$cazy}++ if exists $file_counts{$cazy};
        }
    }

    # Output counts for the current file
    print $out_fh join("\t", map { $file_counts{$_} } nsort keys %file_counts), "\n";

    close $in_fh;
}

close $out_fh;

# Completion message
print "Processing complete. Results written to '$out_file'.\n";


EXAMPLE (simple):
    cal_CAZy.pl -i *CAZ.txt -l cazy.list -o all_cazy.tsv

USAGE

my $hint = "Type cal_CAZy.pl -h (--help) for list of options\n";

# If no arguments are provided, display usage
die "$usage\n$hint\n" unless @ARGV;

# Defining options
my $options = <<'OPTIONS';
OPTIONS:

    -h (--help)    Display this list of options
    -i (--input)   Input files [default: *.txt]
    -l (--list)    Existing CAZyme list [default: cazy.list]
    -o (--out)     Output file name [default: all_cazy.tsv]

OPTIONS

# Initialize variables for command-line options
my $help = 0;
my @files;
my $list_file = 'cazy.list';  # Default value for the list
my $out_file = 'all_cazy.tsv';  # Default output file

GetOptions(
    'h|help' => \$help,
    'i|input=s{1,}' => \@files,
    'l|list=s' => \$list_file,
    'o|out=s' => \$out_file
) or die "Error in command-line arguments\n";

# Display help message if the help option is triggered
if ($help) {
    die "$usage\n$options";
}

# Open output file for writing
open my $out_fh, '>', $out_file or die "Could not open output file '$out_file': $!\n";

# Load CAZyme list into a hash with initial count of 0 for each item
my %hash;
open my $list_fh, '<', $list_file or die "Could not open list file '$list_file': $!\n";
while (my $line = <$list_fh>) {
    chomp $line;
    $hash{$line} = 0;
}
close $list_fh;

# Print header in output file
print $out_fh "\t", join("\t", nsort keys %hash), "\n";

# Process each input file
foreach my $file (@files) {
    open my $in_fh, '<', $file or die "Could not open input file '$file': $!\n";
    
    # Print the file name as the first column
    print $out_fh "$file\t";
    
    # Reset the CAZyme counts before processing each file
    my %file_hash = %hash;  # Copy original hash with counts reset to 0
    
    while (my $line = <$in_fh>) {
        chomp $line;
        
        # Match lines with FUNxxxx and CAZy:xxxx pattern
        if ($line =~ /^(FUN\S+)\s+note\s+CAZy:(\w+)/) {
            my $cazy = $2;
            $file_hash{$cazy}++ if exists $file_hash{$cazy};
        }
    }
    
    # Output the counts for each CAZyme in the file
    print $out_fh join("\t", map { $file_hash{$_} } nsort keys %file_hash), "\n";
    
    close $in_fh;
}

close $out_fh;

print "Processing complete. Results written to '$out_file'.\n";
