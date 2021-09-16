#!/usr/bin/perl
use strict;
use warnings;

# Use this to collect the normalised position of each codon mutation
# Input will need to be a path; want to keep the output file in a central place
# Will append to existing file
# Depends on codon distribution file ()
# Will only give you nt positions, not normalised to codon
# Updated 210916

# Declare vars
my ($input, $outfile); 
my ($normPos, $pos, $length, $mutType);
my ($temp, $row, $geneID);
my @inFile = ();

# Verify input
$input = $ARGV[0];

$outfile = "normalised-pos-ntmuts.txt";
if (!defined $ARGV[0] || !defined $input) {
    &print_usage("\nPlease specify a codon distribution file.");
}

# Build in utility to add header if file is empty

# Execute
if ( -f $input ) { 
    $geneID = $input;
    $geneID =~ s/^[^-]*-//;
    $geneID =~ s/-.*//;

    open F, '<', $input;
    open O, '>>', $outfile;
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/ || $row =~ /^pos/ ) {
                next;
            }
            else {
                @inFile = split "\t", $row;
                $pos = $inFile[0];
                $length = $inFile[11];
                $mutType = $inFile[3];
                until ($temp = $length % 3) { # we need length make sure std
                    $length++;
                } 
                $normPos = ($pos/($length));
                print O "$geneID\t$mutType\t$normPos\n";
            } 
        }
    close O;
    close F;
}
else {
    print "$input does not exist\n";
}

# Print Usage
sub print_usage {
    my ($error) = @_;
    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 [codon-dist.tsv]\n";
    print "\tInput: path to *codon-dist.tsv files)\n";
    print "\tOutput: a table that has normalised positions for all mutations; note that >>\n";
    print "\tSuffix to uncleaned file is: *geneID.macse-codon-dist.tsv\n"
    print "\tSuffix to cleaned file is: *geneID-removed.macse-codon-dist.tsv\n"
    print "\nCheers!\n\n";
    exit 1;
}
