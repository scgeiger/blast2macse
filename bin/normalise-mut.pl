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
my ($input, $outfile, $filename); 
my ($ntPos, $ntPosN, $length, $ntMutID, $aaPos, $aaPosN, $aaMutID, $codIn);
my ($temp, $row, $geneID);
my @inFile = ();
my $header = "#geneID\tntPosN\tntMutID\taaPosN\taaMutID\tcodonIndex\n";

# Verify input
$input = $ARGV[0];

($filename) = ($input =~ /([^.]+)/);
$outfile = $filename . "-normalised-pos-muts.txt";
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
        if ( -z $outfile ) {
            print O $header;
        }
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/ || $row =~ /^pos/ ) {
                next;
            }
            else {
                @inFile = split "\t", $row;
                $length = $length = $inFile[11];
                until ($temp = $length % 3) { # we need length make sure std
                    $length++;
                }

                $ntPos =  $inFile[0];
                $ntPosN = ($ntPos / $length);
                $ntMutID = $inFile[3];
                $aaPos = $inFile[5];
                $aaPosN = ($aaPos / ($length/3));
                $aaMutID = $inFile[10];
                $codIn = $inFile[4];
                print O "$geneID\t$ntPosN\t$ntMutID\t$aaPosN\t$aaMutID\t$codIn\n";
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
    print "\tSuffix to uncleaned file is: *geneID.macse-codon-dist.tsv\n";
    print "\tSuffix to cleaned file is: *geneID-removed.macse-codon-dist.tsv\n";
    print "\nCheers!\n\n";
    exit 1;
}
