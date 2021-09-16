#!/usr/bin/perl
use strict;
use warnings;

# Last edited 210827. Removed excess variables and cleaned file.
# Input: 2 line fasta format (>seqID\nseq\n)
# Output: Unique alleles with seqIDs concatenated as header as (>seqIDs\nseq\n)
# IF USING CLUSTAL may need to allele map these guys. Will cut header length.
# To use allele map, see script ""

# Declare variables
my ($query, $row, $key, $seqID);
my %uniqseqs = ();
my $outfile = "uniqseqs.fasta";

# Verify input
if (!defined $ARGV[0]) {
    &print_usage("\nPlease specify a fasta file.");
}

# Open input file
$query = $ARGV[0];
if (-f $query) {
    open F, '<', $query; 
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }   
            if ($row =~ s/^>//) {
                $seqID = $row;
                $seqID =~ s/^\>//d;
            }   
            else {
                $key = $row;
                $uniqseqs{$key} .= "$seqID ";
            }
        }
    close F;
} 
else {
    print"\n Error in read file\n";
    exit 1 
}

# Generate output
open F, '>', $outfile or die "problem saving output to file\n";

foreach $key (keys %uniqseqs) {
    print F ">$uniqseqs{$key}\n";
    print F "$key\n";
}
close F;

print "\nOutfile is $outfile\n\n";
print "If using outfile for clustal, make sure you allele mapped first.\n";

##########################################
sub print_usage {
    my ($error) = @_; 
    if (defined $error) {
        print STDERR $error, "\n";
    }   

    print "\nUsage: $0 [fasta file]\n"; 
    print "\tInput: 2 line fasta format (>seqID]\\nseq\\n)\n";
    print "\tOutput: Unique alleles with seqIDs concatenated as header as (>seqIDs\\nseq\\n)\n";
    print "\tIf using output for CLUSTAL, it's advised to allele map these seqs first\n";
    print "\nCheers!\n\n";
    exit 1; 
}

