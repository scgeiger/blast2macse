#!/usr/bin/perl
use strict;
use warnings;

# This script is used for a uniqseqs file, in which the first line of a sequence is the sequence id(s) for the next line, containing the sequence. It's assumed that the second line is the entire sequence, and it's not split across multiple lines.
# Output for this glue script is 1. an allele map, stating which allele holds which ids, and 2. a replacement of the input file with the seqids replaced by allele no. 

my $infile = $ARGV[0];
(my $filename) = ($infile =~ /([^.]+)/);
my $outmap_file = "allele-map-" . $filename . ".txt";
my $outseq_file = "am-" . $filename . ".fasta";
my $allele_index = 0;
my $row;
my @allele_sequence = ();
my @allele_seqids = ();
my @ids = ();
my $i;
my $k;
my $seq;

# Verify input
if (!defined $ARGV[0]) {
    die "\n\t Please specify a uniqseq file to be allele-mapped and try again\n";
}

if (-f $infile) {
    open F, '<', $infile;
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            if ($row =~ /^>/) {
                $row =~ s/^>//;
                $allele_index++;
                $allele_seqids[$allele_index] = $row;
            }
            else {
                $seq = $row;
                $allele_sequence[$allele_index] = $seq; 
            }
        }
    close F;
}

open O, '>', $outmap_file;
    print O "#Frequency: ";
    foreach $i (1..$#allele_seqids) {
        @ids = split / /, $allele_seqids[$i];
        print O "$i:" . scalar(@ids) . "; ";
    }
    print O "\n";
    foreach $i (1..$#allele_seqids) {
        @ids = split / /, $allele_seqids[$i];
        foreach $k (@ids) {
            print O "allele_" . $i . "\t" . $k . "\n";
        }
    }
close O;

open O, '>', $outseq_file;
    foreach $i (1..$#allele_sequence) {
        print O ">allele_" . $i . "\n" . "$allele_sequence[$i]" . "\n";
    }
close O;

