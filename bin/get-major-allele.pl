#!/usr/bin/perl
use strict;
use warnings;


my $majority = 0;
my $row;
my $seqid;
my %all_seqs = ();
my @ids = ();
my $major_id;
my $major_seq;
my $infile = $ARGV[0];
$infile =~ s/\..*//;
my $outfile = $infile . ".majority-allele.fasta";
my $total_seqs = 0;
my $major_freq;

# note that, becasue of that, sequences that are all gaps won't count
# input file is *.dealign
if (-f $ARGV[0]) {
    open F, '<', $ARGV[0];
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            if ($row =~ /^>/) {
                $row =~ s/^>//;
                $seqid = $row;
            }
            else {
                $all_seqs{$seqid} = $row; 
            }
        }
    close F;
}

foreach $seqid (keys %all_seqs) {
    @ids = split / /, $seqid;
    $total_seqs = $total_seqs + scalar(@ids);
    if ($majority < scalar(@ids)) {
        $majority = scalar(@ids);
        $major_id = $seqid;
        $major_seq = $all_seqs{$seqid};
        $major_freq = scalar(@ids);
    }
}

open O, '>', $outfile;
    print O "#input: $ARGV[0]; total seqs: $total_seqs; major freq: $major_freq\n";
    print O ">$major_id\n";
    print O "$major_seq";
close O;
