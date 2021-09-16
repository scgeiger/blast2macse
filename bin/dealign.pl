#!/usr/bin/perl
use strict;
use warnings;

# Last edited 210827. Cleaned and added exit 1
# Assumes input is fasta format with >seqid\nseq\n

# Declare variables
my ($row, $header, $seqid, $dealign, $i);
my ($infile, $outfile, $errfile);
my ($checkIn, $checkOut, $checkErr) = (0) x 3;
my @split = ();
my %all_seqs = ();

# Verify input
if (!defined $ARGV[0]) {
    &print_usage("\nPlease specify a sequence file.");
}

$infile = $ARGV[0];
$outfile = $ARGV[0] . ".dealign";
$errfile = $ARGV[0] . ".e-dealign";

# Read input
if (-f $ARGV[0]) {
    open F, '<', $ARGV[0];
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                $header = $row;
            }
            if ($row =~ /^>/) {
                $seqid = $row;
                $checkIn++;
            }
            else {
                $all_seqs{$seqid} = $row;
            }
        }
    close F;
}

# Processing and output
open ERR, '>', $errfile;
open O, '>', $outfile;
    if ($header) {
        print O "$header\n";
    }
    foreach $seqid (keys %all_seqs) {
        $dealign = $all_seqs{$seqid};
        $dealign =~ tr/\-//d;
        $dealign =~ tr/\!//d;
        if (length($dealign) == 0) { #sequence $seqid only has gaps, needs removal
            $checkErr++;
            $seqid =~ s/^>//;
            @split = split / /, $seqid;
            foreach $i (@split) { 
                print ERR "$i\n";
            }
            next;
        }
        else {
            print O "$seqid\n$dealign\n";
            $checkOut++;
        }
    }
close O;
close ERR;

# Clean empty error file
if (-z $errfile) {
    unlink($errfile) or die "\nProblem cleaning empty error file\n";  
}
else {
    print "\nSequence(s) were found to only contain gaps. They are in $errfile\n";
}

# Sanity check
print "Job complete. Summary:\n";
print "\tInput:   $checkIn sequence/alleles\n";
print "\tOutput:  $checkOut sequence/alleles\n";
print "\tRemoved: $checkErr sequence/alleles\n";
print "\nOutfile is $outfile\n\n";

############################
sub print_usage {
    my ($error) = @_;
    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 [fasta file]\n";
    print "\tInput: fasta-formatted >seqid\\nseq\\n file\n";
    print "\tOutput: same format, just without gaps\n";
    print "\tAt this point in time, code specifically removes - and !\n";
    print "\tIt will remove all sequences that only contain gaps\n";
    print "\tCheers!\n\n";
    exit 1;
}
