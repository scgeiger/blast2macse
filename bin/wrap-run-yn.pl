#!/bin/perl
use strict;
use warnings;
#use File::Temp;

# This script is for running swaine's run-yn.pl script on one file
# In the future, build in shortcut so it can be run on uniqueseqs_NT file
#   all seqIDs are in header for each uniq seq; can run calc once and copy over
# Sequences will be de-aligned before running
# Assumes 2ln format
# Outputs will (in final implementation) include: 
#   matrix of all dN dS calculations
#   Rplot showing distribution of pairwise scores
#   Average dN/dS for the whole dataset

# Declare vars
my $inFile = $ARGV[0];
my $i;
my $n;
my $j;
my $nSeq; 
my @inFile = ();
my $filename;
my $row;
my $iSeq;
my $jSeq;
my $iSeqID;
my $jSeqID;
my @tempout = ();
my $tempout;
my $dNdS;
my @matrix;
my $header;
my @out = ();
my $outfile;
my @alldnds = ();
my $meandnds;
my $total;
my $summaryfile;
my $dn_and_ds_are_0; # is marked by -2
my $only_ds_is_0;    # is marked by -1
my $ERRs;            # likely a premature stop codon
my $pairwise_count;  # number of total comparisons

($filename) = ($inFile =~ /([^.]+)/);
$outfile = $filename . "-dnds-matrix.tsv";
$summaryfile = $filename . "-dnds-summary.tsv";

# Take in input file
if (!defined $ARGV[0] || !defined $inFile) {
    &print_usage("\nPlease specify an expanded uniqseqs file (2ln format).");
}

if ( -f $inFile ) {
    open F, '<', $inFile;
        while ($row = <F>) {
            chomp $row;
            $row =~ s/^>//;
            $row =~ s/!/-/g;  # just in case ! is from MACSE
            push @inFile, $row;
        }
    close F;
}

# Initialise important variables
$nSeq = @inFile / 2;
print "nseq is $nSeq\n"; 
$n = 0; # will track progress in row index 
$dn_and_ds_are_0 = 0; # is marked by -2, both dn/ds = 0
$only_ds_is_0 = 0;    # is marked by -1, only ds = 0
$ERRs = 0;            # likely a premature stop codon
$pairwise_count = 0;  # number of total comparisons

open O, '>', $outfile;
    for ($i = 0; $i < $#inFile - 2; $i += 2) {
        print "looking at $n of $nSeq\n";
        @out = ();
        $iSeqID = $inFile[$i];
        $iSeq   = $inFile[$i + 1];
        $header = "";
        $header = "\t$iSeqID";
        $iSeqID =~ s/^\s+//;
        push @out, $iSeqID;
        push @out, "-" for 0..$n;
        $n++;

        for ($j = $i + 2; $j < $#inFile; $j += 2) {
            $jSeqID = $inFile[$j];
            $jSeq   = $inFile[$j + 1];
            $header .= "\t$jSeqID"; #only thru each iteration
    
            open T, '>', "tempseq.fasta";
                print T ">$iSeqID\n$iSeq\n>$jSeqID\n$jSeq\n";
            close T;

            @tempout = `perl /usr/local/bin/run-yn.pl -noalign tempseq.fasta`;
            $pairwise_count++;
            $tempout = $tempout[1];
            @tempout = split("\t", $tempout);
            $dNdS = $tempout[1];
            if ($dNdS eq "-") {
                $dNdS = "ERR";
            }
            push @out, $dNdS; 
            
            if ($dNdS eq "ERR") {
                $ERRs++;
            } 
            elsif ($dNdS == -1) {
                $only_ds_is_0++;
            }
            elsif ($dNdS == -2) {
                $dn_and_ds_are_0++;
            } else { 
                push @alldnds, $dNdS;
            }
        }
        if ( -z $outfile ) {
            print O "$header\n";
        }
        print O join("\t", @out);
        print O "\n";
    }
close O; 

foreach (@alldnds) {
  $total += $_;
}

$meandnds = $total / @alldnds;
my $dndscalcs = @alldnds;

# print summary
open S, '>', $summaryfile;
    print S "nseqs\tpairwisecomparisons\tERRs\tonlyds=0\tdndsboth=0\tmeandnds\tdndscalcs\n";
    print S "$nSeq\t$pairwise_count\t$ERRs\t$only_ds_is_0\t$dn_and_ds_are_0\t$meandnds\t$dndscalcs";
close S;
    
print "Final report: $nSeq Sequences\n";
print "\t$pairwise_count comparisons were made in total\n";
print "\t$ERRs of these were errors (check for premature stop codons)\n";
print "\t$only_ds_is_0 of these had dS = 0, but dN had a value\n";
print "\t$dn_and_ds_are_0 of these had dS = 0 and dN = 0\n";
print "\tThe mean dnds is $meandnds; total is $total, with $dndscalcs values\n";

system("rm -f tempseq.fasta");

# plot distribution of dnds scores using Rscript 

# Subroutines
sub print_usage {
    my ($error) = @_;
    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 [am-*uniqseqs.aln]\n";
    print "\tInput: path to *uniqseqs_NT.aln files)\n";
    print "\nCheers!\n\n";
    exit 1;
}
