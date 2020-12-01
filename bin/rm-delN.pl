#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# The purpose of this script is to remove all sequences that contain 
# insertions or deletions that could interfere with PopGenome calculations.
# Relative to the consensus, these mutations were called in codon-dist. If a sequence contains deletions or Ns, remove those sequences from consideration. 

#input: macse alignment
#input: codon-dist.tsv

#output: seq-rem-PG
my $row;
my $dist_file;
my $aln_file;
my %remove = ();
my $outfile_aln;
my $outfile_rem;
my @row = ();
my $nRem;
my $header;
my $seqID;
my $prefix;

GetOptions (
    'aln=s'   => \$aln_file,
    'dist=s'  => \$dist_file,
);

# make sure inputs are correct
if (!defined $aln_file || !-f $aln_file) {
    &print_usage("\nPlease specify an aligned nucleic acid file for filtering.");
}

if (!defined $dist_file || !-f $dist_file) {
    &print_usage("\nPlease specify a macse codon dist file.");
}
$prefix = $dist_file;
$prefix =~ s/\..*//;

$outfile_aln = $prefix . "-removed-macse.aln";
$outfile_rem = $prefix . "-removed-macse-ids.tsv";

# Open codon dist
if (-f $dist_file) {
    open F, '<', $dist_file;
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            @row = split ( /\t/, $row);
            if ($row[3] eq "Ambiguous" ||
                $row[3] eq "Deletion") {
                $remove{$row[12]}++;
            }
        }
    close F;
}

$nRem = scalar(keys(%remove));
$header = "#aln file: $aln_file; dist file: $dist_file; nRemoved $nRem\n";
# open macse alignment
if (-f $aln_file) {
    open F, '<', $aln_file;
    open ALN, '>', $outfile_aln;
    open SEQ, '>', $outfile_rem;
        print ALN $header;
        print SEQ $header;
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            if ($row =~ s/^>//) {
                $seqID = $row;
                $seqID =~ s/^\>//d;
                if (exists($remove{$seqID})) {
                    print "removing $seqID: has $remove{$seqID} abberations.\n";
                    print SEQ "$seqID\t$remove{$seqID}\n";
                    $seqID = ();
                    next;
                }
            }
            elsif ($seqID) {
                 print ALN ">$seqID\n$row\n"; 
            }
        }
    close SEQ;
    close ALN;
    close F;
}
            
# Save seqid and seq to new file ONLY if seqid is not in hash
# Close

##### usage subroutine
sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 -aln [ALIGN.MACSE] -dist [MACSE CODON DIST]\n";
    print "\tThe script is being tailored to accomodate MACSE nt output\n";
    print "\taln:  macse aln file\n";
    print "\tdist: macse codon dist file\n";
    print "\tThe output can be used to calculate PopGenome stats more reliably\n";
    print "\tMake sure you're using the expanded alignment file.\n";
    print "\nCheers!\n\n";
    exit;
}
