#!/usr/bin/perl
use strict;
use warnings;

# Assumes onelineseq.fmt
# Edited to allow for duplicate gene names
# Same core script as split-refseq-cds.pl except only saves pseudo genes
# Changed it so print happened simultaneous to analysis

# Declare vars
my ($geneID, $seq, $outfile, $inputID);
my ($i, $k);

my $counter = 0;
my $input = $ARGV[0];
my $gene_pattern = "gene";
my $locus_pattern = "locus_tag";
my $protein = "protein";

my @sequence_data = ();
my @temp = ();
my %registered_genes = ();
my @pseudo;

# Verify input
if (!defined $ARGV[0] || !defined $input) {
    &print_usage("\nPlease specify a file of coding sequences.");
}

# Eat input
open F, '<', $input;
    @sequence_data = <F>;
close F;

$inputID = $input;
$inputID =~ s/\..*//;
$outfile = $inputID . "-pseudolist.txt";

# Iterate through input; find seqID then get sequence
open F, '>', $outfile or die "problem saving to file $outfile\n";
    for ($i = 0; $i < scalar(@sequence_data); $i++) {
        if ($sequence_data[$i] =~ m/[0-9]|[a-z]/ && $sequence_data[$i] =~ m/pseudo=/) {
            @temp = split / /, $sequence_data[$i];
            for $k (@temp) {
                next unless ($k =~ m/$gene_pattern/ || $k =~ m/$locus_pattern/);
                $geneID = $k;
                $geneID =~ s/^[^=]*\=//;
                $geneID =~ s /\///;
                $geneID =~ s/\]//;
                if (exists $registered_genes{$geneID}) {
                    $registered_genes{$geneID}++;
                    $geneID = "$geneID" . "($registered_genes{$geneID})";
                }
                else {
                    $registered_genes{$geneID} = 0;
                }
                last;
            }
        #push(@pseudo, $geneID);
        print F "$geneID\n";
        $counter++;
        $geneID = ();
        }
    }
close F;

# Final report
print "\nFiles with duplicate gene names are appended with (n).\n";
print "\t$counter pseudogenes were found.\n\n";

############################
sub print_usage {
    my ($error) = @_;
    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 [2lnseq nt features refseq file]\n";
    print "Will collect sequence IDs for pseudogenes.\n";
    print "Assumes that tags are on one line, seqs on another\n";
    print "Cheers!\n\n";
    exit;
}
