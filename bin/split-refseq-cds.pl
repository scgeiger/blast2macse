#!/usr/bin/perl
use strict;
use warnings;

#Assumes onelineseq.fmt

my $input = $ARGV[0];
my $geneID;
my $seq;
my $outfile;
my $gene_pattern = "gene";
my $locus_pattern = "locus_tag";
my $i;
my $k;
my @sequence_data = ();
my @temp = ();

# Verify input
if (!defined $ARGV[0] || !defined $input) {
    &print_usage("\nPlease specify a file of coding sequences.");
}

open F, '<', $input;
    @sequence_data = <F>;
close F;

for ($i = 0; $i < scalar(@sequence_data); $i++) {
    if ($sequence_data[$i] =~ m/[0-9]|[a-z]/) {
        @temp = split / /, $sequence_data[$i];
        for $k (@temp) {
            next unless ($k =~ m/$gene_pattern/ || $k =~ m/$locus_pattern/);
            $geneID = $k;
            $geneID =~ s/^[^=]*\=//;
            $geneID =~ s /\///;
            $geneID =~ s/\]//;
            print "$geneID\n";
            last;
        }
    }
    else { 
        $seq = $sequence_data[$i];
        $seq =~ s/\R//g;
        $outfile = $geneID . ".nt";
        open F, '>', $outfile or die "problem saving to file $outfile\n";
            print F ">$geneID\n";
            print F "$seq";
        close F;
        $geneID = ();
    }
}

############################
sub print_usage {
    my ($error) = @_;
    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 [2lnseq nt features refseq file]\n";
    print "Will split features into individual files with .nt ext\n";
    print "Cheers!\n\n";
    exit;
}





