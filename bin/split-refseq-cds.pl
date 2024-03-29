#!/usr/bin/perl
use strict;
use warnings;

# Assumes onelineseq.fmt
# Edited to allow for duplicate gene names

# Declare vars
my ($geneID, $seq, $outfile);
my ($i, $k);

my $counter = 0;
my $input = $ARGV[0];
my $gene_pattern = "gene";
my $locus_pattern = "locus_tag";

my @sequence_data = ();
my @temp = ();
my %registered_genes = ();

# Verify input
if (!defined $ARGV[0] || !defined $input) {
    &print_usage("\nPlease specify a file of coding sequences.");
}

# Eat input
open F, '<', $input;
    @sequence_data = <F>;
close F;

# Iterate through input; find seqID then get sequence
for ($i = 0; $i < scalar(@sequence_data); $i++) {
    if ($sequence_data[$i] =~ m/[0-9]|[a-z]/) {
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
    }
    else { 
        $seq = $sequence_data[$i];
        $seq =~ s/\R//g;
        $outfile = $geneID . ".nt";
        open F, '>', $outfile or die "problem saving to file $outfile\n";
            print F ">$geneID\n";
            print F "$seq";
        close F;    
        $counter++;
        $geneID = ();
    }
}

# Final report
print "\nFiles with duplicate gene names are appended with (n).\n";
print "\tchecksum: file had ". scalar(@sequence_data) ." lines, or ". scalar(@sequence_data)/2 ." sequences\n";
print "\t$counter files were printed\n";
if ($counter == scalar(@sequence_data)/2) {
    print "\tThings should be good to go\n\n";
}
else {
    print "\t!!!!Something isn't quite right!!!!\n\n";  
    exit 1;
}

############################
sub print_usage {
    my ($error) = @_;
    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 [2lnseq nt features refseq file]\n";
    print "Will split features into individual files with .nt ext\n";
    print "assumes that tags are on one line, seqs on another\n";
    print "If files are not formatted, please use one-line-seq.pl\n";
    print "Cheers!\n\n";
    exit;
}
