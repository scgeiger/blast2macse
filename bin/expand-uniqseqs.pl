#!/usr/bin/perl
use warnings; 
use strict;

# Last edited 210827. Cleaned and added exit 1
# use this to take a file of extracted uniq alleles and expand so one seq, one line 

# Declare variables
my ($in_file, $seq_file, $outfile, $id, $seqid, $row, $key);
my @key_array = ();
my %in_hash = ();
my %seqs = ();

# Verify input
if (!defined $ARGV[0]) {
    &print_usage("\nPlease specify uniqseqs fasta file.");
} 

$in_file = $ARGV[0];
$outfile = "expanded-". $in_file;

# Get file contents
if (-f $in_file) {
    open I, '<', $in_file;
        while ($row = <I>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            if ($row =~ s/^>//) {
                $seqid = $row;
                $seqid =~ s/^\>//d;
            }
            else {
                $in_hash{$seqid} = $row;
            }
        }
    close I;
} else {&print_usage("Did you specify a uniqseqs output file?")}

# Print Output
open O, '>', $outfile;
    foreach $key (keys %in_hash) {
        @key_array = split (/ /, $key);
        foreach $id (@key_array) {
            print O ">$id\n";
            print O "$in_hash{$key}\n";
        }
     }
close O;

print "\nJob complete. Outfile is $outfile\n\n";

############################
sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 [uniqseqs file]\n";
    print "\tInput: Uniqseqs file with seqIDs in the header\n";
    print "\tOutput: Fasta file with each sequence separated\n";
    print "\nCheers!\n\n";
    exit 1;
}
