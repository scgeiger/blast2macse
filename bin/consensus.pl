#!/usr/bin/perl
use strict;
use warnings;

# Last edited 210827. Cleaned and added exit 1
# Combined loops cleaning data + reading it from file
# Assumes sequences are aligned. 

# Declare Variables
my ($row, $pos, $most_freq, $output, $key, $largest, $length, $yardstick);
my $input_file = $ARGV[0];
my $outfile = "consensus.fasta";

my @input = ();
my @consensus;
my @seq = ();
my @seq_data = ();

# Verify input
if (!defined $ARGV[0]) {
    &print_usage("\nPlease specify an aligned fasta file.");
}

# Get file contents
if (-f $input_file) {
    open F, '<', $input_file;
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            elsif ($row =~ /^>/) {
                next;
            }
            else { 
                $length = length($row);
                if (!defined $yardstick) {
                    $yardstick = length($row);
                }
                if ($length == $yardstick) {
                    @seq = ();
                    @seq = split(//, $row);
                    for ($pos = 0; $pos <= $#seq; $pos++) {
                        $consensus[$pos]->{$seq[$pos]} = 0 unless defined $consensus[$pos]->{$seq[$pos]};
                        $consensus[$pos]->{$seq[$pos]}++;
                    } 
                }
                else {
                    &print_usage("\nLengths of sequences are not equal. Are seqs aligned?\n"); 
                }
            }
        }
    close F;
} else {&print_usage("Please specify a fasta file with aligned sequences")}

# Figure out consensus sequence
foreach $pos (0..$#consensus) {
    $largest = 0;
    $most_freq = ();
    foreach $key (keys %{$consensus[$pos]}) {
        if ($consensus[$pos]->{$key} > $largest) {
            $largest = $consensus[$pos]->{$key};
            $most_freq = $key;
        }
    }
    $output .= $most_freq;
}

# Print Outfile 
open F, '>', $outfile or die "problem saving to file\n";
    print F ">$input_file" ."_consensus\n";
    print F $output;
close F;

print "\nConsensus sequence has been generated, saved in $outfile\n\n";

###########################
sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "Usage: consensus.pl [ALIGN.FASTA]\n";
    print "\tInput: aligned fasta file\n";
    print "\tOutput: consensus sequence for the fasta file\n";
    print "\nThis code can handle DNA or AA sequences, so long as they are aligned\n";
    print "\nCheers!\n\n";
    exit 1;
}
