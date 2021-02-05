#!/usr/bin/perl
use strict;
use warnings;

#declare variables
my $ref;
my $query;
my $start = 0;
my $stop; 
my $buffer;
my @queryfile = ();
my %uniqseqs = ();
my $row;
my $key;
my $seqID;
my $ref_ID;
my $reference_seq;
my $freq;
my $mutnumber;
my $mutation;
my $query_nt;
my $ref_nt;
my $mutstart = "NA";
my $mutstop = "NA";
my $length;
my $mutlength = "NA";
my %mutations = ();
my $index;
my $tmut ;
my $tSNP;
my $tgap;
my $seq_length;
my $c = 0;
my $k;
my $seq;
my @query = ();
my @ref = ();
my %output = ();
my $pos;
my $i;
my $allele_ID = 1;
my $summary_1;
my $summary_2;
my $summary_3;
my $outfile = "uniqseqs.fasta";
my $b_start = 0;
my $b_stop;
my $greedy = "FALSE";

# Get options
$query = $ARGV[0];

                    ##################################
                    # Sequence Input and Preparation #
                    ##################################

# This will save all unique sequences for the specified coordinates as 
# hash keys, and then append each of the sequence identifiers to the end
# of the value associated with the key. The reference will be saved as a 
# single string that is then split into an array. Assumes both inputs contain
# the sequence ID on one line and the sequence itself on the immediate next line,
# only spanning one line in its entirety. Will need to build in a safeguard. 

# Open input file
if (-f $query) {
    open F, '<', $query; 
        while ($row = <F>) {
            chomp $row; #removes newline (unseen) at the end of a line
            if ($row =~ /^#/) {
                next;
            }   
            if ($row =~ s/^>//) {
                $seqID = $row;
                $seqID =~ s/^\>//d;
            }   
            else {
                $key = $row;
                $uniqseqs{$key} .= "$seqID ";
            }
        }
    close F;
} else {print"\n error in read file\n";}

                          ###########################
                          # Generating Output Table #
                          ###########################

open F, '>', $outfile or die "problem saving output to file\n";

foreach $key (keys %uniqseqs) {
    print F ">$uniqseqs{$key}\n";
    print F "$key\n";
}
close F;

##########################################
sub print_usage {
    my ($error) = @_; 

    if (defined $error) {
        print STDERR $error, "\n";
    }   

    print "\nUsage:  uniqallele.pl -query [ALIGN.FASTA] -ref [REF.FASTA] -start [START] -stop [STOP] -buffer [BUFFER] -greedy [TRUE]\n";
    print "\tquery:  aligned fasta files for analysis\n";
    print "\tref:    fasta file for comparison against input\n";
    print "\tstart:  the beginning location (nt) for gene of interest (index 1)\n";
    print "\tstop:   the end location (nt) for gene of interest (index 1)\n";
    print "\tbuffer: distance around gene of interest (nt) for inclusive analysis\n";
    print "\tgreedy: will take entire sequence length, no need for start/stop to be specified\n";
    print "\tnote that all coordinates are indexed by 0 in output.\n";
    print "\nCheers!\n\n";
    exit; 
}

