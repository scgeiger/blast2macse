#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# This code runs a simple analysis, looking at mutations in an indexed manner. Runs rather quickly, and is 
# supposed to make plotting much more rapid. It's been checked against my short.fasta dataset for 
# functionality and accuracy.
# 200617 Does not handle location calling (start stop) properly. 
# 200906 Fixed start/stop issue. Integrated ambiguous calls and aa-typing into this script.

#declare variables
my ($ref, $query, $start, $stop, $greedy, $buffer);  # Input variables
my ($r_nt, $r_aa, $q_nt, $q_aa); # Ref and query scalars
my (%TScoords, %TVcoords, %AMcoords, %INcoords, %DELcoords, %FScoords) = (); #NT mut coords
my (%TScontig, %TVcontig, %AMcontig, %INcontig, %DELcontig, %FScontig) = (); #contiguous muts
my ($key, $row, $pos, $seqID); # Scalar indexes
my (%nt_output, %aa_output) = ();
my (%nt_uniqseqs, %aa_uniqseqs) = ();
my ($nt_reference_seq, $aa_reference_seq);
my ($tmut, $tTS, $tTV, $tAM, $tIN, $tDEL, $tFS, $tMS, $tNS, $tSYN);
my (%MScontig, %NScontig, %SYNcontig);
my (%MScoords, %NScoords, %SYNcoords);
my (@q_codon, @r_codon, @all_q_codons, @all_r_codons) = ();
my ($codon_start, $codon_stop, $codon_b_start, $codon_b_stop);

# Initializations
$start = 0;
my ($mutstart, $mutstop, $mutlength) = "NA";
my @queryfile = ();
my %uniqseqs = ();
my $refID; # not used
my $freq;
my $mutnumber; # checked
my $mutation;
my $length;
my %mutations = ();
my $index;
my $seq_length;
my $k;
my $seq;
my @query = ();
my @ref = ();
my %output = ();
my $i;
my $allele_ID = 1;
my ($summary_1, $summary_2, $summary_3);
my $outfile = "nt-uniq-macse-output.tsv";
my $aa_outfile = "aa-uniq-macse-output.tsv";
my $b_start = 0;
my $b_stop;
$greedy = "FALSE";
my @uniqcoords = ();
my $mut_start;
my $mut_stop;
my @mutant_start = ();
my @sorted_starts = ();

# Get options
GetOptions (
    'ref=s'   => \$ref,     # Reference fasta file
    'query=s' => \$query,   # Query fasta file, aligned
    'start=i' => \$start,   # Start location for gene of interest (aa) (index 1)
    'stop=i'  => \$stop,    # Stop location for gene of interest (aa) (index 1)
    'buffer=i'=> \$buffer,  # Distance around gene for analysis (aa)
    'greedy=s'=> \$greedy,  # Take whole alignment or not
);

# Verify input
if (!defined $ARGV[0] && (!defined $query || !-f $query)) {
    &print_usage("\nPlease specify an aligned macse nt file for analysis.");
}

if (!defined $ARGV[1] && (!defined $ref || !-f $ref)) {
    &print_usage("\nPlease specify a reference file.");
}

if ((!defined $start || !defined $stop) && 
    (!defined $greedy || (defined $greedy && $greedy ne "TRUE"))) {
    print "No start or stop was specified. Please specify start and stop, or greedy as TRUE.\n";
    exit;
}

if (!defined $buffer) {
    $buffer = 0;
}

if (defined $start && defined $stop) {
    $b_start = $start - $buffer - 1;# adjust to buffer, as well as index 0
    $b_stop = $stop + $buffer - 1;  # adjust to buffer, as well as index 0

    if ($b_start < 0) { #to keep the var from dipping into the negatives
        print "Buffered window is out of range. Window minimum is adjusted to 0.\n";
        $b_start = 0;
    }
    if ($b_stop < 0) {
        die "Error: final length for analysis is negative. Something is very wrong.\n";
    }

    $seq_length = $b_stop - $b_start;
}


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
                if ($greedy eq "TRUE") {
                    $b_start = 0;
                    $start = 1; ###SCG this was start = 1, but having indexing issues
                    $b_stop = length($row) - 1; #length takes index 1, b_stop is index 0
                    $stop = length($row);       #length takes index 1 
                    $seq_length = $stop - $start;
                } #greedy clause is correct up to here. The whole seq is digested.

                if (($b_start + $seq_length) > (length($row) - 1)) {
                    print "You've selected a buffer that extends past the sequence end.\n";
                    print "Game over. Please try again.\n";
                    print "Remind me, and I might build in an adjustment.\n";
                    print "If you're sure things are okay, make sure ref has entire sequence on one line\n";
                    die;
                } 
                #$row = uc($row);
                $key = substr $row, $b_start, $seq_length + 1; # +1 keeps last pos
                $nt_uniqseqs{$key} .= "$seqID ";
            }
        }
    close F;
} else {print"\n error in read file\n";}

# Open reference; assumes entire reference sequence is one line
if (-f $ref) {
    open F, '<', $ref;
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }   
            if ($row =~ s/^>//) {
                $refID = $row;
                $refID =~ s/^\>//d;
            }   
            else {
                if ($greedy eq "TRUE") {
                    $b_start = 0; 
                    $start = 1;
                    $b_stop = length($row) - 1; 
                    $stop = length($row); 
                    $seq_length = $b_stop - $b_start;
                }
                if (($b_start + $seq_length) > (length($row) - 1)) {
                    print "You've selected a buffer that extends past the sequence end.\n";
                    print "Game over. Please try again.\n";
                    print "Remind me, and I might build in an adjustment.\n";
                    die;
                } 
#                $row = uc($row);
                $nt_reference_seq = substr $row, $b_start, $seq_length + 1;
            } #tested with greedy. Reference is consumed properly.
        }   
    close F;
}
@ref = split //, $nt_reference_seq;

print "\nsequences have all been collected successfully\n";


                          ##############################
                          # Complete Sequence Analysis #
                          ##############################

# For each unique sequence present in the coordinates, assign values to
# the eventual row in the output table. $index refers to the allele ID 
# number, the start and stop refer to the gene start and stop loci, freq
# takes the appended sequence IDs (assumes no spaces in filenames) and
# enumerates the number of sequences. $mutnumber is initialized to 1, and
# the query is split into an array. If everything matches the reference
# explicitly, it is initialized differently. Should not break down if nothing
# matches the reference. $mutnumber tracks the number of individual mutations
# in a single allele.


#### Nucleotide Sequence Analysis
for $seq (keys %nt_uniqseqs) {
    $index = 0;    # Indexes mutations within an allele
    $nt_output{$seq}->{"Start"}  = $start;    # Index 1
    $nt_output{$seq}->{"Stop"}   = $stop;     # Index 1
    $nt_output{$seq}->{"Bstart"} = $b_start;  # Index 0
    $nt_output{$seq}->{"Bstop"}  = $b_stop;   # Index 0
    $nt_output{$seq}->{"Freq"}   = $nt_uniqseqs{$seq} =~ tr{ }{ };
    $nt_output{$seq}->{"SeqIDs"} = $nt_uniqseqs{$seq};  
    $mutnumber = 1; 
    @query = split //, $seq;
    %TScoords  = ();  # transition
    %TVcoords  = ();  # transversion
    %INcoords  = ();  # insertion 
    %DELcoords = ();  # in-frame deletions
    %FScoords  = ();  # frameshift
    %AMcoords  = ();  # ambiguous
    @uniqcoords   = ();
    @mutant_start = ();
    
    if ($seq eq $nt_reference_seq) { # if the query is same as reference
        $mutations{$seq}->[$index]->{"Mutstart"} = "NA";
        $mutations{$seq}->[$index]->{"Mutstop"}  = "NA";
        $mutations{$seq}->[$index]->{"Mutlength"}= "NA";
        $mutations{$seq}->[$index]->{"Muttype"}  = "NA";
        $mutations{$seq}->[$index]->{"Mutnumber"}= "NA";
        $nt_output{$seq}->{"tAM"}  = 0;
        $nt_output{$seq}->{"tTV"}  = 0;
        $nt_output{$seq}->{"tTS"}  = 0;
        $nt_output{$seq}->{"tIN"}  = 0;
        $nt_output{$seq}->{"tDEL"} = 0;
        $nt_output{$seq}->{"tFS"}  = 0;
        $nt_output{$seq}->{"tmut"} = 0;
        next;
    } 
    # for each position in the unique sequence, check to see if it matches
    # the reference. If it does, then the presense of mutation is F, and 
    # we move on to the next. If there is a difference and query has a - at
    # that position, then this is recorded as a gap.  
    else {
        foreach $pos (0..$#query) { 
        $q_nt = $query[$pos];
        $r_nt = $ref[$pos];
            if ($q_nt eq $r_nt) {
                next;
            }
            if ($q_nt eq "N" || $r_nt eq "N") {
                $AMcoords{$pos} = 1; #Ambiguous
            }
            elsif ($q_nt ne "-" && $q_nt ne "!" && 
                  ($r_nt eq "-" || $r_nt eq "!")) {
                $INcoords{$pos} = 1; #insertion
            } 
            elsif ($q_nt eq "-" && $r_nt ne "-") {
                $DELcoords{$pos} = 1; #deletion
            }
            elsif ($q_nt eq "!") {
                $FScoords{$pos} = 1;
            } 
            elsif (($q_nt eq "A" && $r_nt eq "G") ||   #At this point, should be SNP of some kind
                   ($q_nt eq "G" && $r_nt eq "A") ||
                   ($q_nt eq "C" && $r_nt eq "T") ||
                   ($q_nt eq "T" && $r_nt eq "C")) {
                   $TScoords{$pos} = 1; # Is a transition
            }
            elsif (($q_nt eq "A" && $r_nt eq "C") ||
                   ($q_nt eq "C" && $r_nt eq "A") ||
                   ($q_nt eq "A" && $r_nt eq "T") ||
                   ($q_nt eq "T" && $r_nt eq "A") ||
                   ($q_nt eq "G" && $r_nt eq "T") ||
                   ($q_nt eq "T" && $r_nt eq "G") ||
                   ($q_nt eq "G" && $r_nt eq "C") ||
                   ($q_nt eq "C" && $r_nt eq "G")) {
                   $TVcoords{$pos} = 1; #is transversion
            }
            else { print "Unknown nt case at $pos: ref is $r_nt and query is $q_nt\n";
                   print "seqids: $nt_uniqseqs{$seq}\n";
                   print "this case will be counted as ambiguous\n";
                   $AMcoords{$pos} = 1; #Ambiguous
             };
        }
    }

    %TScontig = ();
    if (%TScoords) {
        @uniqcoords = sort {$a <=> $b } keys %TScoords;
        $mut_start = ();
        foreach $i (0..$#uniqcoords) {
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $start;
            }
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $start;
                    $TScontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $start;
                $TScontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            } 
        }
    }

    %TVcontig = ();
    if (%TVcoords) {
        @uniqcoords = sort {$a <=> $b } keys %TVcoords;
        $mut_start = ();
        foreach $i (0..$#uniqcoords) {
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $start;
            }
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $start;
                    $TVcontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $start;
                $TVcontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }

    %AMcontig = ();
    if (%AMcoords) {
        @uniqcoords = sort {$a <=> $b } keys %AMcoords;
        $mut_start = ();
        foreach $i (0..$#uniqcoords) {
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $start;
            }
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $start;
                    $AMcontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $start;
                $AMcontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }

    %INcontig = ();
    if (%INcoords) {
        @uniqcoords = sort {$a <=> $b } keys %INcoords;
        foreach $i (0..$#uniqcoords) {
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $start;
            }
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $start;
                    $INcontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $start;
                $INcontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }
    
    %DELcontig = ();
    if (%DELcoords) {
        @uniqcoords = sort {$a <=> $b } keys %DELcoords;
        foreach $i (0..$#uniqcoords) {
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $start;
            }
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $start;
                    $DELcontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {   #other code just has an if here...
                $mut_stop = $uniqcoords[$i] + $start;
                $DELcontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }

    %FScontig = ();
    if (%FScoords) {
        @uniqcoords = sort {$a <=> $b } keys %FScoords;
        foreach $i (0..$#uniqcoords) {
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $start;
            }
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $start;
                    $FScontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i];
                $FScontig{$mut_start} = $mut_stop + $start;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }

    @mutant_start = (); 
    push (@mutant_start, keys(%TScontig));
    push (@mutant_start, keys(%TVcontig));
    push (@mutant_start, keys(%AMcontig));
    push (@mutant_start, keys(%INcontig));
    push (@mutant_start, keys(%DELcontig));
    push (@mutant_start, keys(%FScontig));
    @sorted_starts = sort {$a <=> $b} (@mutant_start);

    foreach $i (@sorted_starts) {
        if (exists $TScontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $TScontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($TScontig{$i} - $i) + 1; #counting # nt changed; single poly counted 
            $mutations{$seq}->[$index]->{"Muttype"}  = "TS";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber; 
            $mutnumber++;
            $index++;
        }
        elsif (exists $TVcontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $TVcontig{$i}; 
            $mutations{$seq}->[$index]->{"Mutlength"}= ($TVcontig{$i} - $i) + 1; #counting # nt changed; single poly counted
            $mutations{$seq}->[$index]->{"Muttype"}  = "TV";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        elsif (exists $AMcontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $AMcontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($AMcontig{$i} - $i) + 1;
            $mutations{$seq}->[$index]->{"Muttype"}  = "AM";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        elsif (exists $INcontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $INcontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($INcontig{$i} - $i) + 1;
            $mutations{$seq}->[$index]->{"Muttype"}  = "IN";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        elsif (exists $DELcontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i; 
            $mutations{$seq}->[$index]->{"Mutstop"}  = $DELcontig{$i}; 
            $mutations{$seq}->[$index]->{"Mutlength"}= ($DELcontig{$i} - $i) + 1;
            $mutations{$seq}->[$index]->{"Muttype"}  = "DEL";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        elsif (exists $FScontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $mutstart = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $FScontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($FScontig{$i} - $i) + 1;
            $mutations{$seq}->[$index]->{"Muttype"}  = "FS";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        else { 
            print "couldn't re-trace mutstart $i\n"; 
        }
    }
     
    $_ = 0 for ($tTS, $tTV, $tAM, $tIN, $tDEL, $tFS);
    
    foreach $i (0..$#{$mutations{$seq}}) {
        if ($mutations{$seq}->[$i]->{"Muttype"} eq "TS") {
            $tTS += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "TV") {
            $tTV += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "AM") {
            $tAM += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "IN") {
            $tIN += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "DEL") {
            $tDEL += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "FS") {
            $tFS += $mutations{$seq}->[$i]->{"Mutlength"};
        }
    }
    
    $nt_output{$seq}->{"tTS"} = $tTS;
    $nt_output{$seq}->{"tTV"} = $tTV;
    $nt_output{$seq}->{"tAM"} = $tAM;
    $nt_output{$seq}->{"tIN"} = $tIN; 
    $nt_output{$seq}->{"tDEL"}= $tDEL;
    $nt_output{$seq}->{"tFS"} = $tFS;
    $nt_output{$seq}->{"tmut"}= $tTS+ $tTV+ $tAM + $tIN + $tDEL + $tFS;

} 

                          ###########################
                          # Generating Output Table #
                          ###########################

print "\nfinal nt output generated\n";

open F, '>', $outfile or die "problem saving output to file\n";
    print F "allele\tfreq\ttmut\ttTS\ttTV\ttAM\ttIN\ttDEL\ttFS\tmutnum\ttype\tmstart\tmstop\tmlength\tstart\tbstart\tstop\tbstop\tseqs\n";

foreach $key (sort {$nt_output{$a}->{"Freq"} <=> $nt_output{$b}->{"Freq"} ||
                    $nt_output{$a}->{"tmut"} <=> $nt_output{$b}->{"tmut"}} keys %nt_output)  {
    $summary_1 =$allele_ID . "\t" .
                $nt_output{$key}->{"Freq"} . "\t" .
                $nt_output{$key}->{"tmut"} . "\t" .
                $nt_output{$key}->{"tTS"} . "\t" .
                $nt_output{$key}->{"tTV"} . "\t" .
                $nt_output{$key}->{"tAM"} . "\t" .
                $nt_output{$key}->{"tIN"} . "\t" .
                $nt_output{$key}->{"tDEL"} . "\t" .
                $nt_output{$key}->{"tFS"} . "\t"; 

    $summary_3 =$nt_output{$key}->{"Start"}  . "\t" . 
                $nt_output{$key}->{"Bstart"} . "\t" .
                $nt_output{$key}->{"Stop"}   . "\t" .  
                $nt_output{$key}->{"Bstop"}  . "\t" .
                $nt_output{$key}->{"SeqIDs"};

    for ($i = 0; $i <= $#{$mutations{$key}}; $i++) {
        $summary_2 = $mutations{$key}->[$i]->{"Mutnumber"} . "\t" .
                     $mutations{$key}->[$i]->{"Muttype"}   . "\t" .
                     $mutations{$key}->[$i]->{"Mutstart"}  . "\t" .
                     $mutations{$key}->[$i]->{"Mutstop"}   . "\t" .
                     $mutations{$key}->[$i]->{"Mutlength"} . "\t";
        print F $summary_1 . $summary_2 . $summary_3 . "\n";            
    }
    $allele_ID++;   
}
close F;


                          ##############################
                          # Complete Sequence Analysis #
                          ##############################
print "looking at amino acid mutations now\n";
%mutations = ();
@all_r_codons = ( $nt_reference_seq =~ m/...?/g );
$allele_ID = 1;

for $seq (keys %nt_uniqseqs) {
    @all_q_codons = ( $seq =~ m/...?/g );
    $index = 0;    # Indexes mutations within an allele
    $codon_start = (($start - 1) / 3) + 1;
    $codon_stop = ($stop / 3); 
    $codon_b_start = ($b_start / 3);
    $codon_b_stop = (($b_stop + 1) / 3) - 1;

    $aa_output{$seq}->{"Start"}  = $codon_start;    # Index 1
    $aa_output{$seq}->{"Stop"}   = $codon_stop;     # Index 1
    $aa_output{$seq}->{"Bstart"} = $codon_b_start;  # Index 0
    $aa_output{$seq}->{"Bstop"}  = $codon_b_stop;   # Index 0
    $aa_output{$seq}->{"Freq"}   = $nt_uniqseqs{$seq} =~ tr{ }{ };
    $aa_output{$seq}->{"SeqIDs"} = $nt_uniqseqs{$seq};
    $mutnumber = 1;

    %MScoords  = ();  # missense
    %NScoords  = ();  # nonsense
    %SYNcoords = ();  # transversion
    %AMcoords  = ();  # ambiguous
    %INcoords  = ();  # insertion
    %DELcoords = ();  # in-frame deletions
    %FScoords  = ();  # frameshift
    
    @uniqcoords   = ();
    @mutant_start = ();

    if ($seq eq $nt_reference_seq) { # if the query is same as reference
        $mutations{$seq}->[$index]->{"Mutstart"} = "NA";
        $mutations{$seq}->[$index]->{"Mutstop"}  = "NA";
        $mutations{$seq}->[$index]->{"Mutlength"}= "NA";
        $mutations{$seq}->[$index]->{"Muttype"}  = "NA";
        $mutations{$seq}->[$index]->{"Mutnumber"}= "NA";
        $aa_output{$seq}->{"tMS"}  = 0;
        $aa_output{$seq}->{"tNS"}  = 0;
        $aa_output{$seq}->{"tSYN"} = 0;
        $aa_output{$seq}->{"tAM"}  = 0;
        $aa_output{$seq}->{"tIN"}  = 0;
        $aa_output{$seq}->{"tDEL"} = 0;
        $aa_output{$seq}->{"tFS"}  = 0;
        $aa_output{$seq}->{"tmut"} = 0;
        next;
    }
    else {
        foreach $i (0..$#all_r_codons) { #i is codon position
            if ($all_r_codons[$i] ne $all_q_codons[$i]) {
                #find positions that differ

                @q_codon = split(//, $all_q_codons[$i]);
                @r_codon = split(//, $all_r_codons[$i]);

                # first classify aa mutation
                # Whatif reference is ! or *?
                $q_aa = aa($all_q_codons[$i]);
                $r_aa = aa($all_r_codons[$i]);
#               print "Ref aa: $r_aa\t query aa: $q_aa\n";
                if ($all_r_codons[$i] =~ m/[N]/ ||
                    $all_q_codons[$i] =~ m/[N]/) {
                    $AMcoords{$i} = 1; #Ambiguous
                }
                elsif ($q_aa eq "?" || $r_aa eq "?") {  
                    $AMcoords{$i} = 1;
                } #past this point, assume all aa are not ambiguous
                elsif ($q_aa eq $r_aa) {
                    $SYNcoords{$i} = 1;
                } # past this point, assumed q ne r
                elsif ($q_aa eq "*") {
                    $NScoords{$i} = 1;
                }
                elsif ($q_aa eq "!") {
                    $FScoords{$i} = 1;
                }
                elsif ($q_aa ne "-" && $r_aa ne "-" &&
                    $q_aa ne "*") {
                    $MScoords{$i} = 1;
                }
                elsif ($q_aa eq "-" && $r_aa ne "-") {
                    $DELcoords{$i} = 1;
                }
                elsif ($q_aa ne "-" && $r_aa eq "-") {
                    $INcoords{$i} = 1;
                }
                else { print "Unknown aa case at $i: ref is $r_aa and query is $q_aa\nseqids: $nt_uniqseqs{$seq}\n" };
            }
        }
    }

    %MScontig = ();
    if (%MScoords) {
        @uniqcoords = sort {$a <=> $b } keys %MScoords;
        $mut_start = (); #codon index
        foreach $i (0..$#uniqcoords) {
#        print "uniqcoords is $i\n";
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $codon_start;
#                print "mutstart is $mut_start\n";
            }
#            print "i is $i, end of uniqcoords is $#uniqcoords\n";
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $codon_start;
                    $MScontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $codon_start;
                $MScontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }

    %NScontig = ();
    if (%NScoords) {
        @uniqcoords = sort {$a <=> $b } keys %NScoords;
        $mut_start = (); #codon index
        foreach $i (0..$#uniqcoords) {
#        print "uniqcoords is $i\n";
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $codon_start;
#                print "mutstart is $mut_start\n";
            }
#            print "i is $i, end of uniqcoords is $#uniqcoords\n";
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $codon_start;
                    $NScontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $codon_start;
                $NScontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }
    
    %SYNcontig = ();
    if (%SYNcoords) {
        @uniqcoords = sort {$a <=> $b } keys %SYNcoords;
        $mut_start = (); #codon index
        foreach $i (0..$#uniqcoords) {
#        print "uniqcoords is $i\n";
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $codon_start;
#                print "mutstart is $mut_start\n";
            }
#            print "i is $i, end of uniqcoords is $#uniqcoords\n";
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $codon_start;
                    $SYNcontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $codon_start;
                $SYNcontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }

    %AMcontig = ();
    if (%AMcoords) {
        @uniqcoords = sort {$a <=> $b } keys %AMcoords;
        $mut_start = (); #codon index
        foreach $i (0..$#uniqcoords) {
#        print "uniqcoords is $i\n";
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $codon_start;
#                print "mutstart is $mut_start\n";
            }
#            print "i is $i, end of uniqcoords is $#uniqcoords\n";
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $codon_start;
                    $AMcontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $codon_start;
                $AMcontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }

    %INcontig = ();
    if (%INcoords) {
        @uniqcoords = sort {$a <=> $b } keys %INcoords;
        $mut_start = (); #codon index
        foreach $i (0..$#uniqcoords) {
#        print "uniqcoords is $i\n";
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $codon_start;
#                print "mutstart is $mut_start\n";
            }
#            print "i is $i, end of uniqcoords is $#uniqcoords\n";
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $codon_start;
                    $INcontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $codon_start;
                $INcontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }
    
    %DELcontig = ();
    if (%DELcoords) {
        @uniqcoords = sort {$a <=> $b } keys %DELcoords;
        $mut_start = (); #codon index
        foreach $i (0..$#uniqcoords) {
#        print "uniqcoords is $i\n";
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $codon_start;
#                print "mutstart is $mut_start\n";
            }
#            print "i is $i, end of uniqcoords is $#uniqcoords\n";
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $codon_start;
                    $DELcontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $codon_start;
                $DELcontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }
    
    %FScontig = ();
    if (%FScoords) {
        @uniqcoords = sort {$a <=> $b } keys %FScoords;
        $mut_start = (); #codon index
        foreach $i (0..$#uniqcoords) {
#        print "uniqcoords is $i\n";
            if (!defined $mut_start) {
                $mut_start = $uniqcoords[$i] + $codon_start;
#                print "mutstart is $mut_start\n";
            }
#            print "i is $i, end of uniqcoords is $#uniqcoords\n";
            if ($i < $#uniqcoords) {
                if ($uniqcoords[$i+1] != $uniqcoords[$i] + 1) {
                    $mut_stop = $uniqcoords[$i] + $codon_start;
                    $FScontig{$mut_start} = $mut_stop;
                    $mut_start = ();
                    $mut_stop = ();
                }
            }
            if ($i == $#uniqcoords) {
                $mut_stop = $uniqcoords[$i] + $codon_start;
                $FScontig{$mut_start} = $mut_stop;
                $mut_start = ();
                $mut_stop = ();
            }
        }
    }

    @mutant_start = ();
    push (@mutant_start, keys(%MScontig));
    push (@mutant_start, keys(%NScontig));
    push (@mutant_start, keys(%SYNcontig));
    push (@mutant_start, keys(%AMcontig));
    push (@mutant_start, keys(%INcontig));
    push (@mutant_start, keys(%DELcontig));
    push (@mutant_start, keys(%FScontig));
    @sorted_starts = sort {$a <=> $b} (@mutant_start);

    foreach $i (@sorted_starts) {
        if (exists $MScontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $MScontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($MScontig{$i} - $i) + 1; #counting # nt changed; single poly counted
            $mutations{$seq}->[$index]->{"Muttype"}  = "MS";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        elsif (exists $NScontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $NScontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($NScontig{$i} - $i) + 1; #counting # nt changed; single poly counted
            $mutations{$seq}->[$index]->{"Muttype"}  = "NS";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        elsif (exists $SYNcontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $SYNcontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($SYNcontig{$i} - $i) + 1; #counting # nt changed; single poly counted
            $mutations{$seq}->[$index]->{"Muttype"}  = "SYN";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        elsif (exists $AMcontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $AMcontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($AMcontig{$i} - $i) + 1;
            $mutations{$seq}->[$index]->{"Muttype"}  = "AM";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        elsif (exists $INcontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $INcontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($INcontig{$i} - $i) + 1;
            $mutations{$seq}->[$index]->{"Muttype"}  = "IN";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        elsif (exists $DELcontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $DELcontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($DELcontig{$i} - $i) + 1;
            $mutations{$seq}->[$index]->{"Muttype"}  = "DEL";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        elsif (exists $FScontig{$i}) {
            $mutations{$seq}->[$index]->{"Mutstart"} = $mutstart = $i;
            $mutations{$seq}->[$index]->{"Mutstop"}  = $FScontig{$i};
            $mutations{$seq}->[$index]->{"Mutlength"}= ($FScontig{$i} - $i) + 1;
            $mutations{$seq}->[$index]->{"Muttype"}  = "FS";
            $mutations{$seq}->[$index]->{"Mutnumber"}= $mutnumber;
            $mutnumber++;
            $index++;
        }
        else {
            print "couldn't re-trace mutstart $i\n";
        }
    }

    $_ = 0 for ($tMS, $tNS, $tSYN, $tAM, $tIN, $tDEL, $tFS);

    foreach $i (0..$#{$mutations{$seq}}) {
        if ($mutations{$seq}->[$i]->{"Muttype"} eq "MS") {
            $tMS += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "NS") {
            $tNS += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "SYN") {
            $tSYN += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "AM") {
            $tAM += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "IN") {
            $tIN += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "DEL") {
            $tDEL += $mutations{$seq}->[$i]->{"Mutlength"};
        }
        elsif ($mutations{$seq}->[$i]->{"Muttype"} eq "FS") {
            $tFS += $mutations{$seq}->[$i]->{"Mutlength"};
        }
    }

    $aa_output{$seq}->{"tMS"} = $tMS;
    $aa_output{$seq}->{"tNS"} = $tNS;
    $aa_output{$seq}->{"tSYN"}= $tSYN;
    $aa_output{$seq}->{"tAM"} = $tAM;
    $aa_output{$seq}->{"tIN"} = $tIN;
    $aa_output{$seq}->{"tDEL"}= $tDEL;
    $aa_output{$seq}->{"tFS"} = $tFS;
    $aa_output{$seq}->{"tmut"}= $tMS + $tNS + $tSYN + $tAM + $tIN + $tDEL + $tFS;

}


                          ###########################
                          # Generating Output Table #
                          ###########################

print "\nfinal aa output generated\n";

open F, '>', $aa_outfile or die "problem saving output to file\n";
    print F "allele\tfreq\ttmut\ttMS\ttNS\ttSYN\ttAM\ttIN\ttDEL\ttFS\tmutnum\ttype\tmstart\tmstop\tmlength\tstart\tbstart\tstop\tbstop\tseqs\n";

foreach $key (sort {$aa_output{$a}->{"Freq"} <=> $aa_output{$b}->{"Freq"} ||
                    $aa_output{$a}->{"tmut"} <=> $aa_output{$b}->{"tmut"}} keys %aa_output)  {
    $summary_1 =$allele_ID . "\t" .
                $aa_output{$key}->{"Freq"}. "\t" .
                $aa_output{$key}->{"tmut"}. "\t" .
                $aa_output{$key}->{"tMS"} . "\t" .
                $aa_output{$key}->{"tNS"} . "\t" .
                $aa_output{$key}->{"tSYN"}. "\t" .
                $aa_output{$key}->{"tAM"} . "\t" .
                $aa_output{$key}->{"tIN"} . "\t" .
                $aa_output{$key}->{"tDEL"}. "\t" .
                $aa_output{$key}->{"tFS"} . "\t";

    $summary_3 =$aa_output{$key}->{"Start"}  . "\t" .
                $aa_output{$key}->{"Bstart"} . "\t" .
                $aa_output{$key}->{"Stop"}   . "\t" .
                $aa_output{$key}->{"Bstop"}  . "\t" .
                $aa_output{$key}->{"SeqIDs"};

    for ($i = 0; $i <= $#{$mutations{$key}}; $i++) {
        $summary_2 = $mutations{$key}->[$i]->{"Mutnumber"} . "\t" .
                     $mutations{$key}->[$i]->{"Muttype"}   . "\t" .
                     $mutations{$key}->[$i]->{"Mutstart"}  . "\t" .
                     $mutations{$key}->[$i]->{"Mutstop"}   . "\t" .
                     $mutations{$key}->[$i]->{"Mutlength"} . "\t";
        print F $summary_1 . $summary_2 . $summary_3 . "\n";
    }
    $allele_ID++;
}
close F;


print "Output printed. Remember that start/stop are index 1, bstart/bstop are index 0, and coordinates of mutation are index 0 as well\n";
##########################################
sub print_usage {
    my ($error) = @_; 

    if (defined $error) {
        print STDERR $error, "\n";
    }   

    print "\nUsage: $0 -query [ALIGN.FASTA] -ref [REF.FASTA] -start [START] -stop [STOP] -buffer [BUFFER] -greedy [TRUE]\n";
    print "\tquery:  aligned amino acid files for analysis\n";
    print "\tref:    fasta file for comparison against input\n";
    print "\tstart:  the beginning location (aa) for gene of interest (index 1)\n";
    print "\tstop:   the end location (aa) for gene of interest (index 1)\n";
    print "\tbuffer: distance around gene of interest (aa) for inclusive analysis\n";
    print "\tgreedy: will take entire sequence length, no need for start/stop to be specified\n";
    print "\tnote that all coordinates are indexed by 0 in output.\n";
    print "\nCheers!\n\n";
    exit; 
}

sub aa {
  my ($d) = @_;
  $d = uc $d;
  if ($d =~ /!/)   { return "!"; }
  if ($d eq "TTT") { return "F"; }
  if ($d eq "TTC") { return "F"; }
  if ($d eq "TTY") { return "F"; }
  if ($d eq "TTA") { return "L"; }
  if ($d eq "TTG") { return "L"; }
  if ($d eq "TTR") { return "L"; }
  if ($d eq "TCT") { return "S"; }
  if ($d eq "TCC") { return "S"; }
  if ($d eq "TCA") { return "S"; }
  if ($d eq "TCG") { return "S"; }
  if ($d =~ /^TC/) { return "S"; }
  if ($d eq "TAT") { return "Y"; }
  if ($d eq "TAC") { return "Y"; }
  if ($d eq "TAY") { return "Y"; }
  if ($d eq "TAA") { return "*"; }
  if ($d eq "TAG") { return "*"; }
  if ($d eq "TGT") { return "C"; }
  if ($d eq "TGC") { return "C"; }
  if ($d eq "TGY") { return "C"; }
  if ($d eq "TGA") { return "*"; }
  if ($d eq "TGG") { return "W"; }
  if ($d eq "CTT") { return "L"; }
  if ($d eq "CTC") { return "L"; }
  if ($d eq "CTA") { return "L"; }
  if ($d eq "CTG") { return "L"; }
  if ($d =~ /^CT/) { return "L"; }
  if ($d eq "YTA") { return "L"; }
  if ($d eq "YTG") { return "L"; }
  if ($d eq "YTR") { return "L"; }
  if ($d eq "CCT") { return "P"; }
  if ($d eq "CCC") { return "P"; }
  if ($d eq "CCA") { return "P"; }
  if ($d eq "CCG") { return "P"; }
  if ($d =~ /^CC/) { return "P"; }
  if ($d eq "CAT") { return "H"; }
  if ($d eq "CAC") { return "H"; }
  if ($d eq "CAY") { return "H"; }
  if ($d eq "CAA") { return "Q"; }
  if ($d eq "CAG") { return "Q"; }
  if ($d eq "CAR") { return "Q"; }
  if ($d eq "CGT") { return "R"; }
  if ($d eq "CGC") { return "R"; }
  if ($d eq "CGA") { return "R"; }
  if ($d eq "CGG") { return "R"; }
  if ($d =~ /^CG/) { return "R"; }
  if ($d eq "MGA") { return "R"; }
  if ($d eq "MGG") { return "R"; }
  if ($d eq "MGR") { return "R"; }
  if ($d eq "ATT") { return "I"; }
  if ($d eq "ATC") { return "I"; }
  if ($d eq "ATA") { return "I"; }
  if ($d eq "ATY") { return "I"; }
  if ($d eq "ATW") { return "I"; }
  if ($d eq "ATM") { return "I"; }
  if ($d eq "ATH") { return "I"; }
  if ($d eq "ATG") { return "M"; }
  if ($d eq "ACT") { return "T"; }
  if ($d eq "ACC") { return "T"; }
  if ($d eq "ACA") { return "T"; }
  if ($d eq "ACG") { return "T"; }
  if ($d =~ /^AC/) { return "T"; }
  if ($d eq "AAT") { return "N"; }
  if ($d eq "AAC") { return "N"; }
  if ($d eq "AAY") { return "N"; }
  if ($d eq "AAA") { return "K"; }
  if ($d eq "AAG") { return "K"; }
  if ($d eq "AAR") { return "K"; }
  if ($d eq "AGT") { return "S"; }
  if ($d eq "AGC") { return "S"; }
  if ($d eq "AGY") { return "S"; }
  if ($d eq "AGA") { return "R"; }
  if ($d eq "AGG") { return "R"; }
  if ($d eq "AGR") { return "R"; }
  if ($d eq "GTT") { return "V"; }
  if ($d eq "GTC") { return "V"; }
  if ($d eq "GTA") { return "V"; }
  if ($d eq "GTG") { return "V"; }
  if ($d =~ /^GT/) { return "V"; }
  if ($d eq "GCT") { return "A"; }
  if ($d eq "GCC") { return "A"; }
  if ($d eq "GCA") { return "A"; }
  if ($d eq "GCG") { return "A"; }
  if ($d =~ /^GC/) { return "A"; }
  if ($d eq "GAT") { return "D"; }
  if ($d eq "GAC") { return "D"; }
  if ($d eq "GAY") { return "D"; }
  if ($d eq "GAA") { return "E"; }
  if ($d eq "GAG") { return "E"; }
  if ($d eq "GAR") { return "E"; }
  if ($d eq "GGT") { return "G"; }
  if ($d eq "GGC") { return "G"; }
  if ($d eq "GGA") { return "G"; }
  if ($d eq "GGG") { return "G"; }
  if ($d =~ /^GG/) { return "G"; }
  if ($d eq "---") { return "-"; }
  return "?";
}
