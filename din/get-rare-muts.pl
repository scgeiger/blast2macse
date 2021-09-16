#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# This code will be used to identify rare alleles and rare mutations in 
# those alleles. Previous scripts can identify rare alleles, sure, but
# a major failing is that they can't identify which of the SNPs are rare
# as well, or compare these to a threshold. 

# The inputs to this code will be 1) a sequence alignment (aa) and a ref-
# erence sequence (aa). They will have to be the same length. An optional
# option will be a flag (the cutoff) for rare or no. 

# Output will have seqid/refAA/truepos/qaa/hapindex/hapn/happrop
# /mutn/mutprop/cutoff/rarehap/raremut
# where hapn is the number of seqs with haplotype, and happrop is the pro-
# portion of all haplotyes (hapN/allN) in the population. The mut folllows
# a similar rule.

# This assumes that each sequence is arranged so that >ID is on line0, and
# all of its sequence is on line 1. 

# This code will _not_ identify synonymous mutations. 

###########################################################################

# Init Variables
my ($ref, $query, $cutoff, $mut, $ID, $timestamp);
my ($row, $refID, $refseq, $seqlength, $seqID, $nSeqs, $seq, $i);
my ($alleleProp, $alleleFreq, $rareHap, $rareHapOnly);
my ($mutFreq, $mutProp, $rareMut, $rareMutOnly);
my $outfile;
my %allele_seq = ();
my %all_muts = ();
my %query = ();
my @reference = ();
my @query = ();
my ($alnPos, $truPos);
my @qseq = ();
my @seqID = ();
my @ref = ();

# Get Options
GetOptions (
    'r=s' => \$ref,
    'q=s' => \$query,
    'co=f' => \$cutoff,
    'rh=s' => \$rareHapOnly,
    'rm=s' => \$rareMutOnly,
);

# Verify input
if (!defined $query || !-f $query) {
    &print_usage("\nPlease specify an aligned amino acid file for analysis.");
} else {
    $outfile = $query;
    $outfile =~ s{^.*/|\.[^.]+$}{}g;
    $outfile .= "-mut-pos.tsv";
}

if (!defined $ref || !-f $ref) {
    &print_usage("\nPlease specify a reference sequence (aa ).");
}

if (!defined $cutoff) {
    print "No cutoff was specified. Default is 0.001, or 0.01% of the population\n";
    $cutoff = 0.001;
}

if (!defined $rareHapOnly && !defined $rareMutOnly) {
    print "Neither rareHapOnly or rareMutOnly were specified. All mutations will be printed\n";
    print "please flag -rh = TRUE or -rm = TRUE if wanted\n";
    $rareHapOnly = "FALSE";
    $rareMutOnly = "FALSE";
}
if (!defined $rareHapOnly && (defined $rareMutOnly && $rareMutOnly eq "TRUE")) {
    $rareHapOnly = "TRUE";
}
if (defined $rareHapOnly && $rareHapOnly ne "TRUE") {
    print "Rare haps only was specified, but value was not equal to TRUE. Flag will be ignored.\n";
}
if ((defined $rareHapOnly && $rareHapOnly eq "TRUE") && (!defined $rareMutOnly || $rareMutOnly ne "TRUE")) {
    print "Rare haps only was specified, but rare mut only is not TRUE. Will be false.\n";
    $rareMutOnly = "FALSE";
}
if (defined $rareMutOnly && $rareMutOnly ne "TRUE") {
    print "Rare muts only was specified, but value was not equal to TRUE. Flag will be ignored.\n";
}
if ($rareHapOnly eq "TRUE" && $rareMutOnly ne "TRUE") {
    print "Rare haps only was specified, but rareMutOnly was not. Latter will be ignored.\n";
}
if ($rareMutOnly eq "TRUE" && ($rareHapOnly ne "TRUE" || !defined $rareHapOnly)) {
    print "Rare muts only was specified. By default, RareHapsOnly will be TRUE\n";
    $rareHapOnly = "TRUE";
}

$timestamp = getLoggingTime();

# Open reference file and sort out refseq
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
                $refseq = $row;
                $seqlength = length($refseq)-1; #index 0? Check this
            }
        }
    close F;
} else { print "\n Error: Ref file not reading correctly.\n"; }

@reference = split //, $refseq;

# Open query file and put into hash
if (-f $query) {
    open F, '<', $query;
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            if ($row =~ s/^>//) {
                $seqID = $row;
                $seqID =~ s/^\>//d;
                $nSeqs++;  #index 1
            }
            else {
                if (length($row) != length($refseq)) {
                    # test this v
                    die "Error: Sequences are not the same length as ref\n";
                }
                # %allele_seq is my main datastructure. Simplifies analysis into alleles
                $allele_seq{$row}->{"IDs"} .= "$seqID "; 
                $allele_seq{$row}->{"Freq"}++;
            }
        }
    close F;
} else { print "\n Error: Query file not reading correctly.\n"; }
print "\n All sequences collected successfully\n";

# Collect mutation data for each allele
for $seq (keys %allele_seq) {
    if ($seq eq $refseq) {
        # fill in the table variables
        next;
    }
    if (length($refseq) != length($seq)) {
        print "Query does not match $refID length\n";
        print "Seqids: " . $allele_seq{$seq}->{"IDs"}. "\n";
        next;
    }

    @query = split //, $seq;
    $truPos = 0;  #will be index 1
    foreach $i (0..$#query) { #i is index 0
        if ($reference[$i] eq "-" && $query[$i] eq "-") {
            next; # Another sequence has an insertion here 
        } 
        elsif ($reference[$i] eq $query[$i]) {
            $truPos++;
            next; # Both sequences are identical here
        }
        else {
            $truPos++;
            $mut = $reference[$i] . $truPos . $query[$i];
            $all_muts{$mut} += $allele_seq{$seq}->{"Freq"}; # Track frequency of each mutation
            push(@{$allele_seq{$seq}->{"Muts"}}, $mut);
            $allele_seq{$seq}->{"nMuts"}++;
        }
    }
}

# So by now, we should know where each mutation is located in each allele. 
# Time to synthesize
# Two outputs, rare only and all muts? Rare only can be under a switch. 
# open output
open O, '>', $outfile or die "problem saving output to file";
    print O "#$timestamp; QueryFile: $query; RefFile: $ref; RefID: $refID; Cutoff: $cutoff; nSeqs: $nSeqs";
        if ($rareMutOnly eq "TRUE") {
            print O " RareMutOnly: TRUE\n";
        }
        elsif ($rareHapOnly eq "TRUE") {
            print O " RareHapOnly: TRUE\n"
        } else { print O "\n";}
    print O "seqID\tmut\thapN\thapProp\tmutN\tmutProp\trareHap\trareMut\n";
    
    foreach $seq (sort{$allele_seq{$a}->{"Freq"} <=> $allele_seq{$b}->{"Freq"} ||  
                   $allele_seq{$a}->{"nMuts"} <=> $allele_seq{$b}->{"nMuts"}} keys %allele_seq) { 
        $alleleFreq = $allele_seq{$seq}->{"Freq"}; #n times hap is in pop
        $alleleProp = sprintf("%.4f", $alleleFreq/$nSeqs); # What prop of pop has hap?
        if ($alleleProp <= $cutoff) {
            $rareHap = "YES";
        } else { $rareHap = "NO"; }

        if ($rareHapOnly eq "TRUE" && $rareHap eq "NO") {
            next;
        }

        @seqID = split / /, $allele_seq{$seq}->{"IDs"};
        foreach $ID (@seqID) { 
            if (defined $allele_seq{$seq}->{"Muts"}) {
                foreach $mut (@{$allele_seq{$seq}->{"Muts"}}) {
                    $mutFreq = $all_muts{$mut}; # get freq of mut in pop
                    $mutProp = sprintf("%.4f", $mutFreq/$nSeqs); # get prop of mut in pop
                    if ($mutProp <= $cutoff) {
                        $rareMut = "YES";
                    } else { $rareMut = "NO"; }
                    if ($rareMutOnly eq "TRUE" && $rareMut eq "NO") {
                        next;
                    }
                print O "$ID\t$mut\t$alleleFreq\t$alleleProp\t$mutFreq\t$mutProp\t$rareHap\t$rareMut\n";
                }
            } else {
                next; # no mutations were found relative to ref
            }
        }
    }
close O;

print "All mutations have been printed to file\n";


##############################################################################
sub getLoggingTime {

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d%02d%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}

sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }
    print "\nUsage: $0 -q [ALIGN.FASTA] -r [REF.FASTA] (-co INTEGER) (-rh TRUE) (-rm TRUE)\n";
    print "\tThe script works with single-line alignment format\n";
    print "\tYou can use aa alignments only, make sure ref is aligned with query\n";

    print "\tq:  aligned aa file for analysis\n";
    print "\tr:  sequence comparison against input\n";
    print "\tNote that query and ref must both be the same length\n";
    print "\tco: optional flag specifying cutoff for rareness. 0.001 (0.01%) is default.\n";
    print "\trh/rm: optional flag, will print an extra file only holding rare haplotypes or muts\n\n";
    print "\tNote that all coordinates are indexed by 1 in output.\n";
    print "\nCheers!\n\n";
    exit;
}
