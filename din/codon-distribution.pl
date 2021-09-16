#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $ref;
my $query;
my $nSeqs = 0;
my $r_codon;
my $q_codon;
my @all_q_codons = ();
my @all_r_codons = ();
my $length;
my $row;
my $seqID;
my $refID;
my $refseq;
my @r_array = ();
my @q_array = ();
my $pos;
my @output = ();
my %query = ();
my $out_row;
my $outfile = "macse-codon-distribution.tsv";
my $nt_type, 
my ($q_nt, $r_nt, $q_aa, $r_aa);
my ($i, $j);
my @q_codon = ();
my @r_codon = ();
my $mut_type;
my $codon_index = 0;
my $date;
my $output;
my $aa_type;
my $nSynonymous = 0;
my $nNonsense = 0;
my $nFrameshift_aa = 0;
my $nMissense = 0;
my $nDeletion_aa = 0;
my $nInsertion_aa = 0;
my $nTransition = 0;
my $nTransversion = 0;
my $nInsertion_nt = 0;
my $nFrameshift_nt = 0;
my $nDeletion_nt = 0;
my $nAmbiguous_nt = 0;
my $nAmbiguous_aa = 0;
my $index = 0;      #gives row number
my $codonpos;
my $outline;
my $codon;
my $seqlength;

GetOptions (
    'ref=s'   => \$ref,
    'query=s' => \$query,
);

# Verify input
if (!defined $query || !-f $query) {
    &print_usage("\nPlease specify an aligned nucleic acid file for analysis.");
}

if (!defined $ref || !-f $ref) {
    &print_usage("\nPlease specify a reference file.");
}

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
                $seqlength = length($refseq)-1;
            }
        }
    close F;
}  else { print "\n Error: Ref file not reading correctly.\n"; }

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
                $query{$seqID} = $row;
                if (length($row) != length($refseq)) {
                    # test this v
                    die "Error: Sequences are not the same length as ref\n";
                }
            }
        }
    close F;
} else { print "\n Error: Query file not reading correctly.\n"; }

print "\n All sequences collected successfully\n";

# method 1
@all_r_codons = ( $refseq =~ m/...?/g );


# test print all r codons
#print "\t\t### Start ###\n";
#print "Below are all ref codons\n";
#foreach my $c (0..$#all_r_codons) {
#    print "$c: $all_r_codons[$c]\n";
#}
#print "\n\n";

foreach $seqID (keys (%query)) {
    @all_q_codons = ($query{$seqID} =~ m/...?/g);
    $pos = 0;
    foreach $i (0..$#all_r_codons) {
        if ($all_r_codons[$i] ne $all_q_codons[$i]) {
            #find positions that differ
            @q_codon = split(//, $all_q_codons[$i]);
            @r_codon = split(//, $all_r_codons[$i]);

            # first classify aa mutation
            # Whatif reference is ! or *?
            $q_aa = aa($all_q_codons[$i]);
            $r_aa = aa($all_r_codons[$i]);
            if ($all_r_codons[$i] =~ m/[N]/ ||
                $all_q_codons[$i] =~ m/[N]/) {
                $aa_type = "Ambiguous";
                $nAmbiguous_aa++;
            }
            elsif ($q_aa eq "?" || $r_aa eq "?") {
                $aa_type = "Ambiguous";
                $nAmbiguous_aa++;
            } #past this point, assume all aa are not loqal
            elsif ($q_aa eq $r_aa) {
                $aa_type = "Synonymous";
                $nSynonymous++;
            } # past this point, assumed q ne r
            elsif ($q_aa eq "*") {
                # assumed nonsense mutation
                $aa_type = "Nonsense";
                $nNonsense++;
            }
            elsif ($q_aa eq "!") {
                $aa_type = "Frameshift";
                $nFrameshift_aa++;  
            }
            elsif ($q_aa ne "-" && $r_aa ne "-" &&
                   $q_aa ne "*") {
                $aa_type = "Missense";
                $nMissense++;
            }
            elsif ($q_aa eq "-" && $r_aa ne "-") {
                $aa_type = "Deletion";
                $nDeletion_aa++;
            }
            elsif ($q_aa ne "-" && $r_aa eq "-") {
                $aa_type = "Insertion";
                $nInsertion_aa++;
            }

            foreach $j (0..$#r_codon) { #values will be 0, 1, or 2;
                if ($r_codon[$j] ne $q_codon[$j]) {
#                    $codon_index = $j;
                    $q_nt = $q_codon[$j];
                    $r_nt = $r_codon[$j];
                    $pos = (3 * $i) + $j;  #3 nt / codon * $i codon  
                    $codonpos = $j;
    
                    # Classifying mutational types; nt
                    if ($q_nt =~ m/[N]/ ||
                        $r_nt =~ m/[N]/) {
                        $nt_type = "Ambiguous";
                        $nAmbiguous_nt++;
                    }
                    elsif ($q_nt ne "-" && $r_nt ne "-" && 
                           $q_nt ne "!" && $r_nt ne "!" ) {
                        #assumend to be SNP of some flavor
                        if (($q_nt eq "A" && $r_nt eq "G") ||
                            ($q_nt eq "G" && $r_nt eq "A") ||
                            ($q_nt eq "C" && $r_nt eq "T") ||
                            ($q_nt eq "T" && $r_nt eq "C")) {
                                $nt_type = "Transition";
                                $nTransition++;
                        }
                        if (($q_nt eq "A" && $r_nt eq "C") ||
                            ($q_nt eq "C" && $r_nt eq "A") ||
                            ($q_nt eq "A" && $r_nt eq "T") ||
                            ($q_nt eq "T" && $r_nt eq "A") ||
                            ($q_nt eq "G" && $r_nt eq "T") ||
                            ($q_nt eq "T" && $r_nt eq "G") ||
                            ($q_nt eq "G" && $r_nt eq "C") ||
                            ($q_nt eq "C" && $r_nt eq "G")) {
                                $nt_type = "Transversion";
                                $nTransversion++;
                        }
                    } 
                    elsif ($q_nt ne "-" && $q_nt ne "!" && $r_nt eq "-") {
                        #assumed to be insertion
                        $nt_type = "Insertion";
                        $nInsertion_nt++;
                    }
                    elsif ($q_nt eq "-" && $r_nt ne "-") {
                        #assumed to be a deletion
                        $nt_type = "Deletion";
                        $nDeletion_nt++;
                    }
                    elsif ($q_nt eq "!") {
                        #assumed to be frameshift
                        $nt_type = "Frameshift";
                        $nFrameshift_nt++;
                    }
                    else {
                        $nt_type = "Ambiguous";
                        $nAmbiguous_nt++;
                    }
                    # Done classifying nt mutation
                    # Can assign values to output array here
                     
                    $outline =
                               $pos      . "\t" .
                               $q_nt     . "\t" .
                               $r_nt     . "\t" .
                               $nt_type  . "\t" .
                               $codonpos . "\t" .
                               $i        . "\t" .
                               $all_q_codons[$i]    . "\t" .
                               $all_r_codons[$i]    . "\t" .
                               $q_aa     . "\t" .
                               $r_aa     . "\t" .
                               $aa_type  . "\t" .
                               $seqlength. "\t" .
                               $seqID    . "\t" .
                               $refID    . "\t" .
                               $nSeqs    ;
                    push (@output, $outline);
                }
            }
        }
    }
}

# Printing output
$date = getLoggingTime();
open O, '>', $outfile or die "problem saving output to file\n";
    print O "#Date: $date Query: $query Reference: $ref\n";
    print O "#AA Total Count|| Ambiguous: $nAmbiguous_aa, Synonymous: $nSynonymous, Nonsense: $nNonsense, Frameshift: $nFrameshift_aa, Missense: $nMissense, Insertion: $nInsertion_aa, Deletion: $nDeletion_aa\n"; 
    print O "#NT Total Count|| Ambiguous: $nAmbiguous_nt, Transition: $nTransition, Transversion: $nTransversion, Frameshift: $nFrameshift_nt, Insertion: $nInsertion_nt, Deletion: $nDeletion_nt \n";
    print O "pos\tq_nt\tr_nt\tnt_type\tcodonpos\tcodon_index\tq_codon\tr_codon\tq_aa\tr_aa\taa_type\tlength\tseqid\trefid\tnseqs\n";
    foreach $row (0..$#output) {
        print O "$output[$row]\n";
    }
close O;

    # pos is the position in the alignment of the nt permutation being discussed
    # qnt is the query nt at the site of interest
    # r_nt is the reference nt at the site of interest
    # mut_type is the nucleotide mutation classification
    # codon is the 3nt codon at this site
    # codon_pos is the location in the codon (position 0,1 or 2)
    # q_aa is the translated amino acid fed by the q_codon
    # r_aa is the translated amino acid as fed by the r_codon
    # mutttype is the class of amino acid mutation observed here
    # seqID is the sequence ID being assessed
    # refid is teh reference ID being assessed

print "\nFinal output generated as $outfile\n";

##########################################
sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 -query [ALIGN.FASTA] -ref [REF.FASTA]\n";
    print "\tThe script is being tailored to accomodate MACSE nt output\n";

    print "\tquery:  aligned nucleic acid files for analysis\n";
    print "\tref:    fasta file for comparison against input\n";
    print "\tNote that query and ref must both be the same length\n";

    print "\tNote that all coordinates are indexed by 0 in output.\n";
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

sub getLoggingTime {

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d%02d%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}
