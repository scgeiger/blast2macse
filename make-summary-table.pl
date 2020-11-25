#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Text::CSV;

my ($DATA_ID, $PARENTDIR, $PARENTPATH, $REFPATH, $TIMESTAMP, $GENE);
my $db_nseqs;
my $hit_nseqs;
my $ref_length;
my $aln_length;
my $dataset;
my ($ntTS, $ntTV, $ntIN, $ntDEL, $ntAM, $ntFS, $aaMS, $aaNS, $aaSYN, $aaIN, $aaDEL, $aaFS, $aaAM); 
my ($ntHap, $ntSingHap, $nSegSite, $nSingSNPs, $n_ntSingHap, $n_ntMajorHap);
my ($PG_TajimaD, $PG_RozasR2, $PG_FuliF, $PG_FuliD, $PG_FuFs, $PG_pi, $PG_SegSite, $PG_Haplo, $PG_SingHaps, $PG_nSeqs); 
my ($PG_TajimaD_R, $PG_RozasR2_R, $PG_FuliF_R, $PG_FuliD_R, $PG_FuFs_R, $PG_pi_R, $PG_SegSite_R, $PG_Haplo_R, $PG_SingHaps_R, $PG_nSeqs_R);
my ($seq, $seqID, $row, $nRem);
my $outfile;
my $output;
my $all_file;
my $header;
my @all_genes = ();
my %blast_db = ();
my @temp = ();
my @SUBDIR = ();
my $DIR;
my $command;
my $temp;
my $CSV;
my %hash;
my $pos;
my $no_consensus_err_file;
my $file;

GetOptions (
    'refdir=s'     => \$REFPATH,     # Directory with reference sequences
);

# Verify input
if (!defined $REFPATH) {
    &print_usage("\nPlease specify the directory path holding your reference genes");
}

# Use this code in the ST directory, not in the analysis directory
# Assumes my/path/is/analysis/ST/genes

$TIMESTAMP = getLoggingTime();
$PARENTPATH = cwd();
    print "parent path is $PARENTPATH\n";
($PARENTDIR) = ($PARENTPATH =~ m/([^\/]+)$/);
    print "parent dir is $PARENTDIR\n";

opendir ($DIR, $REFPATH) || die "Error in opening dir $REFPATH\n";
    while ( ($row = readdir($DIR))) {
        if ($row eq "." || $row eq "..") {
            next;
        }
        $row =~ s/\..*//;
        push(@all_genes, $row); #all gene names are saved in array
    }         
closedir $DIR;

opendir ($DIR, $PARENTPATH) || die "Error in opening dir $PARENTPATH\n";
    while ( $row = readdir($DIR)) {
        next if $row =~ /^\.\.?$/ || ! -d $row;
        push(@SUBDIR, $row);
    }
closedir($DIR);

# Making header
$header      =  "dataset"     . "\t" .
                "gene"        . "\t" .
                "db_nseqs"    . "\t" .
                "hit_nseqs"   . "\t" .
                "ref_length"  . "\t" .
                "aln_length"  . "\t" .
                "ntTS"        . "\t" .
                "ntTV"        . "\t" .
                "ntIN"        . "\t" .
                "ntDEL"       . "\t" .
                "ntFS"        . "\t" .
                "ntAM"        . "\t" .
                "aaMS"        . "\t" .
                "aaNS"        . "\t" .
                "aaSYN"       . "\t" .
                "aaIN"        . "\t" .
                "aaDEL"       . "\t" .
                "aaFS"        . "\t" .
                "aaAM"        . "\t" .
                "ntHap"       . "\t" .
                "ntSingHap"   . "\t" .
                "nSegSite"    . "\t" .
                "nSingSNPs"   . "\t" .
                "n_ntSingHap" . "\t" .
                "n_ntMajorHap". "\t" .
                "PG_TajimaD"  . "\t" .
                "PG_RozasR2"  . "\t" .
                "PG_FuliF"    . "\t" .
                "PG_FuliD"    . "\t" .
                "PG_FuFs"     . "\t" .
                "PG_pi"       . "\t" .
                "PG_SegSite"  . "\t" .
                "PG_Haplo"    . "\t" .
                "PG_SingHaps" . "\t" .
                "PG_nSeqs"    . "\t" .
                "nRem"        . "\t" .
                "PG_TajimaD_R". "\t" .
                "PG_RozasR2_R". "\t" .
                "PG_FuliF_R"  . "\t" .
                "PG_FuliD_R"  . "\t" .
                "PG_FuFs_R"   . "\t" .
                "PG_pi_R"     . "\t" .
                "PG_SegSite_R". "\t" .
                "PG_Haplo_R"  . "\t" .
                "PG_SingHaps_R". "\t" .
                "PG_nSeqs_R"  . "\n";

$all_file =  "$PARENTPATH/" . "$PARENTDIR" . "-all-summary.tsv";
open S, '>', $all_file or die "problem saving output to all_file\n";
    print S "#Analysis run on $TIMESTAMP\n";
    print S $header;
 
#$no_consensus_err_file = "no-consensus.err"; 
foreach $GENE (@all_genes) {
    $outfile = "$PARENTPATH/$GENE/" . "$PARENTDIR-$GENE" . "-summary.tsv";
    print "$GENE\n";
#    if ( $GENE =~ /\(/ ) {
#        $GENE =~ s/([\\()])/\\$1/g;
#    }
#    print "$GENE\n";
    $DATA_ID = $PARENTDIR . "-" . $GENE;

    ###############################
    # Check if raw-blast is empty #
    ###############################
    if ( -z "$PARENTPATH/$GENE/$DATA_ID.raw-blast" ) {
    open O, '>', $outfile or die "problem saving output to file\n";
        print O "#Analysis run on $TIMESTAMP, using gene $GENE.\n";
        print "$GENE: Didn't find any BLAST hits. Exiting table with defaults\n";
        $output =  $dataset  . "\t" .
                $GENE         . "\t" .
                $db_nseqs     . "\t" .
                "0"    . "\t" . #hit_nseqs
                "NA"   . "\t" . #ref_length
                $aln_length   . "\t" . #aln_length
                "0"         . "\t" . #ntTS
                "0"         . "\t" . #ntTV
                 "0"   . "\t" . #ntIN
                  "0"  . "\t" . #ntDEL
                 "0"   . "\t" . #ntFS
                 "0"   . "\t" . #ntAM
                 "0"   . "\t" . #aaMS
                 "0"   . "\t" . #aaNS
                 "0"   . "\t" . #aaSYN
                 "0"   . "\t" . #aaIN
                 "0"   . "\t" . #aaDEL
                 "0"   . "\t" . #aaDEL
                 "0"   . "\t" .
                 "0"    . "\t" .
                "NA"   . "\t" .
                "NA"   . "\t" .
                "NA"  . "\t" .
                "NA"  . "\t" .
                "NA" . "\t" .
                "NA"   . "\t" .
                "NA"   . "\t" .
                 "NA"   . "\t" .
                 "NA"   . "\t" .
                 "NA"   . "\t" .
                 "NA"       . "\t" .
                "NA"  . "\t" .
                 "NA". "\t" .
                "NA" . "\t" .
                 "NA"   . "\t" .
                 "NA"       . "\t" .
                "NA" . "\t" .
                "NA" . "\t" .
                "NA"  . "\t" .
                "NA"  . "\t" .
                 "NA"  . "\t" .
                 "NA"  . "\t" .
                "NA".  "\t" .
                "NA"  . "\t" .
                "NA" . "\t" .
                "NA"  . "\n";

            print O $output;
            print S $output;
            #}
            close O;
            next;
    }

    ########################
    # Get reference length #
    ########################
    $temp = (); 
    open F, '<', "$REFPATH/$GENE.nt";
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^>/ || $row =~ /^#/) {
                next;
            }
            else {
                $temp .= $row;
            }
        }
    close F; 
    $ref_length = length($temp);
    #print "ref_length for $GENE is $ref_length\n";

    open O, '>', $outfile or die "problem saving output to file\n";
        print O "#Analysis run on $TIMESTAMP, using gene $GENE.\n";
        print O $header;
        $DIR = ();
        $DIR = $PARENTDIR;
#        foreach $DIR (@SUBDIR) {
#            $DATA_ID = $DIR . "-" . $GENE;
            $DATA_ID = $PARENTDIR . "-" . $GENE;
#            $dataset = "$DIR";
             $dataset = "$PARENTDIR";
            ######################################
            # get number of sequences downloaded #
            ######################################
            @temp = ();
#            print "llama $PARENTPATH/$PARENTDIR-count.txt\n";
            open my $F, '<', "$PARENTPATH/$PARENTDIR-count.txt";
                chomp (@temp = <$F>);
                $db_nseqs = $temp[0];
            close $F;
           
            ############################
            # get number of blast hits #
            ############################
            @temp = ();
            $command = "wc -l ". "$PARENTPATH/$GENE/$DATA_ID.raw-blast";
            if ( $command =~ /\(/ || $command =~ /\)/ || $command =~ /\'/) {
                $command =~ s/([\\()\'])/\\$1/g;
            }
            $command .= " | awk '{ print \$1 }'";
            print "command is $command\n";
            @temp = `$command`;
            chomp $temp[0];
            $hit_nseqs = $temp[0];
#            $hit_nseqs =~ s/\s.*//; #remove everything after number
            #print "hit_nseqs is $hit_nseqs\n"; 
      
            ######################## 
            # get reference length #
            ########################
            # Felt cute might delete later
            #print "ref_length is $ref_length\n";
           
            ######################## 
            # get alignment length #
            ########################
            # Getting length of consensus
            # At this point, have already screened out blast files that are empty
            # If the consensus file does not exist, there is an error, most likely via MACSE alignment
            # Need to exit loop here

            $temp = ();
            open F, '<', "$PARENTPATH/$GENE/$DATA_ID.consensus";
                while ($row = <F>) {
                    chomp $row;
                    if ($row =~ /^>/ || $row =~ /^#/) {
                        next;
                    }
                    else {
                        $temp .= $row;
                    }
                }
            close F; 
            $aln_length = length($temp);
            #print "aln_length is $aln_length\n";

            ###################### 
            # get codon mut data #
            ######################
            @temp = ();
            %hash = ();
            open F, '<', "$PARENTPATH/$GENE/$DATA_ID.macse-codon-dist.tsv";
                while ($row = <F>) {
                    chomp $row;
                    if ($row =~ /^#Date/) { next; } ;
                    if ($row =~ /^#NT/) {
                        $row =~ s/^[^||]*\|\|\s//;
                        @temp = split /\,/, $row;
                        foreach $temp (@temp) {
                            $temp =~ s/^\s+//; #remove leading spaces
                            if ($temp =~ /^Ambiguous/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $ntAM = $temp; 
                                #print "ntAM is $ntAM\n";
                            }
                            if ($temp =~ /^Transition/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $ntTS = $temp;
                                #print "ntTS is $ntTS\n";
                            }
                            if ($temp =~ /^Transversion/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $ntTV = $temp;
                                #print "ntTV is $ntTV\n";
                            }
                            if ($temp =~ /^Frameshift/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $ntFS = $temp;
                                #print "ntFS is $ntFS\n";
                            }
                            if ($temp =~ /^Insertion/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $ntIN = $temp;
                                #print "ntIN is $ntIN\n";
                            }
                            if ($temp =~ /^Deletion/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//; #remove leading spaces
                                $ntDEL = $temp;
                                #print "ntDEL is $ntDEL\n";
                            }
                        }
                    }
                    elsif ($row =~ /^#AA/) {
                        $row =~ s/^[^||]*\|\|\s//;
                        @temp = split /\,/, $row;
                        foreach $temp (@temp) {
                            $temp =~ s/^\s+//; #remove leading spaces
                            if ($temp =~ /^Ambiguous/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $aaAM = $temp;
                                #print "aaAM is $aaAM\n";
                            }
                            if ($temp =~ /^Missense/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $aaMS = $temp;
                                #print "aaMS is $aaMS\n";
                            }
                            if ($temp =~ /^Nonsense/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $aaNS = $temp;
                                #print "aaNS is $aaNS\n";
                            }
                            if ($temp =~ /^Synonymous/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $aaSYN = $temp;
                                #print "aaSYN is $aaSYN\n";
                            }
                            if ($temp =~ /^Frameshift/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $aaFS = $temp;
                                #print "aaFS is $aaFS\n";
                            }
                            if ($temp =~ /^Insertion/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//;
                                $aaIN = $temp;
                                #print "aaIN is $aaIN\n";
                            }
                            if ($temp =~ /^Deletion/) {
                                $temp =~ s/^[^:]*\://;
                                $temp =~ s/^\s+//; #remove leading spaces
                                $aaDEL = $temp;
                                #print "aaDEL is $aaDEL\n";
                            }
                        }
                    }
                    else {
                        @temp = split /\t/, $row;
                        if ($temp[3] eq "Transversion" || $temp[3] eq "Transition") {
                            $pos = $temp[0];
                            $hash{$pos} = 1;
                        }
                    }
                }
                $nSegSite = scalar(keys(%hash));
            close F;

            ####################### 
            # get uniqallele data #
            #######################
            $temp = ();
            @temp = ();
#            my $csv = Text::CSV->new({sep => "\t"});
            my @uniq = ();
            $nSingSNPs = 0;
            $n_ntSingHap = 0;
            $n_ntMajorHap = 0; 
            open F, '<', "$PARENTPATH/$GENE/$DATA_ID.macse-nt-uniq.tsv"; # or die "Could not open uniqallele data\n";
                while ($row = <F>) {
                    chomp $row;
                    push @uniq, $row;
                }
            close F;
            # Getting haplotype numbers
            #$ntHap [last position of column 0] #index 0 
            @temp = split /\t/, $uniq[$#uniq];
            $ntHap = $temp[0];
            #print "ntHap is $ntHap\n";
            %hash = (); 
            #n_ntMajorHap
            #[in final row, column 18 (last column), how many space delmitied items exist?]
            $n_ntMajorHap = $temp[$#temp] =~ tr{ }{ };
            #print "n_ntMajorHap is $n_ntMajorHap\n";
        
            # Iterating through array to get other values
            for $row (1..$#uniq) { #for each row index in uniq output
                @temp = split /\t/, $uniq[$row]; # for each column in 
                $temp = $temp[$#temp] =~ tr{ }{ };
                if ($temp == 1) {
                    $hash{$temp[0]} = 1;
                }
            }
            $n_ntSingHap = $ntSingHap = scalar(keys %hash); 
            #print "ntSingHap is $ntSingHap\n";  
            #print "n_ntSingHap is $n_ntSingHap\n";
            
            ######################
            # get PopGenome data #
            ######################
            # Getting PopGenome Data
            @temp = ();
            open F, '<', "$PARENTPATH/$GENE/$DATA_ID-PopGenome.tsv";
                while ($row = <F>) {
                    chomp $row;
                    push @temp, $row; 
                }
            close F;
            my @temp2 = ();
            @temp2 = split /\t/, $temp[1];
            $PG_TajimaD  = $temp2[0];
            $PG_RozasR2  = $temp2[2];
            $PG_FuliF    = $temp2[3];
            $PG_FuliD    = $temp2[4];
            $PG_FuFs     = $temp2[5];
            $PG_pi       = $temp2[6];
            $PG_SegSite  = $temp2[1];
            $PG_Haplo    = $temp2[8];
            $PG_SingHaps = $temp2[9];
            $PG_nSeqs    = $temp2[7];
            
            ######################
            # Get rm PG data     #
            ######################
            # Getting PopGenome Data with removed seqs
            @temp = ();
            open F, '<', "$PARENTPATH/$GENE/$DATA_ID-removed-macse-PopGenome.tsv";
                while ($row = <F>) {
                    chomp $row;
                    push @temp, $row;
                }
            close F;
            @temp2 = ();
            @temp2 = split /\t/, $temp[1];
            $PG_TajimaD_R  = $temp2[0];
            $PG_RozasR2_R  = $temp2[2];
            $PG_FuliF_R    = $temp2[3];
            $PG_FuliD_R    = $temp2[4];
            $PG_FuFs_R     = $temp2[5];
            $PG_pi_R       = $temp2[6];
            $PG_SegSite_R  = $temp2[1];
            $PG_Haplo_R    = $temp2[8];
            $PG_SingHaps_R = $temp2[9];
            $PG_nSeqs_R    = $temp2[7];

            ######################
            # Get nRem seq       #
            ######################
            # Getting number of seqs removed because of dels or Ns 
            @temp = ();
#            print "path is" . "$PARENTPATH/$GENE/$DATA_ID-removed-macse-ids.tsv\n";
            open F, '<', "$PARENTPATH/$GENE/$DATA_ID-removed-macse-ids.tsv";
                @temp = <F>;
                $temp = $temp[0];
                chomp $temp;
            close F;
            @temp = ();
#            print "temp is $temp\n";
            @temp = split (';', $temp);
#            print "$temp[2]\n";
            $nRem = $temp[2];
#            print "1 $nRem\n";
            $nRem =~ s/^\s+//;
#            print "2 $nRem\n";
            $nRem =~ s/[^ ]* //;
#            print "nrem is $nRem\n";

 
        ###################
        # Printing output #
        ###################
     $output =  $dataset      . "\t" .
                $GENE         . "\t" .
                $db_nseqs     . "\t" .
                $hit_nseqs    . "\t" .
                $ref_length   . "\t" .
                $aln_length   . "\t" .
                $ntTS         . "\t" .
                $ntTV         . "\t" .
                $ntIN         . "\t" .
                $ntDEL        . "\t" .
                $ntFS         . "\t" .
                $ntAM         . "\t" .
                $aaMS         . "\t" .
                $aaNS         . "\t" .
                $aaSYN        . "\t" .
                $aaIN         . "\t" .
                $aaDEL        . "\t" .
                $aaFS         . "\t" .
                $aaAM         . "\t" .
                $ntHap        . "\t" .
                $ntSingHap    . "\t" .
                $nSegSite     . "\t" .
                $nSingSNPs    . "\t" .
                $n_ntSingHap  . "\t" .
                $n_ntMajorHap . "\t" .
                $PG_TajimaD   . "\t" .
                $PG_RozasR2   . "\t" .
                $PG_FuliF     . "\t" .
                $PG_FuliD     . "\t" .
                $PG_FuFs      . "\t" .
                $PG_pi        . "\t" .
                $PG_SegSite   . "\t" .
                $PG_Haplo     . "\t" .
                $PG_SingHaps  . "\t" .
                $PG_nSeqs     . "\t" .
                $nRem        . "\t" .
                $PG_TajimaD_R. "\t" .
                $PG_RozasR2_R. "\t" .
                $PG_FuliF_R  . "\t" .
                $PG_FuliD_R  . "\t" .
                $PG_FuFs_R   . "\t" .
                $PG_pi_R     . "\t" .
                $PG_SegSite_R. "\t" .
                $PG_Haplo_R  . "\t" .
                $PG_SingHaps_R. "\t" .
                $PG_nSeqs_R  . "\n";

    print O $output;
    print S $output;
    #}
    close O;
} 

close S;

print "Summary tables generated\n";
# Remember to track indexes
# Add descriptive header

sub print_usage {
    my ($error) = @_;

    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 -refdir [directory with reference genes]\n";
    print "\tThe script is highly customized\n";

    print "\tBE SURE TO START THIS IN YOUR PARENT DIRECTORY\n";
    print "\trefdir: directory holding reference sequences that were blasted\n";
    print "\tNote that query and ref must both be the same length\n";

    print "\nCheers!\n\n";
    exit;
}

sub getLoggingTime {

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d%02d%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}
