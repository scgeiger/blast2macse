#!/usr/bin/perl
use warnings;
use strict;

# Last edited 210827. Cleaned and added exit 1
# Assumes input sequence file with seqid\tseq\n
# Will put output as fasta format with >seqid\nseq\n

# Declare variables
my ($infile, $filename, $outfile);
my ($row, $seqID);
my ($check_in, $check_out);
my @temp = ();
my %seq_check = ();
my %out_seqs = ();

# Verify input
if (!defined $ARGV[0]) {
    &print_usage("\nPlease specify a sequence file.");
}

$infile = $ARGV[0];
($filename) = ($infile =~ /([^.]+)/);
$outfile = "formatted-" . $filename . ".blast";
$check_in = 0;

# Get file contents
if (-f $infile) {
    open F, '<', $infile;
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                next;
            }
            @temp = split(/\t/, $row);
            $seqID  = ">" . $temp[0];
            if (exists($seq_check{$seqID})) {
                $seq_check{$seqID}++;
                $seqID .= ".$seq_check{$seqID}";
            }
            else {
                $seq_check{$seqID} = 0;
            } 
            $out_seqs{$seqID} = $temp[1];
            $check_in++;
        }
    close F;
} else {&print_usage("Did you specify a file with aligned seqs?")}

# Check that number of seqs in equals number of seqs out
$check_out = scalar(keys %out_seqs);
print "\nFormatting has finished. Sanity check:\n";
print "\tInput: $check_in sequences.\n";
print "\tOutput: $check_out sequences.\n\n";

if ($check_in != $check_out) {
    print "\t!!!!Something isn't quite right!!!!\n";
    return undef;
}

# Print Outfile
open O, '>', $outfile;
    for $seqID (keys %out_seqs) {
        if (defined $out_seqs{$seqID}) {
            print O $seqID . "\n";
            print O $out_seqs{$seqID} . "\n";
        }
        else {
            print "\t$seqID was detected without a sequence.\n";
            print "\tThere may be a formatting problem\n";
            print "\tDeleting output file and returning error\n\n";
            close O;
            unlink($outfile) or die "Couldn't clean up problem file\n";
            exit 1;
        }
    }
close O;

print "\nOutfile is $outfile\n\n";

############################
sub print_usage {
    my ($error) = @_;
    if (defined $error) {
        print STDERR $error, "\n";
    }

    print "\nUsage: $0 [sequence file]\n";
    print "\tAssumes input sequence file with seqid\\tseq\\n\n";
    print "\tThis is blastn outfmt 6 sseqid sseq\n";
    print "\tOutput will be fasta-formatted >seqid\\nseq\\n\n";
    print "\tCheers!\n\n";
    exit 1;
}
