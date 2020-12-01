#!/usr/bin/perl
use strict;
use warnings;

my $row;
my $header;
my $seqid;
my $dealign;
my $outfile = $ARGV[0] . ".dealign";
my $errfile = $ARGV[0] . ".e-dealign";
my %all_seqs = ();
my $infile = $ARGV[0];
my @split = ();
my $i;
my $command;

if (-f $ARGV[0]) {
    open F, '<', $ARGV[0];
        while ($row = <F>) {
            chomp $row;
            if ($row =~ /^#/) {
                $header = $row;
            }
            if ($row =~ /^>/) {
                $seqid = $row;
            }
            else {
                $all_seqs{$seqid} = $row;
            }
        }
    close F;
}

open ERR, '>', $errfile;
open O, '>', $outfile;

if ($header) {
    print O "$header\n";
}

foreach $seqid (keys %all_seqs) {
    $dealign = $all_seqs{$seqid};
    $dealign =~ tr/\-//d;
    $dealign =~ tr/\!//d;
    if (length($dealign) == 0) {
        #print "sequence $seqid only has gaps. It will be removed\n".
        $seqid =~ s/^>//;
        @split = split / /, $seqid;
        foreach $i (@split) { 
            print ERR "$i\n";
        }
        next;
    }
    else {
        print O "$seqid\n$dealign\n";
    }
}

close O;
close ERR;

if (-z $errfile) {
    $command = "rm -f $errfile";
    system $command;
    print "Cleaning empty error file\n";
} 

print "Job complete\n";
