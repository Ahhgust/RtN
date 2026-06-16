#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;

my $inFile = "humanFasta.withReference.fa";
my $outFile = "humans.fa";

GetOptions('i=s' => \$inFile,
           'o=s' => \$outFile);


if (-f $outFile && -s $outFile) {
    die "$outFile exists. Not gonna clobber it!\n";
}

open(FOO, $inFile);
my %seqs;
my $id;
while (<FOO>) {

    $_ =~ s/\s+//g;
    if (substr($_, 0, 1) eq '>') {
        
        $id = substr($_, 1);
        $id =~ s/\s+//g;
    } else {
        my $s = $_;
        $s =~ s/[^GATC]//g;
        if ($s ne '') {
            $seqs{$id} .= $_;
        }
    }
}


my $out =          Bio::SeqIO->new(-file => ">$outFile", format=>'FASTA');

#while (my $seq = $seqio_object->next_seq() ) {
foreach (keys %seqs) {
   my $seq = Bio::Seq->new( -display_id => $_,
                             -seq => $seqs{$_});
   $out->write_seq($seq);

}

$out->close();

