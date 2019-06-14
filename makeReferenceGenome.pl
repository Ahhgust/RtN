#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::Seq;


my $seqio_object = Bio::SeqIO->new(-file => "chrM.circularized.fa", 
                                   -format=>'FASTA');
my $out =          Bio::SeqIO->new(-file => ">humans.fa", -format=>'FASTA');

# reference sequence; on top!
while (my $seq = $seqio_object->next_seq() ) {
    
    my $s = $seq->seq();
    $s =~ tr/acgt/ACGT/;
    $s =~ s/[^GATC]//g;
    my $newSeq = Bio::Seq->new(
        -seq => $s,
        -id => $seq->display_id() );
    
    $out->write_seq($newSeq);
}

$seqio_object = Bio::SeqIO->new(-file => "alg_tot.fasta", -format=>'FASTA');

while (my $seq = $seqio_object->next_seq() ) {
    my $s = $seq->seq();
    $s =~ tr/acgt/ACGT/;
    $s =~ s/[^GATC]//g;
    my $newSeq = Bio::Seq->new(
        -seq => $s . $s,
        -id => $seq->display_id() );
    
    $out->write_seq($newSeq);
}

$out->close();

