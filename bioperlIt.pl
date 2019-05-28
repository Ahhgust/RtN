#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;



my $seqio_object = Bio::SeqIO->new(-file => "humanFasta.withReference.fa");
my $out =          Bio::SeqIO->new(-file => ">humans.fa", format=>'FASTA');

while (my $seq = $seqio_object->next_seq() ) {
    $out->write_seq($seq);
}

$out->close();

