#!/usr/bin/env perl

use strict;
use warnings;

if (@ARGV < 2) {
    die "$0 rCRS_hg38_bam rCRS_bam\n";
}

my $minMapq = 4;

open IN, $ARGV[0] or
    die "Failed to open $ARGV[0] for reading\n";

# BAM file mapped to  the RCRS+hg38
my @s;
my %readsInOne = ();
my $id = '';
while (<IN>) {
    next if  m/^@/; #skip the header
    @s = split /\t/;
    if ($s[2] eq 'chrM' && $s[4] >= $minMapq) { # maps cleanly to the mito
        $id = $s[0];
        $id =~ s/^[^:]+//;
        $readsInOne{$id}=1;
    }
}

close IN;


open IN, $ARGV[1] or
    die "Failed to open $ARGV[1] for reading\n";

# BAM mapped JUST to  rCRS
while (<IN>) {
    if ( m/^@/) {
        print $_;
    } else{
        @s = split /\t/;
        $id = $s[0];
        $id =~ s/^[^:]+//;
        if (!exists $readsInOne{$id}) { # argued to be an incorrect mapping.
            print $_;
        }
    }
    
}

close IN;



