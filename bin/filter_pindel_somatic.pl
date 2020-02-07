#!/usr/bin/perl -w
use strict;

# Removes all variants in a VCF that has the genotype 0/0 for all samples.

die "Usage: filter_pindel_somatic.pl IN_VCF OUT_VCF\n" unless $ARGV[0] and $ARGV[1];
die "$ARGV[0] does not exist or is empty!" unless -s $ARGV[0];

open( UNFILT, $ARGV[0] );
open( FILT, '>'.$ARGV[1] ) or die "Cannot write $ARGV[1]";
while( <UNFILT> ) {
    if( $_ =~ /^#/ ) {
	print FILT $_;
    }
    else {
	my @a = split /\t/;
	my $keep = 0;
	foreach( 9..$#a ) {
	    $keep = 1 if $a[$_] !~ /^0\/0/;
	}
	print FILT $_ if $keep;
    }
}
close UNFILT;
close FILT;
