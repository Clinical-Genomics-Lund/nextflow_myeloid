#!/usr/bin/perl -w
use strict;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use Data::Dumper;
use List::Util qw( max min );

my $vcf = vcf2->new('file'=>$ARGV[0] );

#TODO ADD INFO HEADER STRINGS PROPERLY!
system("zgrep ^## $ARGV[0]");
system("zgrep ^#CHROM $ARGV[0]");

while ( my $v = $vcf->next_var() ) {

    my @filters;


    my( %likelihood, %gl_idx, %genotype, %altobs, %depth );
    my $status = "PASS";
    for my $gt (@{$v->{GT}}) {
        my $type = "T";

	# Fail if GT is 0/0 for tumor
	$status = "FAIL_GT" if $type eq "T" and $gt->{GT} eq "0/0";
	
	my @GL = split /,/, $gt->{GL};
	my @GT = split /\//, $gt->{GT};
	my @AO = split /,/, ($gt->{AO} or "0");

	my $DP = $gt->{DP};
	my $RO = $gt->{RO};

	$depth{$type}      = $DP;
	$altobs{$type}     = \@AO;
	$genotype{$type}   = \@GT;
	$likelihood{$type} = \@GL;
    }


    my $TALT = $genotype{T}->[0];
    if( $TALT eq "0" ) {
	$status = "WARN_NOVAR";
    } 

    $v->{FILTER} = $status;

    vcfstr($v);
}





sub add_info {
    my( $var, $key, $val ) = @_;
    push( @{$var->{INFO_order}}, $key );
    $var->{INFO}->{$key} = $val;
}

sub vcfstr {
    my $v = shift;

    my @all_info;
    print $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

    # Generate and print INFO field
    for my $info_key (@{$v->{INFO_order}}) {
	push @all_info, $info_key."=".$v->{INFO}->{$info_key};
    }
    print join(";", @all_info)."\t";

    # Print FORMAT field
    print join(":", @{$v->{FORMAT}})."\t";

    # Print GT fields for all samples
    for my $gt (@{$v->{GT}}) {
	my @all_gt;
	for my $key ( @{$v->{FORMAT}} ) {
	    push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
	}
	print join(":", @all_gt)."\t";
    }
    print "\n";
}
