#!/usr/bin/perl -w
use strict;
use vcf2;
use Data::Dumper;

my $vcf_filemask = $ARGV[0];

my @files = glob($vcf_filemask);
my %all;
foreach my $fn (@files) {

    my $vcf = vcf2->new('file'=>$fn );

    my $vc = which_variantcaller( $vcf->{meta} );
    
    while ( my $var = $vcf->next_var() ) {
	my($vaf, $vd, $dp) = get_stats($var, $vc);
	next if $dp == 0;
	my $varid = $var->{CHROM}."_".$var->{POS}."_".$var->{REF}."_".$var->{ALT};
	push( @{$all{$varid}}, $vaf);
	#print "$varid\t$vaf\t$vd\t$dp\n";
    }
}

my $num_samples = scalar(@files);

#Summarize pon
foreach my $var ( sort keys %all ) {
    my $num_vars = scalar(@{$all{$var}});
    my $num_germ = count_germline($all{$var});
    my( $mean, $std_dev, $values ) = mean_nongermline($all{$var});
    
    next if $num_germ == $num_vars;
    print $var;
    print "\t".$num_vars;
    print "\t".$num_germ;
    print "\t".$num_samples;
    print "\t".$mean."\t".$std_dev;
    print "\t".$values;
    print "\n";
}
#xprint Dumper(\%all);


    
sub count_germline {
    my $f = shift;

    my $cnt = 0;
    foreach my $frq (@$f) {
	$cnt++ if ($frq > 0.35 and $frq < 0.65) or $frq > 0.85;
    }
    return $cnt;
}

sub mean_nongermline {
    my $f = shift;

    my $cnt = 0;
    my @vals;
    foreach my $frq (@$f) {
	$cnt++;
	push @vals, $frq if $frq < 0.35;
    }
    return (average(@vals), std_dev(@vals), join(",",@vals));
}

sub which_variantcaller{
    my $meta = shift;
    if( $meta->{source} ) {
	return "freebayes" if $meta->{source} =~ /freeBayes/;
	return "mutect2" if $meta->{source} =~ /Mutect2/;
	return "pindel" if $meta->{source} =~ /pindel/;
    }
    if( $meta->{'SentieonCommandLine.TNscope'} ) {
	return "tnscope";
    }
    if( $meta->{INFO}->{MSILEN} ) { # FIXME: Terrible way of detecting vardict VCFs...
	    return "vardict";
    }

    return "unknown";
}


# Returns VAF%, alt count, and DP.
sub get_stats {
    my( $var, $vc ) = @_;

    $var->{FORMAT} = [];

    if( $vc =~ /^(mutect2|tnscope|vardict|pindel)$/ ) {
	for my $gt ( @{$var->{GT}} ) {
	    my ($ref_dp, $alt_dp, $af) = (0,0,0);
	    if( $gt->{AD} ) {
		my @a = split /,/, $gt->{AD};
		$ref_dp = ($a[0] or "0");
		$alt_dp = ($a[1] or "0");
		if( $alt_dp + $ref_dp > 0 ) {
		    $af = $alt_dp / ($alt_dp+$ref_dp);
		}

	    }
	    return( ($gt->{AF} or $af), $alt_dp, $alt_dp+$ref_dp );
	}
    }

    elsif( $vc eq "freebayes" ) {
	for my $gt ( @{$var->{GT}} ) {
	    my( $vaf, $vd ) = (0,0);
	    if( $gt->{AO} and $gt->{AO} ne "." ) {
		$vaf = sprintf "%.4f", $gt->{AO} / $gt->{DP};
		$vd = $gt->{AO};
	    }
	    return( $vaf, $vd, $gt->{DP});

	}
    }
}


sub average {
    my (@values) = @_;

    my $count = scalar @values;
    my $total = 0;
    $total += $_ for @values;

    return $count ? $total / $count : 0;
}


sub std_dev {
    my ($average, @values) = @_;

    my $count = scalar @values;
    my $std_dev_sum = 0;
    $std_dev_sum += ($_ - $average) ** 2 for @values;

    return $count ? sqrt($std_dev_sum / $count) : 0;
}
