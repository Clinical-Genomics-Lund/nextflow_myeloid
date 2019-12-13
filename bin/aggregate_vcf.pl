#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

   

my %opt = ();
my @supported_callers = ('freebayes', 'mutect2', 'tnscope' );

GetOptions( \%opt, 'vcfs=s', 'tumor-id=s', 'normal-id=s', 'fluffify-pindel' );

my @vcfs = check_options( \%opt );

my( $agg_vars, $agg_header ) = aggregate_vcfs( @vcfs );

system("zgrep ^# $vcfs[0]");
foreach my $vid (keys %$agg_vars ) {

    vcfstr($agg_vars->{$vid});
}



sub aggregate_vcfs {
    my @vcfs = @_;

    my @agg;
    my %agg;
    my @headers;
    foreach my $fn ( @vcfs ) {

	my $vcf = vcf2->new('file'=>$fn );

	push @headers, $vcf->{meta};
	
	my $vc = which_variantcaller( $vcf->{meta} );

	while ( my $var = $vcf->next_var() ) {
	    my $simple_id = $var->{CHROM}."_".$var->{POS}."_".$var->{REF}."_".$var->{ALT};

	    if( $agg{$simple_id} ) {
		$agg{$simple_id}->{INFO}->{variant_callers} .= "|$vc";
	    }
	    else {
		add_info( $var, "variant_callers", $vc );
		fix_gt( $var, $vc );
		$agg{$simple_id} = $var;
	    }
 	}
    }
    return( \%agg, \@headers );
}

sub aggregate_headers {
    my $headers = shift;

    my @agg_header = @{$headers->[0]};
    
    
}

sub fix_gt {
    my( $var, $vc ) = @_;

    if( $vc eq "mutect2" or $vc eq "tnscope" ) {
	for my $gt ( @{$var->{GT}} ) {
	    my ($ref_dp, $alt_dp) = (0,0);
	    if( $gt->{AD} ) {
		( $ref_dp, $alt_dp ) = split /,/, $gt->{AD};
	    }
	    add_gt( $var, $gt->{_sample_id}, "VAF", $gt->{AF} );
	    add_gt( $var, $gt->{_sample_id}, "VD", $alt_dp );
	}
    }

    elsif( $vc eq "freebayes" ) {
	for my $gt ( @{$var->{GT}} ) {
	    my $vaf = 0;
	    if( $gt->{AO} and $gt->{AO} ne "." ) {
		$vaf = sprintf ":%.3f", $gt->{AO} / $gt->{DP};
		add_gt( $var, $gt->{_sample_id}, "VAF", $vaf );
	    }
	}
    }
}


	# PINDEL
#	elsif( $caller eq "pindel" ) {
#	    $out_vcf_str .= "\t".$full_info->{GT}->{ $translate_names->{$sample_name} }->{GT};
#	    my( $ref_dp, $alt_dp ) = split /,/, $full_info->{GT}->{ $translate_names->{$sample_name} }->{AD};
	    
#	    $out_vcf_str .= ":".($ref_dp+$alt_dp);
#	    $out_vcf_str .= ":".$alt_dp;
#	    if( $alt_dp == 0 ) {
#		$out_vcf_str .= ":0";
#	    }
#	    else {
#		$out_vcf_str .= sprintf ":%.3f", $alt_dp / ($alt_dp + $ref_dp);
#	    }
#	}

	# MANTA
#	elsif( $caller eq "manta" ) {
#	    $out_vcf_str .= "\t0/1";
#	    my( $ref_dp, $alt_dp ) = split /,/, $full_info->{GT}->{ $translate_names->{$sample_name} }->{PR};
#	    $out_vcf_str .= ":".($ref_dp+$alt_dp);
 #           $out_vcf_str .= ":".$alt_dp;
  #          if( $alt_dp == 0 ) {
   #             $out_vcf_str .= ":0";
    #        }
     #       else {
      #          $out_vcf_str .= sprintf ":%.3f", $alt_dp / ($alt_dp + $ref_dp);
#	    }	    
#	}
	

	# GATK HaplotypeCaller
#	elsif( $caller eq "gatkhc" ) {
#	    $out_vcf_str .= "\t".$full_info->{GT}->{ $translate_names->{$sample_name} }->{GT};
#	    $out_vcf_str .= ":".$full_info->{GT}->{ $translate_names->{$sample_name} }->{DP};
#	    my( $ref_dp, $alt_dp ) = split /,/, $full_info->{GT}->{ $translate_names->{$sample_name} }->{AD};
#	    $out_vcf_str .= ":".$alt_dp;
#	    $out_vcf_str .= ":".($alt_dp/$full_info->{GT}->{ $translate_names->{$sample_name} }->{DP});
#	}
	    

	# STRELKA
#	elsif( $caller eq "strelka" ) {

	    # INDELs
#	    if( $full_info->{GT}->{ $translate_names->{$sample_name} }->{TAR} ) {
#		my $ref_count = ( split ',', $full_info->{GT}->{ $translate_names->{$sample_name} }->{TAR} )[0];
#		my $alt_count = ( split ',', $full_info->{GT}->{ $translate_names->{$sample_name} }->{TIR} )[0];
#		my $dp = $full_info->{GT}->{ $translate_names->{$sample_name} }->{DP};
		
#		my $af = 0;
#		if( ($alt_count + $ref_count) > 0 ) {
#		    $af = $alt_count / ( $alt_count + $ref_count );
#		}
#		my $gt = "0/1";
#		$gt = "0/0" if $af < 0.01;
#		$out_vcf_str .= "\t$gt:".$dp.":".$alt_count.":".$af;    
#	    }

	    # SNVs
#	    else {
#		print STDERR Dumper($full_info);
#		my $REF_FIELD = $full_info->{REF}."U";
#		my $ALT_FIELD = $full_info->{ALT}."U";
		#print STDERR "************* ".$full_info->{REF}."\t".$full_info->{ALT}."\t".$translate_names->{$sample_name}."\n";
#		my $ref_count  = (split ',', $full_info->{GT}->{ $translate_names->{$sample_name} }->{$REF_FIELD} )[0];
#		my $alt_count  = (split ',', $full_info->{GT}->{ $translate_names->{$sample_name} }->{$ALT_FIELD} )[0];
#		my $dp = $full_info->{GT}->{ $translate_names->{$sample_name} }->{DP};
#		my $af = $alt_count / ( $alt_count + $ref_count );
#		my $gt = "0/1";
#		$gt = "0/0" if $af < 0.01;		
#		$out_vcf_str .= "\t$gt:".$dp.":".$alt_count.":".$af;    
#	    }
#	}



sub check_options {
    my %opt = %{ $_[0] };

    help_text() unless $opt{vcfs};
    
    my @files;
    foreach my $opt_key ( sort keys %opt ) {

	if( $opt_key eq "vcfs" ) {
	    my @vcfs = split /,/, $opt{$opt_key};
	    foreach( @vcfs ){
		die "VCF does not exist: $_" unless -s $_;
		push @files, $_;
	    }
	}
    }
    return @files;
}


sub help_text {
    my $error = shift;
    
    print "\n\$ aggregate_vcf.pl --vcfs COMMA_SEP_LIST_OF_VCFS [--fluffify-pindel]\n\n";
    print "   --vcfs        Comma separated list of VCF files to aggregate, in priority order\n";
    print "   Supported callers: ". join(", ", @supported_callers)."\n\n";;
	
    print "   --fluffify-pindel Modify pindel REF/ALT fields to be less than 1000 bp (comply with Manta).\n";
    print "\n";
    exit(0);
}

sub add_info_field {
    my $vcf_str = shift;
    my $data = shift;
    
    my @a = split /\t/, $vcf_str;
    $a[7] .= ";".$data;
    
    return join "\t", @a;
}


sub excl_prefix {
    my( $prefix, $name ) = @_;

    return $name if $prefix eq "DO_NOTHING";

    $name =~ s/^$prefix//;
    return $name;
}

sub is_gzipped {
    my $fn = shift;
    
    my $file_str = `file $fn`;
    return 1 if $file_str =~ /gzip compressed/;
    return 0;
}


sub fluffify_pindel_variants {
    my $vars = shift;

    my $MAX_SIZE = 1000;
    
    foreach my $v ( @$vars ) {
	if( length($v->{REF}) >= $MAX_SIZE or length($v->{ALT}) >= $MAX_SIZE ) {
	    $v->{REF} = substr($v->{REF}, 0, 1);
	    $v->{ALT} = "<".$v->{INFO}->{SVTYPE}.">";
	}
    }
}

sub which_variantcaller{
    my $meta = shift;
    if( $meta->{source} ) {
	return "freebayes" if $meta->{source} =~ /freeBayes/;
	return "mutect2" if $meta->{source} =~ /Mutect2/;
    }
    if( $meta->{'SentieonCommandLine.TNscope'} ) {
	return "tnscope";
    }

    return "unknown";
}

sub add_info {
    my( $var, $key, $val ) = @_;
    push( @{$var->{INFO_order}}, $key );
    $var->{INFO}->{$key} = $val;
}

sub add_gt {
    my( $var, $sample, $key, $val ) = @_;
    push( @{$var->{FORMAT}}, $key ) unless grep /^$key/, @{$var->{FORMAT}};
    foreach my $gt ( @{$var->{GT}} ) {
	if( $gt->{_sample_id} eq $sample ) {
	    $gt->{$key} = $val;
	}
    }
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
