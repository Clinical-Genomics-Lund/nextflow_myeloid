#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

# TODO:
#  * Order sample names in header according to sample-order
#  * Check that all sample-order sample names exist
#  * Check order of sample names if no sample-order given
#  * Remove INFO-fields and old GT-data from both header and variants
#  * Add pindel support


my @supported_callers = ('freebayes', 'mutect2', 'tnscope', 'vardict', 'pindel' );

# Get command line options
my %opt = ();
GetOptions( \%opt, 'vcfs=s', 'tumor-id=s', 'normal-id=s', 'fluffify-pindel', 'sample-order=s', 'paired=s' );
my @vcfs = check_options( \%opt );

my @sample_order;
@sample_order = split /,/, $opt{'sample-order'} if $opt{'sample-order'};

# Aggregate the vcfs
my( $agg_vars, $agg_header, $filters ) = aggregate_vcfs( @vcfs );

# Output final vcf
print_header($filters, $vcfs[0]);
foreach my $vid (keys %$agg_vars ) {
    vcfstr($agg_vars->{$vid}, \@sample_order);
}


sub print_header {
    my $filters = shift;
    my $file = shift;
    
    print "##fileformat=VCFv4.2\n";
    print "##origin=".join(",", @vcfs)."\n";
    print '##INFO=<ID=variant_callers,Number=.,Type=String,Description="List of variant callers which detected the variant">'."\n";
    print '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'."\n";
    print '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">'."\n";
    print '##FORMAT=<ID=VAF,Number=1,Type=Integer,Description="ALT allele observation fraction">'."\n";
    print '##FORMAT=<ID=VD,Number=1,Type=Integer,Description="ALT allele observation count">'."\n";
    print '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'."\n";
    print '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles; If unknown, will be -1">'."\n";
    print '##INFO=<ID=END,Number=.,Type=Integer,Description="END position set to start position for insertions">'."\n";
    foreach( @$filters ) {
	print "##FILTER=<ID=$_,Description=\"$_\">\n";
    }
    if( @sample_order ) {
	print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
	print join("\t", @sample_order)."\n";
    }
    else {
	system("zgrep ^#CHROM $file");
    }
	    
}


sub aggregate_vcfs {
    my @vcfs = @_;

    my @agg;
    my %agg;
    my @headers;

    my %all_filters;
    my %filters;
    
    foreach my $fn ( @vcfs ) {

	my $vcf = vcf2->new('file'=>$fn );

	push @headers, $vcf->{meta};
	
	my $vc = which_variantcaller( $vcf->{meta} );

	
	while ( my $var = $vcf->next_var() ) {
	    my $simple_id = $var->{CHROM}."_".$var->{POS}."_".$var->{REF}."_".$var->{ALT};
        my $type = $var->{INFO}->{SVTYPE};
        my $end = $var->{INFO}->{END};
        my $len = $var->{INFO}->{SVLEN};
		if ($vc eq 'MELT' || $vc eq 'Delly') {$len = $end - $var->{POS} -1; }
	    $var->{INFO} = {};
	    $var->{INFO_order} = [];
	    
	    # Collect all filters for each variant
	    if( $var->{FILTER} ) {
		my @vc_filters = split /;/, $var->{FILTER};
		foreach( @vc_filters ) {
		    $filters{$simple_id}->{$_} = 1;
		    $all_filters{$_} = 1;
		}
	    }
        next if $type eq 'BND';
        add_info($var, "SVTYPE", $type);
        add_info($var, "SVLEN", $len);
        add_info($var, "END", $end);

	    next if $vc eq "freebayes" and is_weird_freebayes($var);
	    
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

    foreach my $id (keys %agg) {
	$agg{$id}->{FILTER} = summarize_filters( keys %{$filters{$id}} )
    }
    
    return( \%agg, \@headers, [keys %all_filters] );
}


sub is_weird_freebayes {
    my $var = shift;
    return 1 unless defined $var->{GT}->[0]->{AD};
    return 0;
}

sub summarize_filters {
    my @filters = @_;

    my @non_pass;
    my $pass = 1;
    foreach (@filters ) {
	push( @non_pass, $_ ) if $_ ne "PASS" and $_ ne ".";
	$pass = 0 if $_ =~ /FAIL/;
    }
    unshift @non_pass, "PASS" if $pass;
    
    if( @non_pass ) {
	return join ";", @non_pass;
    }
    else {
	return "PASS";
    }	
}

sub fix_gt {
    my( $var, $vc ) = @_;

    $var->{FORMAT} = [];
    
    if( $vc =~ /^(manta)$/ ) {
        for my $gt ( @{$var->{GT}} ) {
            my ($ref_dp,$alt_dp) = split /,/, $gt->{PR};
            my $af = 0;
            unless ($alt_dp == 0) {
                $af = sprintf "%.3f", $alt_dp / ($alt_dp + $ref_dp);
            }
            my $gt_s = "0/1";
            if ($af == 0) { $gt_s = "0/0";}
            add_gt( $var, $gt->{_sample_id}, "GT", $gt_s);
            add_gt( $var, $gt->{_sample_id}, "VAF", $af );
            add_gt( $var, $gt->{_sample_id}, "VD", $alt_dp );
            add_gt( $var, $gt->{_sample_id}, "DP", $alt_dp+$ref_dp );
        }
    }

    elsif( $vc eq "CNVkit" ) {
		if ($opt{paired} eq 'paired') {
			my %tmp = ( 'GT' => '0/0', 'CN'=> '2', '_sample_id'=> 'dummy');
			my @tmp = split /\./, $var->{GT}->[0]->{_sample_id};
			
			if ( $tmp[1] eq $opt{'tumor-id'} ) {
				$tmp{_sample_id} = $opt{'normal-id'}; 
				$var->{GT}->[0]->{_sample_id} = $opt{'tumor-id'};
			}
			else {
				$tmp{_sample_id} = $opt{'tumor-id'};
				$var->{GT}->[0]->{_sample_id} = $opt{'normal-id'};
			}
			push @{$var->{GT}},\%tmp;
		}
        
		for my $gt ( @{$var->{GT}} ) {
            my $cnum = "0.5";
            if ($gt->{CN}) { $cnum = $gt->{CN}; }
			add_gt( $var, $gt->{_sample_id}, "GT", $gt->{GT});
			add_gt( $var, $gt->{_sample_id}, "VAF", $cnum );
			add_gt( $var, $gt->{_sample_id}, "VD", "0");
			add_gt( $var, $gt->{_sample_id}, "DP", "0");

		}
    }
	elsif( $vc eq "MELT") {
		if ($opt{paired} eq 'paired') {
			my %tmp = ( 'GT' => '0/0');
			if ($var->{GT}->[0]->{_sample_id} eq $opt{'tumor-id'} ) {$tmp{_sample_id} = $opt{'normal-id'}; }
			else {$tmp{_sample_id} = $opt{'tumor-id'}; }
			push @{$var->{GT}},\%tmp;
		}
		for my $gt ( @{$var->{GT}} ) {

			add_gt( $var, $gt->{_sample_id}, "GT", $gt->{GT});
			add_gt( $var, $gt->{_sample_id}, "VAF", "0" );
			add_gt( $var, $gt->{_sample_id}, "VD", "0");
			add_gt( $var, $gt->{_sample_id}, "DP", "0");

		}
	}
    elsif( $vc eq "Delly") {
		
		for my $gt ( @{$var->{GT}} ) {
            my $af = 0;
            my $altc = 0;
            my $depth = 0;
            if ($gt->{RV} == 0) {
                $af = sprintf "%.3f", $gt->{DV} / ($gt->{DR}+$gt->{DV});
                $altc = $gt->{DV};
                $depth = ($gt->{DV}+$gt->{DR});
            }
            else {
                $af = sprintf "%.3f", $gt->{RV} / ($gt->{RV}+$gt->{RR});
                $altc = $gt->{RV};
                $depth = ($gt->{RV}+$gt->{RR});
            }
            
			add_gt( $var, $gt->{_sample_id}, "GT", $gt->{GT});
			add_gt( $var, $gt->{_sample_id}, "VAF", $af );
			add_gt( $var, $gt->{_sample_id}, "VD", $altc );
			add_gt( $var, $gt->{_sample_id}, "DP", $depth );

		}
	}
}

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
    
    print "   --sample-order    Comma separated list of order of samples names in output VCF.\n";
    print "                     If not given, all input vcfs must have same order.\n\n";
    print "   --fluffify-pindel Modify pindel REF/ALT fields to be less than 1000 bp (comply with Manta).\n";
    print "\n";
    exit(0);
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
	return "CNVkit" if $meta->{source} =~ /CNVkit/;
	return "manta" if $meta->{source} =~ /GenerateSVCandidates/;
	return "MELT" if $meta->{source} =~ /MELT/;
    #return "Delly" if $meta->{source} =~ //;
    }
    else {
        return "Delly" if $meta->{FILTER}->{PASS}->{Description} =~ /All filters passed/;
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
    my( $v, $sample_order ) = @_;

    my @all_info;
    print $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

    # Generate and print INFO field
    for my $info_key (@{$v->{INFO_order}}) {
	push @all_info, $info_key."=".$v->{INFO}->{$info_key};
    }
    print join(";", @all_info)."\t";

    # Print FORMAT field
    print join(":", @{$v->{FORMAT}})."\t";

    
    my %order;
    my $i=0;
    if( @$sample_order > 0 ) {
	$order{$_} = $i++ foreach @{$sample_order};
    }
    else {
	$order{$_->{_sample_id}} = $i++ foreach @{$v->{GT}};
    }

    # Print GT fields for all samples
    for my $gt ( sort {$order{$a->{_sample_id}} <=> $order{$b->{_sample_id}}} @{$v->{GT}}) {
	my @all_gt;
	for my $key ( @{$v->{FORMAT}} ) {
	    push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
	}
	print join(":", @all_gt)."\t";
    }
    print "\n";
}
