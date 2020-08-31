#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;


# Get command line options
my %opt = ();
GetOptions( \%opt, 'vcfs=s', 'tumor-id=s', 'normal-id=s');
##my @vcfs = check_options( \%opt );
my @vcfs = split/,/,$opt{vcfs};




my $vcf = vcf2->new('file'=>$vcfs[0]);
### which caller ###
my $vc = which_variantcaller( $vcf->{meta} );
print $vc."\n";
while ( my $var_t = $vcf->next_var() ) {

    my $ct = $var_t->{CHROM}; # tumor chrom
    my $st = $var_t->{POS}; # start tumor
    my $et = $var_t->{INFO}->{POS}; #end tumor

    my $vcf2 = vcf2->new('file'=>$vcfs[1]);
    while ( my $var_n = $vcf2->next_var() ) {
        
        my $cn = $var_n->{CHROM}; # normal chrom
        my $sn = $var_n->{POS}; # start normal
        my $en = $var_n->{INFO}->{POS}; #end normal
        my $gtn = $var_n->{GT}->[0]->{GT};
        if ($ct eq $cn) {
            if ($vc eq 'melt' ) {
                if (  abs($st - $sn) <= 5 ){
                    print $var_t->{vcf_str}."\t$gtn";
                }
                else {
                    print $var_t->{vcf_str}."\t0/0\n";
                }
            }
            else {
               # if ($st )
                print $var_t->{vcf_str}."\t$gtn";
            }
        }
        else {
            next;
        }
            

    }

}



    






sub which_variantcaller{
    my $meta = shift;
    if( $meta->{source} ) {
	return "melt" if $meta->{source} =~ /MELT/;
	return "manta" if $meta->{source} =~ /GenerateSVCandidates/;
	return "cnvkit" if $meta->{source} =~ /CNVkit/;
    }

    return "unknown";
}

