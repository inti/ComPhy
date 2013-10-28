#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;
use Getopt::Long;
use Pod::Usage;
use PDL;
use PDL::Matrix;
use PDL::NiceSlice;
use Data::Dumper;

use constant E_CONSTANT => log(10);

# local modules
use PhyloStratiphytUtils;

our (   $help, $man, $out, $score_matrix, $annot, $collapse, $summary_annotation, $bootstrap );

GetOptions(
        'help' => \$help,
        'man' => \$man,
        'out|o=s' => \$out,
        'score_matrix|m=s' => \$score_matrix,
        'annot|a=s' => \$annot,
        'collapse|s=s' => \$collapse,
        'summary_annotation' => \$summary_annotation,
        'bootstrap|b=i' => \$bootstrap,

) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);
pod2usage(-exitstatus => 2, -verbose => 2) if (not defined $annot and not defined $score_matrix);

defined $bootstrap or $bootstrap = 2000;

# read genes annotation
print_OUT("Reading gene's annotation from [ $annot ]");
my %annotation = ();
my %gene_sets = ();
open (ANNOT,"$annot") or die $!;
while(my $line = <ANNOT>){
    chomp($line);
    my ($set,$gene) = split(/\t/,$line);
    push @{ $gene_sets{$set}->{genes} }, $gene;
    push @{ $annotation{$gene}->{sets} }, $set;
}
close(ANNOT);
print_OUT("    ... done ...");

# gene gene's scores
# sum number of genes for each stratum of each gene_set 
print_OUT("Reading gene's phylostratiphraphic scores");

open(MATRIX,"$score_matrix") or die $!;
my @phylostratum = ();
while (my $line = <MATRIX>){
    chomp($line);
    if ($. == 1) {
        @phylostratum = split(/\t/,$line);
        shift @phylostratum if ($phylostratum[0] eq 'ID');
        next;
    }
    my ($gene,@scores) = split(/\t/,$line);
    $annotation{$gene}->{scores} = \@scores;
}

my %gene_index = ();
my $c = 0;
foreach my $gene (sort {$a cmp $b } keys %annotation){
    $gene_index{$gene} = $c++;
}
my $phylostratum_scores = mpdl map {  $annotation{$_}->{scores};  } sort { $gene_index{$a} <=> $gene_index{$b} } keys %annotation;

## define variables for the z-score calculation
my $collapse_idx = -1;
if (defined $collapse){
    my $c = 0;
    foreach my $stratum (@phylostratum){
        $collapse_idx = $c if ($stratum eq $collapse);
        $c++;
    }
    if ($c == -1){
        print_OUT("Could not find [ $collapse ] among the phylostratum. Continuing without executing this option");
        $collapse= undef;
    } else {
        $collapse_idx++;
        print_OUT("Will generate a new phylostratum with levels after [ $collapse ] called [ After_$collapse ]");
    }
}

print_OUT("Calculating enrichment z-scores");
my $phylostratum_totals = $phylostratum_scores->sumover->flat;
my $N = $phylostratum_totals->sum;
my $out_string_z_score = join "\t", ("GO",@phylostratum);
if (defined $collapse){
    $out_string_z_score .= "\tAfter_$collapse";
}
$out_string_z_score .= "\n";
my $out_string_sets = join "\t", ("GO",@phylostratum);
if (defined $collapse){
    $out_string_sets .= "\tAfter_$collapse";
}
$out_string_sets .= "\n";


foreach my $set (keys %gene_sets){
    my $profile = zeroes scalar @phylostratum;
    foreach my $gene (@{$gene_sets{$set}->{genes}}){
        $profile += pdl $annotation{$gene}->{scores};
    }
    my $R = $profile->sum();
    my @z_scores = ();
    for (my $i = 0; $i < $profile->nelem; $i++){
        my ($r) = $profile($i)->list;
        my ($n) = $phylostratum_totals($i)->list;
        if ($n == 0){
            push @z_scores, "NaN";
            next;
        }
        # Values are
        # N=total number of genes on genome,
        # R = total number of genes on the gene-set.
        # r = number of gene-set's genes on a phylostratum and
        # n = total number of genes  on phylostratum of interest.
        #
        push @z_scores, hygeometric_dist_z_score($N,$n,$R,$r);
    }
    if (defined $collapse){
        my ($n) = $phylostratum_totals($collapse_idx:$profile->nelem-1)->flat->sum;
        my ($r) = $profile($collapse_idx:$profile->nelem-1)->flat->sum;
        push @z_scores, hygeometric_dist_z_score($N,$n,$R,$r);
#        print "@phylostratum\n",$profile,"\n",$profile($collapse_idx:$profile->nelem-1)->flat->sum,"\n","N $N n $n R $R r $r z ",hygeometric_dist_z_score($N,$n,$R,$r);;
#        getc;
        
    }

    $out_string_z_score .= "$set\t" . join "\t", @z_scores;
    $out_string_z_score .= "\n";
    $out_string_sets .= "$set\t" . join "\t", $profile->list;
    if (defined $collapse){
        $out_string_sets .= "\t" . $profile($collapse_idx:$profile->nelem-1)->flat->sum;
    }
    $out_string_sets .= "\n";
}
print_OUT("    ... done ...");

print_OUT("Writting results to [ $out.z_score.txt ]");
open(OUT,">$out.z_score.txt") or die $!;
print OUT $out_string_z_score;
close(OUT);

print_OUT("Writting results to [ $out.hard_score.txt ]");
open(OUT,">$out.hard_score.txt") or die $!;
print OUT $out_string_sets;
close(OUT);

if (defined $summary_annotation){
    my @yes_annot = ();
    my @not_annot = ();
    foreach  my $gene (keys %annotation){
        push @yes_annot, $annotation{$gene}->{scores} if (defined $annotation{$gene}->{sets});
        push @not_annot, $annotation{$gene}->{scores} if (not defined $annotation{$gene}->{sets});
    }
    my $y_annot = pdl @yes_annot;
    $y_annot = $y_annot->xchg(0,1);
    my $y_annot_sums = $y_annot->sumover;
    my $n_annot = pdl @not_annot;
    $n_annot = $n_annot->xchg(0,1);
    my $n_annot_sums = $n_annot->sumover;
    my $N_yes = scalar @yes_annot;
    my $N_not = scalar @not_annot;
    print_OUT("Calculating boostrap confidence intervals for annotation summaries");
    my %bootstrap_indexes = ();
    for (my $b = 0 ; $b < $bootstrap; $b++){
        $bootstrap_indexes{$b} = pdl get_index_sample($N_yes,$N_yes,1);
    }    
    my $YesBootstrap = zeroes scalar @phylostratum, $bootstrap;
    while (my ($b,$idx) = each %bootstrap_indexes){
        $YesBootstrap(,$b) .= $y_annot($idx,)->sumover->flat;
    }

    %bootstrap_indexes = ();
    for (my $b = 0 ; $b < $bootstrap; $b++){
        $bootstrap_indexes{$b} = pdl get_index_sample($N_not,$N_not,1);
    }
    my $NotBootstrap = zeroes scalar @phylostratum, $bootstrap;
    while (my ($b,$idx) = each %bootstrap_indexes){
        $NotBootstrap(,$b) .= $n_annot($idx,)->sumover->flat;
    }
    print_OUT("    ... done ...");
    print_OUT("Printing annotation summary per phylostratum to [ $out.annot_summary.txt ]");
    open (SUMMARY,">$out.annot_summary.txt")or die $!;
    print SUMMARY "RANK\tPhylostratum\twith_annotation\twith_0.025\twith_0.975\twithout_annotation\twithout_0.025\twithout_0.975\tannot_freq\tannot_freq_0.025\tannot_freq_0.975\n";
    for (my $i = 0; $i < scalar @phylostratum; $i++){
        my $total = $n_annot_sums($i) + $y_annot_sums($i);
        my $freq = 0;
        if ($total > 0){
            ($freq) = list $y_annot_sums($i)/$total;
        }
        my $low_yes = quantile($YesBootstrap($i,)->flat,'0.025',7);
        my $up_yes = quantile($YesBootstrap($i,)->flat,'0.975',7);
        
        my $low_not = quantile($NotBootstrap($i,)->flat,'0.025',7);
        my $up_not = quantile($NotBootstrap($i,)->flat,'0.975',7);
        
        my $low_freq  = 0;
        $low_freq = $low_yes/($low_yes + $low_not) if ($low_yes + $low_not > 0);
        my $up_freq = 0;
        $up_freq = $up_yes/($up_yes + $up_not) if ($up_yes + $up_not > 0);
        print SUMMARY $i,"\t";        
        print SUMMARY $phylostratum[$i],"\t";
        print SUMMARY $y_annot_sums($i)->list,"\t";
        print SUMMARY $low_yes,"\t",$up_yes,"\t";
        print SUMMARY $n_annot_sums($i)->list,"\t";
        print SUMMARY $low_not,"\t",$up_not,"\t";
        print SUMMARY $freq,"\t";
        print SUMMARY $low_freq,"\t",$up_freq;
        print SUMMARY "\n";
    }
}

print_OUT("Finished");

exit;


sub hygeometric_dist_z_score {
    my ($N,$n,$R,$r) = @_;
    my $R_over_N = $R/$N;
    my $z = ($r - $n*$R_over_N)/sqrt($n*($R_over_N)*(1-$R_over_N)*(1 - ($n - 1)/( $N - 1)));
    return($z);
}


__END__

=head1 NAME

 Calculation of enrichment scores of gene-sets of phylostratum.

=head1 DESCRIPTION

B<This program> will perform use hypegeometric distribution statistics to calculate the z-score enrichment statistic for enrichment of a gene-set among the genes of a phylostratum.

Comments, suggestions or complains you be addressed to IntiPedroso@gmail.com

=head1 SYNOPSIS

    ./z_score_phylostratum_enrichment.pl -annot annotation_file -score_matrix phylostratiphy_score_matrix -out output_file_name

    General options
    -h, --help		print help message
    -m, --man		print complete documentation

    Input
    -annot, -a      Gene annotation file.
    -score_matrix, -m   Hard-score matrix output from phylostratigraphy.pl
 
    Output
    -out, --o      Name of output files
 
    Additional Options
    -collapse, -s   Collapse counts of phylostratums. Need to provide id of phylostratum after which to collaps the counts.


=head1 OPTIONS

=over 8

=item B<-help>

    Print help message

=item B<-man>

    Print complete documentation

=item B<-annot, -a>

    Gene annotation file. Format is tab separed with two columns. First column has the gene-set id and second column the gene-id.

=item B<-score_matrix, -m>

    Hard-score matrix output from phylostratigraphy.pl

=item B<-o, --out>

    Name of output file
 
=intem B<-collapse, -s>
 
 Collapse counts of phylostratums. Need to provide id of phylostratum after which to collaps the counts.

=back

=head2 DESCRIPTION

 The reported z-score corresponds to (Latex expression): z-sco\mbox{re}=\frac{r\; -\; n\frac{R}{N}}{\sqrt{n\; \cdot \; \left( \frac{R}{n} \right)\left( 1-\frac{R}{n} \right)\left( 1-\frac{n-1}{N-1} \right)}}
 where R = is the total number of genes on the gene-set, N=total number of genes on genome, r = number of gene-set's genes on a phylostratum and n = total number of genes  on phylostratum of interest.
 
=head3 Examples

=over 8

=item B<1. Get taxonomy data>

TODO


 
=back
 
 
 
=cut
