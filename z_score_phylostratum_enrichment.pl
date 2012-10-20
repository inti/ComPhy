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

our (   $help, $man, $out, $score_matrix, $annot );

GetOptions(
        'help' => \$help,
        'man' => \$man,
        'out|o=s' => \$out,
        'score_matrix|m=s' => \$score_matrix,
        'annot|a=s' => \$annot,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);
pod2usage(-exitstatus => 2, -verbose => 2) if (not defined $annot and not defined $score_matrix);



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
print_OUT("Calculating enrichment z-scores");
my $phylostratum_totals = $phylostratum_scores->sumover->flat;
my $N = $phylostratum_totals->sum;
my $out_string_z_score = join "\t", ("GO",@phylostratum);
$out_string_z_score .= "\n";
my $out_string_sets = join "\t", ("GO",@phylostratum);
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
            push @z_scores, 0;
            next;
        }
        push @z_scores, hygeometric_dist_z_score($N,$n,$R,$r);
    }
    $out_string_z_score .= "$set\t" . join "\t", @z_scores;
    $out_string_z_score .= "\n";
    $out_string_sets .= "$set\t" . join "\t", $profile->list;
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
