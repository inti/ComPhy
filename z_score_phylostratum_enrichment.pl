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
my $out_string = join "\t", ("GO",@phylostratum);
$out_string .= "\n";
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
    $out_string .= "$set\t" . join "\t", @z_scores;
    $out_string .= "\n";
}
print_OUT("    ... done ...");

print_OUT("Writting results to [ $out ]");
open(OUT,">$out") or die $!;
print OUT $out_string;
close(OUT);
print_OUT("Finished");

exit;


sub hygeometric_dist_z_score {
    my ($N,$n,$R,$r) = @_;
    my $R_over_N = $R/$N;
    my $z = ($r - $n*$R_over_N)/sqrt($n*($R_over_N)*(1-$R_over_N)*(1 - ($n - 1)/( $N - 1)));
    return($z);
}
