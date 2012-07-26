#!/usr/bin/perl -w
use strict;
use Bio::TreeIO; 
use Bio::DB::Taxonomy;
use Bio::Phylo::IO 'parse';
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use PDL;
use PDL::Matrix;
use PDL::NiceSlice;
use MLDBM qw(DB_File Storable);
use Fcntl;

use constant E_CONSTANT => log(10);
# local modules

use PhyloStratiphytUtils;

our (   $help, $man, $tax_folder, $blast_out, $blast_format, $user_provided_query_taxon_id, $out,
        $use_coverage, $virus_list, $hard_threshold,$gi_tax_id_info,$blastdbcmd,
        $seq_db);

GetOptions(
    'help' => \$help,
    'man' => \$man,
    'blast=s@' => \$blast_out,
    'tax_folder=s' => \$tax_folder,
    'blast_format=s' => \$blast_format,
    'use_coverage' => \$use_coverage,
    'query_taxon=s' => \$user_provided_query_taxon_id,
    'virus_list=s' => \$virus_list,
    'out|o=s' => \$out,
    'hard_threshold|hard=f' => \$hard_threshold,
    'gi_tax_id=s' => \$gi_tax_id_info,  # dysbindin.tax_info.csv
    'blastdbcmd=s' => \$blastdbcmd,
    'seq_db|db=s' => \$seq_db,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);


defined $blastdbcmd or $blastdbcmd = `which blastdbcmd`;
chomp($blastdbcmd);
defined $seq_db or $seq_db = "nr"; # assuming proteins and that path to dbs is on a enviromental variable

print_OUT("Parsing taxonomy information");

print_OUT("Mapping sequence ids to taxonomy ids");


# load tree of life information, both node' connections and names of nodes.
print_OUT("   '-> Reading phylogenetic tree and species information");
my $nodesfile = $tax_folder . "nodes.dmp";
my $namefile = $tax_folder . "names.dmp";
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => $nodesfile, -namesfile => $namefile);
my $tree_functions = Bio::Tree::Tree->new(); # load some tree functions

### define some variables to start storing the results

my %S = (); # hash will store to score for each species.

print_OUT("Starting to parse blast output");

my %target_taxons = ();
my $seq_counter = 0;
my %hits_gis = (); # store the gis of the hits
my %largest_hit_values = ('score' => -1, 'e_value' => 10000, 'p_value' => 1 );
# loop over blast results.
# for each hit we will store the log(p-value) of the blast hit for each taxon of the tartget sequence.
foreach my $file (@$blast_out){
    #print_OUT("   '-> [ $file ]");
    open (FILE,$file)or die $!;
    my %fields = ();
    while (my $line = <FILE>){
        # if line has the header of the table. split the header and keep the columns names.
        if ($line =~ m/^#/){
            if ($line =~ m/^# Fields:/){
                # clean up a bit the line
                $line =~ s/^# Fields: //;
                $line =~ s/\s$//;
                $line =~ s/\.//g;
                $line =~ s/\%/percent/g;
                # split the fields names on the comas.
                my @fs = split(/\,\s/,$line);
                # store the field names and positions on a hash
                for (my $i = 0; $i < scalar @fs; $i++) {
                    $fs[$i] =~ s/\s+/_/g;
                    $fields{$fs[$i]} = $i;
                }
                $seq_counter++;
            }
            next;
        }
        chomp($line);
        my @data = split(/\t+/,$line);
        # remove the frame of the db hit, we only need to know about the frame of the query seq
        if (exists $fields{'query/sbjct_frames'}){
            $data[ $fields{'query/sbjct_frames'} ] =~ s/\/\w+$//;
        }
        # extract the indetifiers of the target sequence
        my @subject_id = split(/\|/,$data[ $fields{'subject_id'}]);
        # define coverage as the fraction of the target sequence (subject) covered by the query sequence
        my $coverage = 1;
        if (defined $use_coverage) {
            $coverage = ($data[ $fields{'s_end'}] - $data[ $fields{'s_start'}])/$data[ $fields{'query_length'}]; 
        }
        # if using hard threshold then remove hits by e-value
        if (defined $hard_threshold){
            next if ($data[ $fields{'evalue'}] > $hard_threshold);
        }
        # get the p-value for the hit from the e-value.
        my $e_value = pdl $data[ $fields{'evalue'}];
        my $p_value = pdl E_CONSTANT**(-$e_value);
        $p_value = pdl 1 - $p_value;
        $p_value = pdl $e_value if ($p_value == 0);
        my $score = -1*($p_value->log) * $coverage;
        push @{ $S{ $data[ $fields{'query_id'}] }  }, { 'subject_id' => $subject_id[3], 'score' => $score, 'p_value' => $p_value, 'e_value' => $e_value };
        if ($e_value > 0 and $score >  $largest_hit_values{'score'}){
            %largest_hit_values = ('score' => $score->list, 'e_value' => $e_value->list, 'p_value' => $p_value->list );
        }
        $hits_gis{$subject_id[1]} = '';
    }
}

# check of entries equal to infinite. this happens when the e-value is equal to 0.
# on this cases we will replace the inf for the largest score available and its corresponding e-value. 

print_OUT("Finished processing blast output: [ $seq_counter ] sequences of which [ " . scalar (keys %S) . " ] have hits");

my ($seq_to_tax_id,$target_taxons) = fetch_tax_ids_from_blastdb([keys %hits_gis],$blastdbcmd,$seq_db,$out,$gi_tax_id_info, );

## remove hits to unwanted taxons
if (defined $virus_list){
    print_OUT("Parsing list of viral taxa to exclude");
    open(VL,$virus_list) or die $!;
    while(my $id = <VL>){
        chomp($id);
        delete($target_taxons->{$id});
    }
    close(VL);    
}


# get taxon information for the taxon of the query sequences.
my $query_taxon = $db->get_taxon(-taxonid => $user_provided_query_taxon_id);

print_OUT("Starting to calculate PhyloStratum Scores");
print_OUT("Identifiying last common ancestors between [ " . $query_taxon->scientific_name . " ] and [ " . scalar (keys %$target_taxons) . " ] target taxons");
print_OUT("   '-> Building tree with query and target species");

print_OUT("   '-> Finding LCAs");
# tax counter is number of taxon minus 1 because the taxon counter starts from 0;
my $taxon_counter = (scalar (keys %$target_taxons) ) - 1;
my %qry_ancestors = (); 
my $matrix_number = 0;
foreach my $taxon ( reverse $tree_functions->get_lineage_nodes($query_taxon) ){
    # here matrix numbers increase from tip to root of tree.
    $qry_ancestors{ $taxon->id } = { 'taxon' => $taxon, 'matrix_number' => $matrix_number++, 'scientific_name' => $taxon->scientific_name };
}

print_OUT("Finishing to calculate scores");

my $M = zeroes scalar (keys %S), scalar (keys %qry_ancestors);

# here store the lca between qry and target taxons.
my %LCA = ();
my %TAXON = ();
my %NODES = ();
my $gene_counter = 0;
foreach my $qry_gene (keys %S){
    my $oldest_lca_matrix_pos = -1;
    foreach my $hit (@{ $S{$qry_gene} }){
        # skip if info about the hits was not recovered from the sequence db
        next if (not exists $seq_to_tax_id->{ $hit->{'subject_id'} } );
        # get some info about the target taxon and the sequence
        my $subject_taxid = $seq_to_tax_id->{$hit->{'subject_id'}}->{'taxid'};
        my $subject_gi = $seq_to_tax_id->{$subject_taxid}->{'gi'};
        # next if hits is with qry taxon.
        next if ( $subject_taxid == $query_taxon->id); # CHECK THIS IS RIGHT
        # do not query for taxons from db more than once
        my $subject_taxon_from_db;
        if (not exists $NODES{ $subject_taxid } ){
            $subject_taxon_from_db = $db->get_taxon(-taxonid => $subject_taxid);
            $TAXON{ $subject_taxid } = $subject_taxon_from_db;
        } else {
            $subject_taxon_from_db = $TAXON{ $subject_taxid };
        }
        
        if (not exists $LCA{ $subject_taxid }){
            my $lca = $tree_functions->get_lca( ($subject_taxon_from_db,$query_taxon) );
            if (not defined $lca){
                my @lineage = $tree_functions->get_lineage_nodes($subject_taxon_from_db);
                if ($lineage[0] ne "cellular organisms") {
                    $LCA{ $subject_taxid } = 'none';
                    next;
                }
            }
            $LCA{ $subject_taxid } = $lca;
        }
        next if ($LCA{ $subject_taxid } eq 'none'); # skip if previously stablished that there is not an LCA.

        my $lca_matrix_pos = $qry_ancestors{ $LCA{ $subject_taxon_from_db->id }->id }->{'matrix_number'} ; 
        next if (not defined $lca_matrix_pos);
        if (defined $hard_threshold) {
            $oldest_lca_matrix_pos = $lca_matrix_pos if ($oldest_lca_matrix_pos < $lca_matrix_pos); # record the oldest LCA of the hits of this gene.
        } else {
            if ($hit->{'score'} eq "inf"){
                $M($gene_counter,$lca_matrix_pos) += $largest_hit_values{'score'};
            } else {
                $M($gene_counter,$lca_matrix_pos) += $hit->{'score'};
            }
        }
    }
    # score as 1 the oldest LCA of this gene.
    if (defined $hard_threshold){
        next if ($oldest_lca_matrix_pos == -1);
        $M($gene_counter,$oldest_lca_matrix_pos) .= 1;
    }
    $gene_counter++;
}

print_OUT("Printing query taxan Phylostratum Scores results to [ $out.qry_node_phylostratumscores.txt ]");

## normalise the scores by the sum of their logs, i.e., product of their probabilities.
unless (defined $hard_threshold) { # do not do it if using hard_threshold because the matrix has a single entry per gene equal to 1.
    $M /=$M->xchg(0,1)->sumover ; 
    $M->inplace->setnantobad->inplace->setbadtoval(0);
}

my $qry_node_PhylostratumScores = $M->sumover;

open(OUT,">$out.qry_node_phylostratumscores.txt") or die $!;
print OUT join "\t", map { $qry_ancestors{$_}->{"scientific_name"};} sort { $qry_ancestors{$a}->{'matrix_number'} <=> $qry_ancestors{$b}->{'matrix_number'} }  keys %qry_ancestors;
print OUT "\n";
print OUT join "\t", $qry_node_PhylostratumScores->list;
print OUT "\n";
close(OUT);


print_OUT("Printing results to [ $out.txt ]");
open(OUT,">$out.txt") or die $!;
# print colnames of output file

print OUT join "\t", map { $qry_ancestors{$_}->{"scientific_name"};} sort { $qry_ancestors{$a}->{'matrix_number'} <=> $qry_ancestors{$b}->{'matrix_number'} }  keys %qry_ancestors;
print OUT "\n";
$gene_counter = 0;
foreach my $qry_gene (keys %S){
    print OUT $qry_gene,"\t";
    print OUT join "\t", $M($gene_counter,)->xchg(0,1)->list;
    print OUT "\n";
    $gene_counter++;
}
close(OUT);


print_OUT("Done");


exit;

__END__

=head1 NAME
 
 Perl implementation of PAGE: parametric analysis of gene set enrichment. As bonus includes strategies to correct for gene clusters.
 
=head1 DESCRIPTION

B<This program> will perform a Phylostratigraphy analysis. It provides a implementation and extension of the original methodology of Domazet-Loso et. al. (2003).
 Please see readme for information on Perl module requirments and input and output files.
 
 For comments or complains please contact intipedroso@gmail.com
 
 References:
  Domazet-Loso, T., and Tautz, D. (2003). An evolutionary analysis of orphan genes in Drosophila. Genome Res. 13, 2213-2219.

=head1 SYNOPSIS

    ./phylostratiphy.pl [options]

    General options
    -h, --help		print help message
    -m, --man		print complete documentation
    Input Files
    -blast         Output from blast
    -tax_folder    Folder with taxonomy information
    -blast_format  Format of blast output
    -virus_list    List with ids of viral taxons
    -query_taxon   Taxnomy identifier of query specie.
 
    Output Files
    -o, --o        Name of output files

    Analysis options
    -hard, --hard_threshold    Switches to original methodology of Domazet-Loso et. al. (2003)
    -use_coverage  Means that scores are weightes by coverage of sequence alignment. Not compatible with -hard.
 
 
=head1 OPTIONS

=over 8

=item B<-help>

Print help message

=item B<-man>

Print complete documentation

=item B<-blast>
 Output from blast
 
=item B<-tax_folder>
 
 Folder with taxonomy information
 
=item B<-blast_format>  
 
 Format of blast output

=item B<-virus_list>    
 
 List with ids of viral taxons
 
=item B<-query_taxon>   
 
 Taxnomy identifier of query specie.
 
=item B<-o, --out>
 
 Name of output files

=item B<-hard, --hard_threshold>
 
 Switches to original methodology of Domazet-Loso et. al. (2003)

=item B<-use_coverage>
 
 Means that scores are weightes by coverage of sequence alignment. Not compatible with -hard.

=back
 
=head2 DESCRIPTION

TODO

=head3 Examples

=over 8

=item B<1. Basic Analysis>

To run the analyses with a small example do

>my_perl phylostratiphy.pl -blast dysbindin.blast_out.txt -tax_folder data/ -query_taxon 9606 -virus_list tmp.virus.txt -out test_phylostratiphy

The example consistes of human sequences and so -query_taxon 9606 corresponds to the human tax id.
MORE TO COME
 
 
=back



=cut

