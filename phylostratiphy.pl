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

our (   $help, $man, $tax_folder, $blast_out, $blast_format, $query_taxon, $out,
        $nucl_only, $prot_only, $seq_to_gi, $tax_info, $use_coverage, $virus_list);

GetOptions(
    'help' => \$help,
    'man' => \$man,
    'blast=s@' => \$blast_out,
    'tax_folder=s' => \$tax_folder,
    'blast_format=s' => \$blast_format,
    'nucl_only' => \$nucl_only,
    'prot_only' => \$prot_only,
    'seq_to_gi=s' => \$seq_to_gi,
    'tax_info=s' => \$tax_info,
    'use_coverage' => \$use_coverage,
    'query_taxon=s' => \$query_taxon,
    'virus_list=s' => \$virus_list,
    'out|o=s' => \$out,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);

$seq_to_gi ||= "seq_to_gi";
$tax_info ||= "tax_info";

print_OUT("Parsing taxonomy information");

# define the ids presents of the taxo-to-id mapping files
my $seq_id_files = 'both';
if ($nucl_only) {
    $seq_id_files = "nucl";
    print_OUT("Only working with protein sequences");
} elsif ($prot_only) {
    $seq_id_files = "prot";        
    print_OUT("Only working with nucleic acid sequences");
}
print_OUT("Mapping sequence ids to taxonomy ids");


# load tree of life information, both node' connections and names of nodes.
print_OUT("   '-> Reading phylogenetic tree and species information");
my $nodesfile = $tax_folder . "nodes.dmp";
my $namefile = $tax_folder . "names.dmp";
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => $nodesfile, -namesfile => $namefile);


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
        # get the p-value for the hit from the e-value.
        my $e_value = pdl $data[ $fields{'evalue'}];
        my $p_value = pdl E_CONSTANT**(-$e_value);
        $p_value = pdl 1 - $p_value;
        $p_value = pdl $e_value if ($p_value == 0);
        my $score = -1*($p_value->log) * $coverage;
        push @{ $S{ $data[ $fields{'query_id'}] }  }, { 'subject_id' => $subject_id[1], 'score' => $score, 'p_value' => $p_value, 'e_value' => $e_value };
        if ($e_value > 0 and $score >  $largest_hit_values{'score'}){
            %largest_hit_values = ('score' => $score->list, 'e_value' => $e_value->list, 'p_value' => $p_value->list );
        }
        $hits_gis{$subject_id[1]} = '';
    }
}

# check of entries equal to infinite. this happens when the e-value is equal to 0.
# on this cases we will replace the inf for the largest score available and its corresponding e-value. 

print_OUT("Finished processing blast output: [ $seq_counter ] sequences of which [ " . scalar (keys %S) . " ] have hits");
# for easy operation store the score of each gene on each specie on a matrix. Later internal nodes of the tree will be added as additional columns, that will make calculation of scores for new columns (internal nodes) faster.


my ($seq_to_tax_id,$target_taxons) = fetch_tax_ids_from_blastdb([keys %hits_gis] );


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
my $main_taxon = $db->get_taxon(-taxonid => $query_taxon);

print_OUT("Starting to calculate PhyloStratum Scores");
print_OUT("Identifiying last common ancestors between [ " . $main_taxon->scientific_name . " ] and [ " . scalar (keys %$target_taxons) . " ] target taxons");
print_OUT("   '-> Building tree with query and target species");
# get the tree using the taxonomy ids.
my $tree = get_tree_from_taxids($db,[$main_taxon->id,keys %$target_taxons]);
# remove redundant nodes, i.e., those with only one ancestor AND one descentdant.
$tree->contract_linear_paths;

# get the node for the taxon of interest
my $qry_node = $tree->find_node(-id => $main_taxon->id);
# get the root of the tree
my $tree_root = $tree->get_root_node;

print_OUT("   '-> Finding LCAs");
my %PATHS = ();
# tax counter is number of taxon minus 1 because the taxon counter starts from 0;
my $taxon_counter = (scalar (keys %$target_taxons) ) - 1;
# add nodes in lineafe of query node
@{ $PATHS{$qry_node->id} }= reverse $tree->get_lineage_nodes($qry_node);
foreach my $ancestor (@{ $PATHS{$qry_node->id} }){
    if (not exists $target_taxons->{$ancestor} ){
        $target_taxons->{$ancestor->id}->{'matrix_number'} = ++$taxon_counter;
        $target_taxons->{$ancestor->id}->{'scientific_name'} = $ancestor->scientific_name;
	}
    $target_taxons->{$ancestor->id}->{'scientific_name'} = $ancestor->scientific_name;
}   

foreach my $tree_leaf (keys %$target_taxons){
    # skip if the leave is the taxon of interest
    next if ($tree_leaf == $main_taxon->id);
    my $leaf_node = $tree->find_node(-id => $tree_leaf);
    if (not defined $leaf_node){
        print_OUT("Taxon not found for [ $tree_leaf ]");
        next;
    }
    # get LCA between leaf and query taxon
    my $lca = $tree->get_lca(($leaf_node,$qry_node));
    # get the path to the root of the target taxon
    @{ $PATHS{ $tree_leaf} }= reverse $tree->get_lineage_nodes($leaf_node);
    # check that the taxons have a taxon id to use for the matrix operations later
    foreach my $ancestor (@{ $PATHS{$tree_leaf} }){
        if (not exists $target_taxons->{$ancestor->id} ){
            $target_taxons->{$ancestor->id}->{'matrix_number'} = ++$taxon_counter;
        }
        $target_taxons->{$ancestor->id}->{'scientific_name'} = $ancestor->scientific_name;
    }
    $target_taxons->{$leaf_node->id}->{'scientific_name'} = $leaf_node->scientific_name;
}


print_OUT("Finishing to calculate scores");

my $M = zeroes scalar (keys %S), scalar (keys %$target_taxons);

my $gene_counter = 0;
foreach my $qry_gene (keys %S){
    foreach my $hit (@{ $S{$qry_gene} }){
        next if (not exists $seq_to_tax_id->{$hit->{'subject_id'}});
        my $subject_taxid = $seq_to_tax_id->{$hit->{'subject_id'}}->{'taxid'};
        my @taxon_idx = ();
        push  @taxon_idx, $target_taxons->{ $subject_taxid }->{'matrix_number'};
        foreach my $taxon (@{ $PATHS{$subject_taxid} }){ 
            push  @taxon_idx, $target_taxons->{ $taxon->id }->{'matrix_number'}; 
        }
        next if (scalar @taxon_idx == 0);        
        my $idx = pdl @taxon_idx;
        if ($hit->{'score'} eq "inf"){
            $M($gene_counter,$idx) += $largest_hit_values{'score'};
        } else {
            $M($gene_counter,$idx) += $hit->{'score'};
        }
    }
    $gene_counter++;
}

close(OUT);
# nomalize the scores.  

print_OUT("Printing query taxan Phylostratum Scores results to [ $out.qry_node_phylostratumscores.txt ]");
# get the phylostratum scores for the query node.
my $qry_node_ancestestors = pdl map { $target_taxons->{$_->id}->{'matrix_number'};} @{ $PATHS{$qry_node->id} };

# get taxon objects for the ancestors
my @qry_ancestors_array =  @{ $PATHS{$qry_node->id}};

for (my $idx = 1; $idx  < scalar @qry_ancestors_array; $idx++){
    my $present_taxa  = (@qry_ancestors_array)[$idx - 1];
    my $previous_taxa = (@qry_ancestors_array)[$idx];
    my $present_idx  = $target_taxons->{$present_taxa->id}->{'matrix_number'};
    my $previous_idx = $target_taxons->{$previous_taxa->id}->{'matrix_number'};
    $M(,$present_idx) -= $M(,$previous_idx);
}

for my $idx ( 0 .. scalar (keys %S) - 1){
    my $min_score = $M($idx,$qry_node_ancestestors)->min;
    if ($min_score < 0){
        $M($idx,$qry_node_ancestestors) -= $M($idx,$qry_node_ancestestors)->min;
    }
    $M($idx,$qry_node_ancestestors) /= $M($idx,$qry_node_ancestestors)->sum;
}

# replace nan to 0.
$M->inplace->setnantobad->inplace->setbadtoval(0);
my $qry_node_PhylostratumScores = $M(,$qry_node_ancestestors)->sumover;

open(OUT,">$out.qry_node_phylostratumscores.txt") or die $!;
print OUT join "\t", map { $target_taxons->{$_->id}->{'scientific_name'};} @{ $PATHS{$qry_node->id} };
print OUT "\n";
print OUT join "\t", $M(,$qry_node_ancestestors)->sumover->list;
print OUT "\n";
close(OUT);


print_OUT("Printing results to [ $out.txt ]");
open(OUT,">$out.txt") or die $!;
# print colnames of output file

print OUT join "\t", map { $target_taxons->{$_->id}->{'scientific_name'};} @{ $PATHS{$qry_node->id} };
print "\n";
$gene_counter = 0;
foreach my $qry_gene (keys %S){
    print OUT $qry_gene,"\t";
    print OUT join "\t", $M($gene_counter,$qry_node_ancestestors)->xchg(0,1)->list;
    print OUT "\n";
    $gene_counter++;
}
close(OUT);


print_OUT("Done");


exit;

