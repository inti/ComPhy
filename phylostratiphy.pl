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


# local modules

use PhyloStratiphytUtils;

our (   $help, $man, $tax_folder, $blast_out, $blast_format, $query_taxon,
        $nucl_only, $prot_only, $seq_to_gi, $tax_info, $use_coverage);

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
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);

$seq_to_gi ||= "seq_to_gi";
$tax_info ||= "tax_info";

print_OUT("Parsing taxonomy information");

my $seq_id_files = 'both';
if ($nucl_only) {
    $seq_id_files = "nucl";
    print_OUT("Only working with protein sequences");
} elsif ($prot_only) {
    $seq_id_files = "prot";        
    print_OUT("Only working with nucleic acid sequences");
}
print_OUT("Mapping sequence ids to taxonomy ids");

my (%seq_to_tax_db,%tax_info_db,%tax_tree_db);
if (defined $tax_folder){
    if (not -e "$tax_folder/$seq_to_gi.db") {
        print_OUT("Moving into [ $tax_folder ] to create db");
        chdir($tax_folder);
        unlink("$seq_to_gi.db");
        %seq_to_tax_db  = build_database($seq_to_gi,$tax_folder,$seq_id_files);   
        chdir("../");
    }
    print_OUT("Openning databases");
    print_OUT("   '-> [ $tax_folder/$seq_to_gi.db ]");
    tie %seq_to_tax_db, "DB_File", "$tax_folder/$seq_to_gi.db" or die "Cannot open db file [ $tax_folder/$seq_to_gi.db ]: $!\n";
    
} else {
    print_OUT("Taxonomy information will downloaded and stored on folder [ data ]");
    $tax_folder = "data";
    unless (-d $tax_folder){ mkdir($tax_folder);}
    chdir($tax_folder);
    %seq_to_tax_db  = build_database($seq_to_gi,$tax_folder,$seq_id_files);   
    chdir("../");
}

print_OUT("   '-> Reading phylogenetic tree and species information");
my $nodesfile = $tax_folder . "nodes.dmp";
my $namefile = $tax_folder . "names.dmp";
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile', -nodesfile => $nodesfile, -namesfile => $namefile);


### define some variables to start storing the results

my %S = (); # hash will store to score for each species.


print_OUT("Starting to parse blast output");
$seq_to_tax_db{"322792145"} = 13686; # this is to be removed 

my %target_taxons = ();
my $seq_counter = 0;
foreach my $file (@$blast_out){
    #print_OUT("   '-> [ $file ]");
    open (FILE,$file)or die $!;
    my %fields = ();
    while (my $line = <FILE>){
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
            }
            $seq_counter++;
            next;
        }
        chomp($line);
        my @data = split(/\t+/,$line);
        # remove the frame of the db hit, we only need to know about the frame of the query seq
        if (exists $fields{'query/sbjct_frames'}){
            $data[ $fields{'query/sbjct_frames'} ] =~ s/\/\w+$//;
        }
        my @subject_id = split(/\|/,$data[ $fields{'subject_id'}]);
        next if (not exists $seq_to_tax_db{$subject_id[1]} ); # to be removed later
        next if (not defined $db->get_taxon(-taxonid => $seq_to_tax_db{$subject_id[1]})); # to avoid crash later when recovering the taxon objects.
        my $coverage = 1;
        if (defined $use_coverage) {
            $coverage = ($data[ $fields{'s_end'}] - $data[ $fields{'s_start'}])/$data[ $fields{'query_length'}]; 
        }
        if ($data[ $fields{'evalue'}] == 0){
            $S{ $data[ $fields{'query_id'}] }{ $seq_to_tax_db{$subject_id[1]} } += 500;
        } else {
            next if ((-log ($data[ $fields{'evalue'}])) < 0); # to be removed when using p-value instead of e-value.
            $S{ $data[ $fields{'query_id'}] }{ $seq_to_tax_db{$subject_id[1]} } += (-log ($data[ $fields{'evalue'}])) * $coverage;
        }
        $target_taxons{$seq_to_tax_db{$subject_id[1]}} = "";
    }
}
print_OUT("Finished processing blast output: [ $seq_counter ] sequences of which [ " . scalar (keys %S) . " ] have hits");

# for easy operation store the score of each gene on each specie on a matrix. Later internal nodes of the tree will be added as additional columns, that will make calculation of scores for new columns (internal nodes) faster.
# print score matrix;

my $S_g_f = [];
my %S_g_f_gene_idx = ();
my %S_g_f_taxon_idx = ();
my $gene_counter = 0;
my $taxon_counter = 0;

foreach my $spc (sort {$a cmp $b} keys %target_taxons){
    $S_g_f_taxon_idx{$spc} = $taxon_counter++;
}

while (my ($gene, $target_species) = each %S){
    $S_g_f_gene_idx{$gene} = $gene_counter++;
    foreach my $spc (sort {$a cmp $b} keys %target_taxons)  {
        if (not exists $target_species->{$spc}) { 
            $S_g_f->[ $S_g_f_gene_idx{$gene}  ] [ $S_g_f_taxon_idx{$spc} ] = 0;
        } else {
            $S_g_f->[ $S_g_f_gene_idx{$gene}  ] [ $S_g_f_taxon_idx{$spc} ] = $target_species->{$spc};
        }
    }
}

# create matrix in PDL format
$S_g_f = mpdl $S_g_f;
$S_g_f /=  $S_g_f->xchg(0,1)->sumover;

# get taxon information for the query taxon.
my $main_taxon = $db->get_taxon(-taxonid => $query_taxon);


print_OUT("Starting to calculate PhyloStratum scores");
print_OUT("Identifiying last common ancestors between [ " . $main_taxon->scientific_name . " ] and [ " . scalar (keys %target_taxons) . " ] target species");

# get target species and add the query specie
my @species_names = map { $db->get_taxon(-taxonid => $_)->scientific_name;  } keys %target_taxons;
# obtain a tree containing only the species of interest
my $tree = $db->get_tree((@species_names,$main_taxon->scientific_name));
# remove redundant nodes, i.e., those with only one ancestor AND one descentdant.
$tree->contract_linear_paths;

# get the node for the taxon of interest
my $qry_node = $tree->find_node($main_taxon->id);

my %LCA = (); # this hash stores the the last-common ancestor between the query species and the target specie. This node represent oldest node to which the score of the species needs to be addedd.

# get the root of the tree
my $tree_root = $tree->get_root_node;
# loop over each of the tree leaves
foreach my $node_id (@{ return_all_Leaf_Descendents($tree_root) }){
    # skip if the leave is the taxon of interest
    next if ($node_id->id == $main_taxon->id);
    # get the last common ancestor of the target and taxon of interest
    my $lca = $tree->get_lca(($node_id,$qry_node));
    $LCA{$lca->id} = $lca;
    # get the path to the root of the target taxon
    my @path = reverse $tree->get_lineage_nodes($node_id);
    # loop from the target taxon to the LCA and add the target taxon score to each ancestor
    foreach my $ancestor ( @path){
        if (exists $S_g_f_taxon_idx{$node_id->ancestor->id}){
            my $matrix_idx = $S_g_f_taxon_idx{$node_id->ancestor->id};
            $S_g_f(,$matrix_idx) += $S_g_f(,$matrix_idx);
        } else {
            # record its position on the matrix
            $S_g_f_taxon_idx{$node_id->ancestor->id} = $taxon_counter++;
            # create new column
            $S_g_f = $S_g_f->glue(1, $S_g_f(,$S_g_f_taxon_idx{$node_id->id}) );
        }
        last if ($ancestor->id == $lca->id);
    }
}



my $last_recorded_lca = '';

my @path_to_taxon = reverse $tree->get_lineage_nodes($main_taxon);
for (my $i = 0; $i < scalar @path_to_taxon; $i++){
    $last_recorded_lca = $i if (exists $S_g_f_taxon_idx{ $path_to_taxon[$i]->id() });
    next if ($last_recorded_lca eq '');
    if (not exists $S_g_f_taxon_idx{ $path_to_taxon[$i]->id() }){
        $S_g_f_taxon_idx{ $path_to_taxon[$i]->id() } = $S_g_f_taxon_idx{ $path_to_taxon[ $last_recorded_lca ]->id };
    }
}

my $I_g_f;
for (my $i = 0; $i < scalar @path_to_taxon; $i++){
    next if (not exists $S_g_f_taxon_idx{ $path_to_taxon[$i]->id() });
    last if (not defined $path_to_taxon[$i]->ancestor);
    my $idx_current = $S_g_f_taxon_idx{ $path_to_taxon[$i]->id() };
    my $idx_previous = $S_g_f_taxon_idx{ $path_to_taxon[$i]->ancestor->id() };
    my $tmp_col = $S_g_f(,$idx_current) - $S_g_f(,$idx_previous);
    push @{ $I_g_f }, [$tmp_col->list];
}
$I_g_f = pdl $I_g_f;
my $S_f = $I_g_f->sumover; # here store the cummulative score for each stratum over all genes.

print "I_g_f\n";
print join " ",$I_g_f->dims;
print "S_f\n";
print $S_f;
getc;

print_OUT("Done");


exit;


