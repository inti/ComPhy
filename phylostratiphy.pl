#!/usr/bin/perl -w
use strict;
use Bio::DB::Taxonomy;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;
use Bio::DB::EUtilities;
use Getopt::Long;
use Pod::Usage;
use PDL;
use PDL::Matrix;
use PDL::NiceSlice;
use Data::Dumper;

use constant E_CONSTANT => log(10);
my $EMAIL = 'intipedroso@gmail.com';
# local modules

use PhyloStratiphytUtils;

our (   $help, $man, $tax_folder, $blast_out, $blast_format, $user_provided_query_taxon_id, $out,
        $use_coverage, $hard_threshold,$gi_tax_id_info,$blastdbcmd,
        $seq_db);

GetOptions(
    'help' => \$help,
    'man' => \$man,
    'blast=s@' => \$blast_out,
    'tax_folder=s' => \$tax_folder,
    'blast_format=s' => \$blast_format,
    'use_coverage' => \$use_coverage,
    'query_taxon=s' => \$user_provided_query_taxon_id,
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
my $db = "";
my $tree_functions = Bio::Tree::Tree->new(); # load some tree functions

my $taxNCBI = Bio::LITE::Taxonomy::NCBI->new( db=>"NCBI", names=> $namefile, nodes=>$nodesfile, dict=>"$tax_folder/gi_taxid_prot.bin");
my @ql = $taxNCBI->get_taxonomy( $user_provided_query_taxon_id);

### define some variables to start storing the results

my %S = (); # hash will store to score for each species.

print_OUT("Starting to parse blast output");

my %gi_taxData = ();
my %lineages = (); # each entry has a 'lca_with_qry' and a 'lineage'
my %target_taxons = ();
my $seq_counter = 0;
my %hits_gis = (); # store the gis of the hits
my %largest_hit_values = ('score' => -1, 'e_value' => 10000, 'p_value' => 1 );
my %ids_not_found = ();
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
        # get the taxid and lineage for the sequence.
        if (not defined $gi_taxData{ $subject_id[1] }){
            my $sbjct_taxid = $taxNCBI->get_taxid( $subject_id[1] );
            if ($sbjct_taxid == 0){
                push @{ $ids_not_found{ $data[ $fields{'subject_id'}] }}, $data[ $fields{'query_id'}];
            } else {
                $gi_taxData{ $subject_id[1] } = { 'lineage' => [], 'tax_id' => $sbjct_taxid, 'lca_with_qry' => ""};
            }
        }
        if (not exists $lineages{ $gi_taxData{ $subject_id[1] }->{'tax_id'} }){
            my @sbjct_lineage = $taxNCBI->get_taxonomy( $gi_taxData{ $subject_id[1] }->{'tax_id'} );
            if (scalar @sbjct_lineage == 0){
                push @{ $ids_not_found{ $data[ $fields{'subject_id'}] }}, $data[ $fields{'query_id'}];
            } else {
                my $lca = get_lca_from_lineages(\@sbjct_lineage,\@ql); # need double checking on the ones that do not give match
                next if ($lca eq "diff_root"); # exclude those that have a different root to cell organisms.
                $gi_taxData{ $subject_id[1] }->{'lineage'} = \@sbjct_lineage;
                $gi_taxData{ $subject_id[1] }->{'lca_with_qry'} = $lca;
                $lineages{ $gi_taxData{ $subject_id[1] }->{'tax_id'} }->{'lca_with_qry'} = $lca;
                $lineages{ $gi_taxData{ $subject_id[1] }->{'tax_id'} }->{'lineage'} = \@sbjct_lineage;
            }
        } else {
            $gi_taxData{ $subject_id[1] }->{'lineage'} = $lineages{ $gi_taxData{ $subject_id[1] }->{'tax_id'} }->{'lineage'} ;
            $gi_taxData{ $subject_id[1] }->{'lca_with_qry'} = $lineages{ $gi_taxData{ $subject_id[1] }->{'tax_id'} }->{'lca_with_qry'} ;
        }
        # get the p-value for the hit from the e-value.
        my $e_value = pdl $data[ $fields{'evalue'}];
        my $p_value = pdl E_CONSTANT**(-$e_value);
        $p_value = pdl 1 - $p_value;
        $p_value = pdl $e_value if ($p_value == 0);
        my $score = -1*($p_value->log) * $coverage;
        push @{ $S{ $data[ $fields{'query_id'}] }  }, { 'subject_id' => $subject_id[1],
                                                        'score' => $score,
                                                        'p_value' => $p_value,
                                                        'e_value' => $e_value };
        if ($e_value > 0 and $score > $largest_hit_values{'score'}){
            %largest_hit_values = ( 'score' => $score->list,
                                    'e_value' => $e_value->list,
                                    'p_value' => $p_value->list );
        }
        $hits_gis{$subject_id[1]} = '';
    }
}

my @accs = ();
my %accs_to_gi = ();
foreach my $hidden (keys %ids_not_found){
    my @subject_id = split( /\|/,$hidden );
    push @accs, ($subject_id[3] =~ m/(.*)\.\d+$/);
    $accs_to_gi{ $accs[-1] } = $subject_id[1];
}
if (scalar @accs > 0){
    print_OUT("There were [ " . scalar @accs . " ] ids unmatched. Using webservices to get taxonomy information on them.");
    my $factory = Bio::DB::EUtilities->new( -eutil => 'esearch',
                                            -email => $EMAIL,
                                            -db    => 'protein',
                                            -retmax     => 10*(scalar @accs),
                                            -term  => join(',',@accs),
                                            -usehistory => 'y');
    
    my $hist = $factory->next_History || die print_OUT("No history data returned from esearch");
    my @uids = $factory->get_ids;
    $factory->set_parameters(   -eutil => 'esummary',
                                -db => 'protein',
                                -id => \@uids,
                                -email => $EMAIL,
                                -history => $hist);
#    $factory->reset_parameters(-eutil => 'esummary', -db => 'protein', -id => \@uids,-email => $EMAIL);
    while (my $ds = $factory->next_DocSum) {
        my ($taxid) = $ds->get_contents_by_name("TaxId");
        my ($caption) = $ds->get_contents_by_name("Caption");
        if (not defined $lineages{ $taxid }){
            my @sbjct_lineage = $taxNCBI->get_taxonomy( $taxid );
            my $lca = get_lca_from_lineages(\@sbjct_lineage,\@ql); # need double checking on the ones that do not give match
            next if ($lca eq "diff_root"); # exclude those that have a different root to cell organisms.
            $gi_taxData{ $accs_to_gi{ $caption } }->{'lineage'} = \@sbjct_lineage;
            $gi_taxData{ $accs_to_gi{ $caption } }->{'lca_with_qry'} = $lca;
            $lineages{ $taxid }->{'lineage'} =  \@sbjct_lineage;
            $lineages{ $taxid }->{'lca_with_qry'} = $lca;
        }
    }
}


print_OUT("Finished processing blast output: [ $seq_counter ] sequences of which [ " . scalar (keys %S) . " ] have hits");

my %qry_ancestors_ranks = ();
my $c = 0;
map { $qry_ancestors_ranks{$_} = $c++;   } @ql;

my %PhyloStratum_scores = ();
foreach my $qry_seq (keys %S) {
    my @ranks = sort {$a <=> $b} map { $qry_ancestors_ranks{ $gi_taxData{$_->{'subject_id'}}->{'lca_with_qry'} }; } @{$S{$qry_seq}};
    my @phyloScores = list zeroes scalar @ql;
    $phyloScores[$ranks[0]] = 1;
    $PhyloStratum_scores{$qry_seq} = \@phyloScores;
}

# Print out the PhyloStratumScores as hardcoded
open(OUT,">$out.qry_node_phylostratumscores.txt") or die $!;
print OUT "ID \t",join "\t", @ql;
print OUT "\n";
print OUT join "\t", list sumover mpdl values %PhyloStratum_scores;
print OUT "\n";
close(OUT);
exit;

# TODO:
#1 print out output for hard coded analysis.
#2 solve issues with gi not found due to changes on the sequence db
#3 implement soft-coded



my ($seq_to_tax_id,$target_taxons) = fetch_tax_ids_from_blastdb([keys %hits_gis],$blastdbcmd,$seq_db,$out,$gi_tax_id_info);

# get taxon information for the taxon of the query sequences.
my $query_taxon = $db->get_taxon(-taxonid => $user_provided_query_taxon_id);

print_OUT("Starting to calculate PhyloStratum Scores");
print_OUT("   '-> Will identifiying last common ancestors between [ " . $query_taxon->scientific_name . " ] and [ " . scalar (keys %$target_taxons) . " ] target taxons");
# tax counter is number of taxon minus 1 because the taxon counter starts from 0;
my $taxon_counter = (scalar (keys %$target_taxons) ) - 1;
my %qry_ancestors = (); 
my $matrix_number = 0;
# add query species first
foreach my $taxon ( ($query_taxon,reverse $tree_functions->get_lineage_nodes($query_taxon)) ){
    # here matrix numbers increase from tip to root of tree.
    $qry_ancestors{ $taxon->id } = { 'taxon' => $taxon, 'matrix_number' => $matrix_number++, 'scientific_name' => $taxon->scientific_name };
}

print_OUT("   '-> Finishing to calculate scores");

my $M = zeroes scalar (keys %S), scalar (keys %qry_ancestors);

# here store the lca between qry and target taxons.
my %LCA = ();
my %TAXON = ();
my %NODES = ();
my $gene_counter = 0;
my $n_genes = scalar (keys %S);
foreach my $qry_gene (keys %S){
    # print counter
    print scalar localtime,"\t",progress_bar($gene_counter,$n_genes);
    
    my $oldest_lca_matrix_pos = -1;
    foreach my $hit (@{ $S{$qry_gene} }){
        # skip if info about the hits was not recovered from the sequence db
        next if (not exists $seq_to_tax_id->{ $hit->{'subject_id'} } );
        # get some info about the target taxon and the sequence
        my $subject_taxid = $seq_to_tax_id->{$hit->{'subject_id'}}->{'taxid'};
        my $subject_gi = $seq_to_tax_id->{$subject_taxid}->{'gi'};
        my $subject_name = $seq_to_tax_id->{$subject_taxid}->{'scientific_name'};
        # do not query for taxons from db more than once
        my $subject_taxon_from_db;
        if (not exists $NODES{ $subject_taxid } ){
            $subject_taxon_from_db = $db->get_taxon(-taxonid => $subject_taxid);
            $TAXON{ $subject_taxid } = $subject_taxon_from_db;
        } else {
            $subject_taxon_from_db = $TAXON{ $subject_taxid };
        }
        if (not defined $subject_taxon_from_db) {
	    $subject_taxon_from_db = $db->get_taxon(-name => $subject_taxid);
            $TAXON{ $subject_taxid } = $subject_taxon_from_db;
	} 
	if (not defined $subject_taxon_from_db) {
		print_OUT("\nDid not find a taxon for target specie with id [ $subject_taxid ] and name [ $subject_name ]");
		next;
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
print "\n"; # to finish the counter above
print_OUT("   '-> Printing query taxan Phylostratum Scores results to [ $out.qry_node_phylostratumscores.txt ]");

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


print_OUT("   '-> Printing results to [ $out.txt ]");
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

if (not defined $gi_tax_id_info){
    print_OUT("   '-> Sequence to taxonomy information has bee stored at [ $out.gi_tax_id.csv ]");
}


print_OUT("Done");


exit;

__END__

=head1 NAME
 
 Perl implementation of the PhyloStratigraphy comparative genomics methods.
 
=head1 DESCRIPTION

B<This program> will perform a PhyloStratigraphy analysis. It provides a implementation and extensions of the original methodology of Domazet-Loso et. al. (2003).
 Please see README for information on Perl module requirments and some additional information.
 
 Comments, suggestions or complains you be addressed to IntiPedroso@gmail.com
 
 References:
  Domazet-Loso, T., and Tautz, D. (2003). An evolutionary analysis of orphan genes in Drosophila. Genome Res. 13, 2213-2219.

=head1 SYNOPSIS

    ./phylostratiphy.pl [options]

    General options
    -h, --help		print help message
    -m, --man		print complete documentation
 
    Input Files
    -blast         Output from blast (Required)
    -tax_folder    Folder with taxonomy information (Required)
    -blast_format  Format of blast output (Not yet functional)
    -query_taxon   Taxnomy identifier of query specie. (Required)
    -seq_db        Path or name to sequence db in blast format. (Optional)
    -blastdbcmd    Path to blastdbcmd executable (Optional)
    -gi_tax_id     Tabe with taxonomy information for the target sequences of the blast output (Optional)
 
    Output Files
    -o, --o        Name of output files (Required)

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
 
 Taxnomy identifier of query specie. User must use the id provided by the NCBI Taxonomy database, e.g, human = 9606 and Acromyrmex echinatior = 103372.
 
=item B<-seq_db> 
 
 Path or name to sequence db in blast format. If nor provided it will be set to 'nr' and it will be assumed that you have set the  ~/.ncbirc as indicated on the BLAST documentation with the path to the folder where the sequence databased are stored, i.e., the $BLASTDB enviromental variable has been set.

=item B<-blastdbcmd>

 Path to blastdbcmd executable. If not provided it will be assumed that it can found with the command `which blastdbcmd`

=item B<-gi_tax_id>
 
 Tabe with taxonomy information for the target sequences of the blast output. This table has the format obtained by extrating the taxonomy information the sequence db with the following command (as it is done internally by the software)
 >blastdbcmd -outfmt "%a,%g,%T,%L,%S" -entry_batch file_with_sequence_ids.txt -db seq_db -out output_file.csv
 where seq_db is whatever the -seq_db as been set to and output_file.csv is the file that you need to provide to the -gi_tax_id option. If this option is not provided this command will be run internally and you will get the file as "out".gi_tax_id.csv, where "out" is the whatever you provide to the option -o, --out

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

