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
use NCBI_PowerScripting;

our (   $help, $man, $tax_folder, $blast_out, $blast_format, $user_provided_query_taxon_id, $out,
        $use_coverage, $hard_threshold, $soft_threshold, $gi_tax_id_info,$blastdbcmd,
        $seq_db, $not_use_ncbi_entrez, $guess_qry_specie, $ncbi_entrez_batch_size, $max_eutils_queries,$EMAIL );

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
    'soft_threshold|soft' => \$soft_threshold,
    'gi_tax_id=s' => \$gi_tax_id_info,  # dysbindin.tax_info.csv
    'blastdbcmd=s' => \$blastdbcmd,
    'seq_db|db=s' => \$seq_db,
    'no_ncbi_entrez' => \$not_use_ncbi_entrez,
    'ncbi_entrez_batch_size=i' => \$ncbi_entrez_batch_size,
    'guess_qry_specie' => \$guess_qry_specie,
    'max_eutils_queries=i' => \$max_eutils_queries,
    'email' => \$EMAIL,
) or pod2usage(0);

pod2usage(0) if (defined $help);
pod2usage(-exitstatus => 2, -verbose => 2) if (defined $man);
pod2usage(-exitstatus => 2, -verbose => 2) if (not defined $EMAIL);


#### DEFINE SOME DEFAULT VALUES #############
defined $blastdbcmd or $blastdbcmd = `which blastdbcmd`;
chomp($blastdbcmd);
defined $seq_db or $seq_db = "nr"; # assuming proteins and that path to dbs is on a enviromental variable
defined $hard_threshold or $hard_threshold = 1e-3;
defined $soft_threshold and $hard_threshold = undef;
defined $blast_format or $blast_format = 'table';
defined $ncbi_entrez_batch_size or $ncbi_entrez_batch_size = 500;
$ncbi_entrez_batch_size = 500 if ($ncbi_entrez_batch_size > 500);
defined $max_eutils_queries or $max_eutils_queries = 0;


print_OUT("Starting to parse blast output");

### define some variables to start storing the results
my %S = (); # hash will store to score for each species.
# loop over blast results.
# for each hit we will store the log(p-value) of the blast hit for each taxon of the tartget sequence.
foreach my $file (@$blast_out){
    #print_OUT("   '-> [ $file ]");
    if ($blast_format eq 'table'){
        my $parsed_blast_out = parse_blast_table($file);
        @S{keys %{$parsed_blast_out}} = values %{$parsed_blast_out};
    } elsif ($blast_format eq 'paralign'){
        my $parsed_blast_out = parse_paralign_table($file);
        @S{keys %{$parsed_blast_out}} = values %{$parsed_blast_out};
    }
}
print_OUT("Finished processing blast output with results for [ " . scalar (keys %S) . " ] sequences.");

#### LOAD TAXONOMY DB
print_OUT("Reading taxonomy information");

# load tree of life information, both node' connections and names of nodes.
print_OUT("   '-> Reading phylogenetic tree and species information");
my $nodesfile = $tax_folder . "nodes.dmp";
my $namefile = $tax_folder . "names.dmp";
my $taxNCBI = Bio::LITE::Taxonomy::NCBI->new( db=>"NCBI", names=> $namefile, nodes=>$nodesfile, dict=>"$tax_folder/gi_taxid_prot.bin");


######## FILTER BLAST HITS AND GET TAXONMY INFORMATION FOR TARGET SEQUENCES AND SPECIES.
my %gi_taxData = (); # store taxonomy information for each target gi
my %lineages = (); # each entry has a 'lca_with_qry' and a 'lineage'
my $seq_counter = 0;
my %ids_not_found = (); # store ids of target sequences without taxonomy information on local DBs

#### if requested guess the query specie identity
my %specie_guess = ();
if (defined $guess_qry_specie){
    print_OUT("Query specie id not provided. Guessing query specie based on blast results.");
    while (my ($qry_seq,$blast_subjects) = each %S){
        next if (scalar @{$blast_subjects} == 0);
        # loop over blast results for this query sequence
        my $target_counter = -1; # set to -1 so that first item sets it to 0.
        foreach my $target_seqs (@{$blast_subjects}){
            if (ref($target_seqs) eq 'ARRAY'){
                next if (scalar @{$target_seqs} == 0);
            }
            next if ($target_seqs->{'evalue'} > $hard_threshold);
            my $sbjct_taxid = $taxNCBI->get_taxid( $target_seqs->{'subject_id'} );
            if ($sbjct_taxid == 0){
                push @{ $ids_not_found{ $target_seqs->{'subject'} }}, $qry_seq;
            } else {
                $specie_guess{$sbjct_taxid}++ if ( $target_seqs->{'percent_identity'} == 100);
                $gi_taxData{ $target_seqs->{'subject_id'} } = { 'lineage' => [],
                                                                'tax_id' => $sbjct_taxid,
                                                                'lca_with_qry' => ''};
            }
        }
    }
    ($user_provided_query_taxon_id) = sort {$specie_guess{$b} <=> $specie_guess{$a} } keys %specie_guess;
    if (not defined $user_provided_query_taxon_id){
        print_OUT("I could not guess the query specie.");
        print_OUT("Bye :(.");
        exit;
    }
    print_OUT("   '-> I am gessing that query specie NCBI Taxonomy id is [ $user_provided_query_taxon_id ].");
}

#### Obtain lineage of query specie.
my @ql = $taxNCBI->get_taxonomy( $user_provided_query_taxon_id);

# First pass over the data to get the taxonomy id of each target sequence
while (my ($qry_seq,$blast_subjects) = each %S){
    # loop over blast results for this query sequence
    my $target_counter = -1; # set to -1 so that first item sets it to 0.
    foreach my $target_seqs (@{$blast_subjects}){
        $target_counter++;
        if (ref($target_seqs) eq 'ARRAY'){
            next if (scalar @{$target_seqs} == 0); # skip if this query did not have blast hits.
        }
        ## Filter blast results.
        # if using hard threshold then remove hits by e-value
        if (defined $hard_threshold){
            if ($target_seqs->{'evalue'} > $hard_threshold){
                # flag this target sequence as not passing the filter
                $S{ $qry_seq }->[ $target_counter ]->{'pass_hard_thresold'} = 0;
                # if skip the taxonomy of this seq unless we are interested in doing soft-threshold method.
                unless (defined $soft_threshold ){
                    next;
                }
            } else {
                $S{ $qry_seq }->[ $target_counter ]->{'pass_hard_thresold'} = 1;
            }
        }
        if (not defined $gi_taxData{ $target_seqs->{'subject_id'} }){
            my $sbjct_taxid = $taxNCBI->get_taxid( $target_seqs->{'subject_id'} );
            if ($sbjct_taxid == 0){
                push @{ $ids_not_found{ $target_seqs->{'subject'} }}, $qry_seq;
            } else {
                $gi_taxData{ $target_seqs->{'subject_id'} } = { 'lineage' => [],
                                                                'tax_id' => $sbjct_taxid,
                                                                'lca_with_qry' => ''};
                if (not exists $lineages{ $sbjct_taxid }){
                    my @sbjct_lineage = $taxNCBI->get_taxonomy( $sbjct_taxid );
                    if (scalar @sbjct_lineage == 0){
                        push @{ $ids_not_found{ $target_seqs->{'subject'} }}, $qry_seq;
                    } else {
                        my $lca = get_lca_from_lineages(\@sbjct_lineage,\@ql);
                        if ($lca eq "diff_root") {# exclude those that have a different root to cell organisms.
                            $gi_taxData{ $target_seqs->{'subject_id'} }->{'lineage'} = \@sbjct_lineage;
                            $gi_taxData{ $target_seqs->{'subject_id'} }->{'lca_with_qry'} = $lca;
                            next;
                        }
                        $gi_taxData{ $target_seqs->{'subject_id'} }->{'lineage'} = \@sbjct_lineage;
                        $gi_taxData{ $target_seqs->{'subject_id'} }->{'lca_with_qry'} = $lca;
                        $lineages{ $sbjct_taxid }->{'lca_with_qry'} = $lca;
                        $lineages{ $sbjct_taxid }->{'lineage'} = \@sbjct_lineage;
                    }
                } else {
                    $gi_taxData{ $target_seqs->{'subject_id'} }->{'lineage'} = $lineages{ $sbjct_taxid }->{'lineage'} ;
                    $gi_taxData{ $target_seqs->{'subject_id'} }->{'lca_with_qry'} = $lineages{ $sbjct_taxid }->{'lca_with_qry'} ;
                }
            }
        } else {
            my $sbjct_taxid = $gi_taxData{ $target_seqs->{'subject_id'} }->{'tax_id'};
            if (not exists $lineages{ $sbjct_taxid }){
                my @sbjct_lineage = $taxNCBI->get_taxonomy( $sbjct_taxid );
                if (scalar @sbjct_lineage == 0){
                    push @{ $ids_not_found{ $target_seqs->{'subject'} }}, $qry_seq;
                } else {
                    my $lca = get_lca_from_lineages(\@sbjct_lineage,\@ql);
                    if ($lca eq "diff_root") {# exclude those that have a different root to cell organisms.
                        $gi_taxData{ $target_seqs->{'subject_id'} }->{'lineage'} = \@sbjct_lineage;
                        $gi_taxData{ $target_seqs->{'subject_id'} }->{'lca_with_qry'} = $lca;
                        next;
                    }
                    $gi_taxData{ $target_seqs->{'subject_id'} }->{'lineage'} = \@sbjct_lineage;
                    $gi_taxData{ $target_seqs->{'subject_id'} }->{'lca_with_qry'} = $lca;
                    $lineages{ $sbjct_taxid }->{'lca_with_qry'} = $lca;
                    $lineages{ $sbjct_taxid }->{'lineage'} = \@sbjct_lineage;
                }
            } else {
                $gi_taxData{ $target_seqs->{'subject_id'} }->{'lineage'} = $lineages{ $sbjct_taxid }->{'lineage'} ;
                $gi_taxData{ $target_seqs->{'subject_id'} }->{'lca_with_qry'} = $lineages{ $sbjct_taxid }->{'lca_with_qry'} ;
            }        
        }
    }

}


######## FIND TAXONOMY INFORMATION FOR IDS WITH INCOSISTENT INFORMATION ON LOCAL FILES ################
# compile all ids for which taxonomy information was not found with local files.
my $try_again_eutils = 1;
my $last_unmatched = 0;
TRYAGAIN:
my @accs = ();
my %accs_to_gi = ();
foreach my $hidden (keys %ids_not_found){
    my @subject_id = split( /\|/,$hidden );
     my $accn = $subject_id[3];
($accn) =~ s/\.\d+$//;
    push @accs, $subject_id[1];
    #push @accs, $accn;
    $accs_to_gi{ $accn } = { 'accn' => $accn, 'gi' => $subject_id[1], 'full_id' => $hidden};
    if (defined $subject_id[4]){
	$accs_to_gi{ $subject_id[4] } = { 'accn' => $subject_id[4], 'gi' => $subject_id[1], 'full_id' => $hidden};
	push @accs, $subject_id[4];
    }
    #$accs_to_gi{ $subject_id[1] } = { 'accn' => $accn, 'gi' => $subject_id[1], 'full_id' => $hidden};
}
print_OUT("There were [ " . (scalar @accs) . " ] subject ids from blast search unmatched with local taxonomy DBs");
if (defined $not_use_ncbi_entrez){
    print_OUT("   '-> Printing this ids to [ $out.ids_without_taxonomy_information.txt ].");
    open(OUT,">$out.ids_without_taxonomy_information.txt") or die $!;
    print OUT join "\n", @accs;
    close(OUT);
} else {
    if (scalar @accs > 0){
        print_OUT("   '-> Using NCBI webservices to get taxonomy information on them. Trial [ $try_again_eutils ]");
        print_OUT("   '-> It will use approximately [ " . round_up((scalar @accs)/$ncbi_entrez_batch_size) . " ] queries.");
        while (scalar @accs > 0){
            my $batch_size = 1000;
            $batch_size = scalar @accs if (scalar @accs < 5000);
            my @tmp_accs_list = splice(@accs,0, $batch_size);
            my $id_list = join ",",@tmp_accs_list;
            my %params = (
                            'db' => 'protein',
                            'term' => $id_list,
                            'id' => $id_list,
                            'email' => $EMAIL,
                            'tool' => $0,
                            'usehistory' => 'y',
                            'WebEnv' => "",            
                            );
            # search on data for accessions
            %params = esearch(%params);
            # Uncomment this section to try with epost
            #my %epost_result = epost_uids(%params);
            #$params{'WebEnv'} = $epost_result{'WebEnv'};
            #$params{'query_key'} = $epost_result{'query_key'};
            #$params{count} = $epost_result{count};
            # end section for epost
            print_OUT("Identified [ " . $params{count} . " ] entries.");
            next if ($params{count} == 0);
            # define parameters needed for efetch
            $params{num} = $params{count};
            $params{outfile} = "$0.$$.tmp.sequence.file.txt";
            $params{retmode} = 'text';
            $params{rettype} = 'gb';
            $params{batch} = $ncbi_entrez_batch_size;
            
            # get the sequence files.
            efetch_batch(%params);
            my %matched_ids = ();
            my $seqio = Bio::SeqIO->new( -file => $params{outfile}, -format => 'genbank' );
            while( my $seq = $seqio->next_seq ) {
		my $acc  = $seq->accession();
		next if (not exists $accs_to_gi{ $acc }->{'gi'});
                my %features = ();
                for my $feat_object ($seq->get_SeqFeatures) {
                    if ($feat_object->has_tag("db_xref")){
                        my ($id) = $feat_object->get_tag_values("db_xref");
                        my ($type,$value) = split(/:/,$id);
                        $features{$type} = $value;
                    }
                }
                next if (not defined $features{'taxon'});
                my $taxid  = $features{'taxon'};
                if (not defined $lineages{ $taxid }){
                    my @sbjct_lineage = reverse $seq->species->classification;
                    my $lca = get_lca_from_lineages(\@sbjct_lineage,\@ql); # need double checking on the ones that do not give match
                    next if ($lca eq "diff_root"); # exclude those that have a different root to cell organisms.
                    $gi_taxData{ $accs_to_gi{ $acc }->{'gi'} }->{'lineage'} = \@sbjct_lineage;
                    $gi_taxData{ $accs_to_gi{ $acc }->{'gi'} }->{'lca_with_qry'} = $lca;
                    $lineages{ $taxid }->{'lineage'} =  \@sbjct_lineage;
                    $lineages{ $taxid }->{'lca_with_qry'} = $lca;
                } else {
                    $gi_taxData{ $accs_to_gi{ $acc }->{'gi'} }->{'lineage'} = $lineages{ $taxid }->{'lineage'};
                    $gi_taxData{ $accs_to_gi{ $acc }->{'gi'} }->{'lca_with_qry'} = $lineages{ $taxid }->{'lca_with_qry'};
                }
                # remove it from the list of ids not solved
                delete( $ids_not_found{$accs_to_gi{ $acc }->{'full_id'} }  );
            }
            unlink($params{outfile});
        }
        if (scalar (keys %ids_not_found) > 0){
            if ($try_again_eutils < $max_eutils_queries) {
                $try_again_eutils++;
                goto(TRYAGAIN) if ($last_unmatched == 0);
                if (scalar (keys %ids_not_found) < $last_unmatched){
                    $last_unmatched = scalar (keys %ids_not_found);
                    goto(TRYAGAIN);
                    sleep(20);
                } else {
                    goto(PRINTUNMATCHED);
                }
            } else {
PRINTUNMATCHED:
                print_OUT("      ... there are still [ " . scalar (keys %ids_not_found) . " ] without match on DB. They will be printed to [ $out.unmatched_ids.txt ]");
                open (BADIDS,">$out.unmatched_ids.txt") or die $!;
                print BADIDS "UNMATCHED_IDS\tQUERY_IDS\n";
                while (my ($id,$qry_ids) = each %ids_not_found) {
                    print BADIDS "$id\t", join ",",@{$qry_ids};
                    print BADIDS "\n";
                }
                close(BADIDS);
            }
        }
        print_OUT("   '-> done with NCBI webservices.")
    }
}

############## PROCEED TO FINISH CALCULATIONS ##################
# Get the oldest stratum for every gene
my %qry_ancestors_ranks = ();
my $c = 0;
map { $qry_ancestors_ranks{$_} = $c++;   } @ql;
my %PhyloStratum_scores = ();
foreach my $qry_seq (keys %S) {
    my $oldest_stratum = scalar @ql - 1;
    foreach my $target_seq ( @{ $S{ $qry_seq } } ){
        #print Dumper($target_seq);
        #        getc;
        if (ref($target_seq) eq 'ARRAY'){
            if (scalar @{$target_seq} == 0){
                goto(WITHOUT_HITS);
            }
        }
        next if ($target_seq->{'pass_hard_thresold'} == 0);
        next if (not exists $gi_taxData{ $target_seq->{'subject_id'}}); # this ids have been printed to a file already.
            
        my $lca = $gi_taxData{ $target_seq->{'subject_id'} }->{'lca_with_qry'};
        next if ($lca eq "diff_root");
        my $target_seq_stratum = $qry_ancestors_ranks{ $lca };
        if ($target_seq_stratum < $oldest_stratum) {
            $oldest_stratum = $target_seq_stratum ;
        }
    }
WITHOUT_HITS:
    my @phyloScores = list zeroes scalar @ql;
    $phyloScores[ $oldest_stratum ] = 1;
    $PhyloStratum_scores{$qry_seq} = \@phyloScores;
}
print_OUT("Finished calculating hard coded scores");
print_OUT("Printing summary scores to [ $out.qry_node_phylostratumscores.txt ].");

# Print out the PhyloStratumScores as hardcoded
open(OUT,">$out.qry_node_phylostratumscores.txt") or die $!;
print OUT join "\t", @ql;
print OUT "\n";
print OUT join "\t", list sumover mpdl values %PhyloStratum_scores;
print OUT "\n";
close(OUT);

print_OUT("Printing gene phylostratum scores to [ $out.txt ].");

# Make table with scores for each gene
my $phyloScoresTable = join "\t", ("ID",@ql);
$phyloScoresTable .= "\n";

while (my ($qry_id,$qry_scores) = each %PhyloStratum_scores) {
    $phyloScoresTable .= join "\t", ($qry_id,@$qry_scores);
    $phyloScoresTable .= "\n";
}
open(OUT,">$out.txt") or die $!;
print OUT $phyloScoresTable;
close(OUT);

print_OUT("Done");

exit;

# TODO:
# 3- implement soft-coded
# 4- guess specie


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
    -soft, --soft_threshold    Uses a softhreshols strategy that aims to account for coverage and uncertaintiy of blast results.
    -use_coverage  Means that scores are weightes by coverage of sequence alignment. Not compatible with -hard.
    -no_ncbi_entrez     Do not use NCBI Entrez API to get information on sequence ids without data on local DBs.
    -ncbi_entrez_batch_size Number of IDs to submit to the NCBI API at time. DO NOT SET IT TO MORE THAN 500 (defualt 500).
    -guess_qry_specie   use blast result to guess query species.


 
 
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

=item B<-soft, --soft_threshold>
 
 Uses a softhreshols strategy that aims to account for coverage and uncertaintiy of blast results.

=item B<-use_coverage>
 
 Means that scores are weightes by coverage of sequence alignment. Not compatible with -hard.
 
=item B<-no_ncbi_entrez>
 
 Do not use NCBI Entrez API to get information on sequence ids without data on local DBs.

=item B<-ncbi_entrez_batch_size> 
 
 Number of IDs to submit to the NCBI API at time. DO NOT SET IT TO MORE THAN 500 (defualt 500).

=item B<-guess_qry_specie>
 
 Use blast result to guess query species. This done by selecting the most common specie with 100% identify blast hits. This (may) will only work if the specie of interest is actually on the DB use for the blast searches.

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

