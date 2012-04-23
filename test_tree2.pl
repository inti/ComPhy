use Bio::TreeIO; 
use Bio::DB::Taxonomy;
use Bio::Phylo::IO 'parse';
use Data::Dumper;
use IO::String;

my $nodesfile = "~/NOBACKUP/DATA/MySoft/ComparaPhyloStratiphy/data/nodes.dmp";
my $namefile = "~/NOBACKUP/DATA/MySoft/ComparaPhyloStratiphy/data/names.dmp";
my $db = Bio::DB::Taxonomy->new(-source => 'flatfile'
                              -nodesfile => $nodesfile,
                              -namesfile => $namefile);

my @species_names = ("Bombus terrestris", "Apis mellifera");
my $tree = $db->get_tree(@species_names);
my $project = parse(
    '-format'     => 'newick',
    '-string'     => $tree,
    '-as_project' => 1,
);


print Dumper($project);

getc;
my @nodes = $tree->get_nodes();
my $lca = $tree->get_lca(-nodes => \@nodes);
print $lca->id, "\t",$lca->description,"\n";
my $out = new Bio::TreeIO(-fh => \*STDOUT, -format => 'newick');
$out->write_tree($tree);

print ">",$tree,"<";

