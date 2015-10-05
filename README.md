# Intro

This is a Perl implementation of the PhyloStratigraphy comparative genomics methods. This program will perform a PhyloStratigraphy analysis. It provides a implementation and extensions of the original methodology of Domazet-Loso et. al. (2003).
 
	Please see README for information on Perl module requirments and some additional information.
 
## Command line options

For information on commands and options do > perl phylostratiphy.pl -man

## Installation

The following Perl modules are needed as requierment beyond those distributed with perl:
- PDL 
- Bio::DB::Taxonomy # part of BioPerl distribution
- Bio::LITE::Taxonomy 
- Bio::LITE::Taxonomy::NCBI 
- Bio::LITE::Taxonomy::NCBI::Gi2taxid

Both modules are available at CPAN. For additional information on installing PDL see pdl.perl.org

## Contact

Comments, suggestions or complains should be addressed to Inti Pedroso.
Inti Pedroso: IntiPeroso at gmail dot com


## References
	
- Pedroso I, Brown MJF, Sumner S. Detecting gene innovations for phenotypic diversity across multiple genomes. http://arxiv.org/abs/1212.3827
	
- Domazet-Loso, T., and Tautz, D. (2003). An evolutionary analysis of orphan genes in Drosophila. Genome Res. 13, 2213-2219.
