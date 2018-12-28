#!/usr/bin/perl -w
# the first line of perl code has to be as above
#
# Author: Alejandro Schaffer
#
#
# Code to find the taxonomy ancestors of taxa at a specified level (e.g., order) 
# 
# 
#
# Usage: find_taxonomy_ncestors.pl \
#        --input_summary <input file> \
#        --input_taxa <input file with NCBI's taxonomy tree as a list of edges>  \
#        --input_level <taxonomy level for which we want the output> \
#        --outfile <output file> 

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";


# input/output file names
my $input_file;                #input file that summarizes vecscreen matches 
my $input_taxa_file;           #input file with taxonomy in a tab-delimited four-column format 
                               #(taxid, parent taxid, rank, depth where the depth of the root is 1)
my $input_level;               #input taxonomy level, e.g., order
my $output_file;               #output file with two columns added showing the number and the name of taxonomy ancestor at the specified level 

# variables used in processing input
my $nextline;  #one line of the input summary file
my $m;         #loop index
my @fields;    #entries in one row
my $accession; #accession part of an identifier  

# The following five column indices are for the $input_file with 
# vecscreen matches passed in as the argument --input_summary
my $ACC_COLUMN          = 0; #the column with the identifier
my $TAXID_COLUMN        = 1; #the column with the taxid of the sequence represented by the identifier
my $NAME_COLUMN         = 2; #the column with the name of the taxid

# special taxon constants
my $TAXONOMY_ROOT        = 1;     #taxon to use as a default when no other taxon fits the situation
my $BACTERIA_TAXON       = 2;     #taxon that is the root of all Bacteria

# variables used in constructing output
my $output_taxid;        # taxonomic ancestor of interest as a number
my $output_taxon_name;   # name of output_taxid

# data structures for keeping track of information on matches
my %genome_match_taxon;        # hash of results seen before


# hashes mapping taxon to other information
my %taxonomy_parent;              # hash, maps a taxon to its parent
my %taxonomy_level;               # hash, maps a taxon to its level (a.k.a. depth) in the tree, where the level of the root is 1
my %rank_hash;                    # hash, maps each taxon to its taxonomic rank (e.g., phylum)
my %belongs_in_bacteria;          # hash, maps each taxon to 0 or 1 depending on whether it belongs in bacteria (1)

# variables used for processing the main input file, the vecscreen summary file
my $i;                      # loop index
my $num_columns_input = 3;  # number of columns in input 
my $one_taxid;              # taxid (in the taxonomy tree) of one query sequence in vecscreen summary
my $one_name;               # name of taxid

my @taxonomy_ranks = ("superkingdom","kingdom","phylum","class","order","family","tribe","genus"); #these are the ranks for which matches to UniVec entries have been collected

my $num_taxonomy_ranks = scalar(@taxonomy_ranks); # size of @taxonomy_ranks
my $rank_index;        #loop index for iteration in @taxonomy_ranks
my $one_rank;          #entry in @taxonomy_ranks

# variables used in special debugging mode
my $DEBUG             = 0;  # set to 1 to enter debug mode

# variables for dealing with options (epn-options.pm)
my %opt_HH           = ();
my @opt_order_A      = ();
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
$opt_group_desc_H{"1"} = "basic options";
#       option                         type      default group   requires incompat preamble-outfile                                                                help-outfile
opt_Add("-h",                          "boolean",0,          0,    undef, undef,   undef,                                                                          "display this help",                                       \%opt_HH, \@opt_order_A);
opt_Add("--input_summary",             "string", undef,      1,    undef, undef,   "input file of identifiers",                                                    "File name <s> of identifiers",                      \%opt_HH, \@opt_order_A);
opt_Add("--input_taxa",                "string", undef,      1,    undef, undef,   "input file with NCBI's taxonomy",                                              "File name <s> with NCBI's taxonomy",         \%opt_HH, \@opt_order_A);
opt_Add("--input_level",               "string", undef,      1,    undef, undef,   "input level of desired ancestor",                                              "String <s> of desired taxonomy level", \%opt_HH, \@opt_order_A);
opt_Add("--outfile",                   "string", undef,      1,    undef, undef,   "output file to create",                                                        "Name <s> of output file to create",                       \%opt_HH, \@opt_order_A);


my $synopsis = "find_taxonomy_ancestors.pl: Find the ancestors of specifie taxids at a specified taxonomy level\n";
my $usage    = "Usage:\n\n"; 
$usage      .= "\tfind_taxonomy_ancestors.pl \\\n";
$usage      .= "\t--input_summary <input file of identifiers> \\\n"; 
$usage      .= "\t--input_taxa <input file with NCBI's taxonomy> \\\n";
$usage      .= "\t--input_level <desired taxonomy level> \\\n";
$usage      .= "\t--outfile <output file>\n\n";
$usage      .= "\nFor example:\n";
$usage      .= "find_taxonomy_ancestors.pl --input_summary test_accessions.txt "; 
$usage      .= "--input_taxa taxonomy.txt ";
$usage      .= "--input_level order ";
$usage      .= "--outfile output.txt\n";

my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.16";
my $releasedate   = "December 2018";

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay =
    &GetOptions('h'                            => \$GetOptions_H{"-h"},
                'input_summary=s'              => \$GetOptions_H{"--input_summary"},
                'input_taxa=s'                 => \$GetOptions_H{"--input_taxa"},
                'input_level=s'                => \$GetOptions_H{"--input_level"},
                'outfile=s'                    => \$GetOptions_H{"--outfile"});

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) {
    opt_OutputHelp(*STDOUT, $synopsis . $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
    if(! $options_okay) { die "ERROR, unrecognized option;"; }
    else                { exit 0; } # -h, exit with 0 status
}

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# define file names
$input_file        = opt_Get("--input_summary", \%opt_HH);
$input_taxa_file           = opt_Get("--input_taxa", \%opt_HH);
$input_level               = opt_Get("--input_level", \%opt_HH);
$output_file               = opt_Get("--outfile", \%opt_HH);

# die if any of the required options were not used
my $errmsg = undef;
if(! defined $input_file)        { $errmsg .= "ERROR, --input_summary option not used.\n"; }
if(! defined $input_taxa_file)           { $errmsg .= "ERROR, --input_taxa option not used.\n"; }
if(! defined $input_level)               { $errmsg .= "ERROR, --input_level option not used.\n"; }
if(! defined $output_file)               { $errmsg .= "ERROR, --outfile option not used.\n"; }
if(defined $errmsg) { 
  die $errmsg . "\n$usage\n";
}


# open output files
open(SUMMARY, "<", $input_file) or die "Cannot open 1 $input_file for input\n"; 
open(OUTPUT,  ">", $output_file)        or die "Cannot open 2 $output_file for output\n"; 

# process input files and store the relevant information 
process_taxonomy_tree($input_taxa_file);
recognize_rank($input_level);


######################################################################
#
# Process the summary input file, one line at a time. 
#
$nextline = <SUMMARY>; #skip header
while(defined($nextline = <SUMMARY>)) {
    chomp($nextline);
    # default initialization of three output fields
    $output_taxid = $TAXONOMY_ROOT;
    @fields = split /\t/, $nextline;
    $accession   = $fields[$ACC_COLUMN];
    $one_taxid   = $fields[$TAXID_COLUMN];
    $one_name = $fields[$NAME_COLUMN];

    $output_taxid = find_ancestor($one_taxid,$input_level);
    
    # output original information from input summary file ($input_file)
    for($i = 0; $i < $num_columns_input; $i++) {
	print OUTPUT "$fields[$i]\t";
    }
    print OUTPUT "$output_taxid";  
    print OUTPUT "\n";
}
close (SUMMARY);
close (OUTPUT);

################################################
# SUBROUTINES
################################################
# 
#
# List of subroutines:
# process_taxonomy_tree();
# 
# recognize_rank():              helper function for process_genome_match_info() 
#
# Other subroutines:
# find_ancestor(): return ancestor at a specified rank for a taxon
#
################################################
# Subroutine: process_taxonomy_tree()
# Synopsis: reads a file that includes NCBI's taxonomy information in four columns (taxon, parent taxon, rank, level)
#
# Args: $taxonomy_information_file
#
# Returns: nothing

sub process_taxonomy_tree {

    my $sub_name = "process_taxonomy_tree()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_taxonomy_file) = @_;
    my $local_nextline; #one line of taxonomy information
    my @local_fields;   #split of one line of taxonomy information
    my $local_TAXID_COLUMN       = 0;
    my $local_PARENT_COLUMN      = 1;
    my $local_FORMAL_RANK_COLUMN = 2;
    my $local_LEVEL_COLUMN       = 3;
    my $local_BACTERIA_COLUMN       = 4;

    open(TAXONOMY, "<", $local_taxonomy_file) or die "Cannot open $local_taxonomy_file for input in $sub_name\n"; 
    
    while(defined($local_nextline = <TAXONOMY>)) {
	chomp($local_nextline);
	@local_fields = split /\t/, $local_nextline;
	$taxonomy_parent{$local_fields[$local_TAXID_COLUMN]} = $local_fields[$local_PARENT_COLUMN];
	$taxonomy_level{$local_fields[$local_TAXID_COLUMN]}  = $local_fields[$local_LEVEL_COLUMN];
	$belongs_in_bacteria{$local_fields[$local_TAXID_COLUMN]}  = $local_fields[$local_BACTERIA_COLUMN];
        $rank_hash{$local_fields[$local_TAXID_COLUMN]} = $local_fields[$local_FORMAL_RANK_COLUMN];
    }
    close(TAXONOMY);
}


# Subroutine: recognize_rank()
# Synopsis: takes a candidate rank as its only argument;
#           returns 1 if the rank is recognized and 0 otherwise.
#
# Args: $rank
#
# Returns: 1 if the rank is recognized and 0 otherwise

sub recognize_rank {

    my $sub_name = "recognize_rank()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_rank) = @_;
    my $local_r;
    
    for($local_r = 0; $local_r < $num_taxonomy_ranks; $local_r++) {
	if ($local_rank eq $taxonomy_ranks[$local_r]) { 
	    return 1;
	}
    }
    printf "$local_rank is not a recognized taxonomy rank\n";
    return 0;
}


# Subroutine: find_ancestor()
# Synopsis: given a taxon, returns its ancestor that is at a specific rank such as "genus";
#           it returns 1 (for root) if there is no ancestor at the specified target rank
#
# Args: $local_taxon
#       $local_target_rank
#
# Returns: the ancestor taxon that is at the specified rank or 1 (for root)

sub find_ancestor {

    my $sub_name = "find_ancestor()";
    my $nargs_exp = 2;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_taxon, $local_target_rank) = @_;
    my $local_ancestor; #one line of taxonomy information

    if ((!defined($local_taxon)) || (!defined($taxonomy_parent{$local_taxon}))) {
        print STDERR "Taxid $local_taxon is not in the taxonomy tree, the tree likely needs to be updated\n";
        return(1)
    }
    $local_ancestor = $local_taxon;
    while ((1 != $local_ancestor) && ((!(defined($rank_hash{$local_ancestor}))) || (!($rank_hash{$local_ancestor} eq $local_target_rank)))) {
        $local_ancestor = $taxonomy_parent{$local_ancestor};
	if ($DEBUG && defined($rank_hash{$local_ancestor}) ) {
	    print "Setting local_ancestor to $local_ancestor of rank $rank_hash{$local_ancestor} \n"; 
	}
    }            
    return($local_ancestor);
}        


    
