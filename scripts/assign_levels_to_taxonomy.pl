#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer
# Code to assign levels to nodes in the taxonomy tree and recognize which taxids are descendants of Bacteria
# The level of a node is one more than the number of edges 
# between that node and the root (1); the root is at level 1.
# Usage: assign_levels_to_taxonomy.pl \ 
#        --input_taxonomy <input file with NCBI's taxonomy tree as a list of edges> \ [REQUIRED]
#        --outfile <taxonomy tree with level added as the rightmost column>         \ [REQUIRED]

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";

# variables related to input/output files 
my $input_taxa_file; # input file with taxa of vecscreen matches
my $output_file;     # output summary file with two columns added showing
                     # the taxon of the sequence and whether it is a genome match

# the following column indices are for the taxonomy_files 
my $SELF_COLUMN      = 0; #the column with the taxid
my $PARENT_COLUMN    = 1; #the column with the parent
my $RANK_COLUMN      = 2; #the column with the taxonomy rank
my $LEVEL_COLUMN     = 3; #the column with the taxonomy level of this taxid (in output file only)
my $BACTERIA_TAXON = 2; #taxon that is the root of all bacteria
my $ARTIFICIAL_TAXON = 81077; # taxon to use for artificial sequences
my $TAXONOMY_ROOT    = 1; # taxon to use as a default when no other taxon fits the situation

# variables related to taxa
my %taxonomy_parent_H; # maps a taxon to its parent
my %taxonomy_rank_H;   # taxonomic formal rank
my %taxonomy_level_H;  # maps a taxon to its level
my %belongs_in_bacteria_H; #maps a taxon to 0 or 1
my @taxids_A;          # array of all taxa
my $i;                 # loop index
my $num_taxids;        # number of taxids

# variables related to command line options, see epn-options.pm
my %opt_HH = ();
my @opt_order_A = ();
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
#       option          type       default  group  requires incompat  preamble-outfile                               help-outfile
opt_Add("-h",           "boolean", 0,           0,    undef, undef,   undef,                                         "display this help",                                         \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"1"} = "required options";
opt_Add("--input_taxa", "string",  undef,       1,    undef, undef,   "input file mapping vescreen matches to taxa", "REQUIRED: file name <s> mapping vecscreen matches to taxa", \%opt_HH, \@opt_order_A);
opt_Add("--outfile",    "string",  undef,       1,    undef, undef,   "output file",                                 "REQUIRED: file name <s> of output file to create",          \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();

my $all_options_recognized =
    &GetOptions('h'            => \$GetOptions_H{"-h"},
                'input_taxa=s' => \$GetOptions_H{"--input_taxa"},
                'outfile=s'    => \$GetOptions_H{"--outfile"});

my $synopsis = "assign_levels_to_taxonomy.pl :: assign levels to all taxids in NCBI's taxonomy tree\n\n";
my $usage    = "Usage: perl assign_levels_to_taxonomy.pl ";

my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.16";
my $releasedate   = "December 2018";

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);


# define file names
$input_taxa_file = opt_Get("--input_taxa", \%opt_HH);
$output_file     = opt_Get("--outfile",    \%opt_HH);

# We die if any of: 
# - non-existent option is used
# - any of the required options are not used. 
# - -h is used
my $reqopts_errmsg = "";
if(! defined $input_taxa_file) { $reqopts_errmsg .= "ERROR, --input_taxa not used. It is required.\n"; }
if(! defined $output_file)     { $reqopts_errmsg .= "ERROR, --outfile option not used. It is required.\n"; }

if($GetOptions_H{"-h"}) {
  opt_OutputHelp(*STDOUT, $synopsis . $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  exit 0;
}
if(($reqopts_errmsg ne "") || (! $all_options_recognized)) { 
  if   ($reqopts_errmsg ne "") { die $reqopts_errmsg; }
  else                         { die "ERROR, unrecognized option;"; }
}

# process the taxonomy file
$num_taxids = process_taxonomy_tree($input_taxa_file);
# print "Returned $num_taxids taxids\n";

# assign levels
assign_levels($num_taxids);

# create output file
open(OUTPUT, ">$output_file") or die "Cannot open 2 $output_file\n"; 
for($i = 0; $i < $num_taxids; $i++) {
    print OUTPUT "$taxids_A[$i]\t$taxonomy_parent_H{$taxids_A[$i]}\t$taxonomy_rank_H{$taxids_A[$i]}\t$taxonomy_level_H{$taxids_A[$i]}\t$belongs_in_bacteria_H{$taxids_A[$i]}\n";
}
close(OUTPUT);

#################################################################
# SUBROUTINES
#################################################################

#################################################################
# Subroutine: process_taxonomy_tree()
#
# Synopsis: Reads a file that includes NCBI's taxonomy information 
#           in four columns (taxon, parent taxon, rank, level).
#
# THIS SUBROUTINE IS SIMILAR BUT NOT IDENTICAL TO ONE OF THE 
# SAME NAME IN add_taxonomy.pl AND compare_vector_matches_wtaxa.pl
#
# Args: $taxonomy_information_file: the NCBI taxonomy file
#
# Returns: Number of taxids read/processed.
#
# Dies: If file is in unexpected format: one or more lines does
#       not have exactly 5 columns.
# 
##################################################################
sub process_taxonomy_tree {
  my $sub_name = "process_taxonomy_tree()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($local_taxonomy_file) = @_;
  my $local_nextline; #one line of taxonomy information
  my @local_fields_A; #split of one line of taxonomy information
  my $local_num_taxids;
  
  open(TAXONOMY, "<$local_taxonomy_file") or die "Cannot open 1 $local_taxonomy_file\n"; 
  
  $local_num_taxids = 0;
  while(defined($local_nextline = <TAXONOMY>)) {
    chomp($local_nextline);
    @local_fields_A = split /\t/, $local_nextline;
    $taxids_A[$local_num_taxids] = $local_fields_A[$SELF_COLUMN];
    $taxonomy_parent_H{$local_fields_A[0]} = $local_fields_A[$PARENT_COLUMN];
    $taxonomy_rank_H{$local_fields_A[0]} = $local_fields_A[$RANK_COLUMN];
    $belongs_in_bacteria_H{$local_fields_A[0]} = 0;
    $local_num_taxids++;
  }
  $belongs_in_bacteria_H{$BACTERIA_TAXON} = 1;
  close(TAXONOMY);
  return($local_num_taxids);
}

#################################################################
# Subroutine: assign_levels()
#
# Synopsis: Assign levels to the taxids, one level at a time.
#
# Args: $num_taxids: number of taxids in $taxonomy_level hash.
#
# Returns: nothing.
#
# Dies: Never.
# 
##################################################################
sub assign_levels {
  my $sub_name = "assign_levels()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($local_num_taxids) = @_;
  my $local_i; #loop index
  my $local_new_assignment; #did we assign a new value in this round
  
  for($local_i = 0; $local_i < $local_num_taxids; $local_i++) {
    $taxonomy_level_H{$taxids_A[$local_i]} = 0;
  }
  
  $taxonomy_level_H{$TAXONOMY_ROOT} = 1;
  $local_new_assignment = 1;
  while ($local_new_assignment) {
    $local_new_assignment = 0;
    # print "In while loop\n";
    for($local_i = 0; $local_i < $local_num_taxids; $local_i++) {
      # if (10239 == $taxids_A[$local_i]) {
      # printf "Testing 10239 in for loop and its parent is $taxonomy_parent_H{$taxids_A[$local_i]} with level $taxonomy_level{$taxonomy_parent_H{$taxids_A[$local_i]}}\n"
      # }
      if ((0 == $taxonomy_level_H{$taxids_A[$local_i]}) && (0 < $taxonomy_level_H{$taxonomy_parent_H{$taxids_A[$local_i]}})) {
        $taxonomy_level_H{$taxids_A[$local_i]} = $taxonomy_level_H{$taxonomy_parent_H{$taxids_A[$local_i]}} + 1;
	if ($belongs_in_bacteria_H{$taxonomy_parent_H{$taxids_A[$local_i]}}) {
	    $belongs_in_bacteria_H{$taxids_A[$local_i]} = 1;
	}
        $local_new_assignment = 1;
      }
    }
  }
}
