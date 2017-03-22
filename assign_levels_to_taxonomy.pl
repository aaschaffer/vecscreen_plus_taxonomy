#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer
# Code to assign levels to nodes in the taxonomy tree
# The level of a node is one more than the number of edges between that node and the root (1); the root is at level 1.
# Usage: assign_levels_to_taxonomy.pl --input_taxonomy <input file with NCBI's taxonomy tree as a list of edges>  --outfile <taxonomy tree with level added as the rightmost column>

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";

my $input_taxa_file; #input file with taxa of vecscreen matches
my $output_file; #output summary file with two columns added showing the taxon of the sequence and whether it is a genome match
my $nextline; #one line of the file
my $line_index; #counter for lines;
my @fields; #entries in one row;
#The following column indices are for the taxonomy_files 
my $SELF_COLUMN =0; #the column with the taxid
my $PARENT_COLUMN =1; #the column with the parent
my $RANK_COLUMN = 2; #the column with the taxonomy rank
my $LEVEL_COLUMN = 3; #the column with the taxonomy level of this taxid (in output file only)
my $ARTIFICIAL_TAXON = 81077; #taxon to use for artificial sequences
my $TAXONOMY_ROOT = 1; #taxon to use as a default when no other taxon fits the situation
my %taxonomy_parent; #maps a taxon to its parent
my %taxonomy_rank; #taxonomic formal rank
my %taxonomy_level; #maps a taon to its level
my @taxids; #array of all taxa
my $i; #loop index
my $num_taxids; #number of taxids


my %opt_HH = ();
my @opt_order_A = ();
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
$opt_group_desc_H{"1"} = "basic options";
#       option          type       default               group   requires incompat   preamble-outfile   help-outfile
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,     undef,            "display this help",                  \%opt_HH, \@opt_order_A);
opt_Add("--input_taxa",           "string", undef,                        1,    undef, undef,     "input file mapping vescreen matches to taxa",  "File name <s> mapping vecscreen matches to taxa",     \%opt_HH, \@opt_order_A);
opt_Add("--outfile",           "string", undef,                        1,    undef, undef,     "output file",  "Name <s> of output file",     \%opt_HH, \@opt_order_A);

# Usage: assign_levels_to_taxonomy.pl --input_taxonomy <input file with NCBI's taxonomy tree as a list of edges>  --outfile <taxonomy tree with level added as the rightmost column>

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "assign_levels_to_taxonomy.pl: assign levels to all taxids in NCBI's taxonomy tree\n";
$usage   .= "Usage:  --input_taxa <input file with NCBI's taxonomy>   --outfile <output file>\n";
$usage   .= "\nFor example:\nassign_levels_to_taxonomy.pl  --input_taxa taxonomy_tree.txt --outfile taxonomy_wlevels.txt\n";

my $options_okay =
    &GetOptions('h'            => \$GetOptions_H{"-h"},
# basic options
                'input_taxa=s'            => \$GetOptions_H{"--input_taxa"},
                'outfile=s'            => \$GetOptions_H{"--outfile"});


# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) {
    opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
    if(! $options_okay) { die "ERROR, unrecognized option;"; }
    else                { exit 0; } # -h, exit with 0 status
}


# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

$input_taxa_file = opt_Get("--input_taxa", \%opt_HH);
$output_file = opt_Get("--outfile", \%opt_HH);

$num_taxids = process_taxonomy_tree($input_taxa_file);
print "Returned $num_taxids taxids\n";
assign_levels($num_taxids);

open(OUTPUT, ">$output_file") or die "Cannot open 2 $output_file\n"; 
for($i = 0; $i < $num_taxids; $i++) {
    print OUTPUT "$taxids[$i]\t$taxonomy_parent{$taxids[$i]}\t$taxonomy_rank{$taxids[$i]}\t$taxonomy_level{$taxids[$i]}\n";
}
close(OUTPUT);



# Subroutine: process_taxonomy_tree()
# Synopsis: reads a file that includes NCBI's taxonomy information in three columns (taxon,parent taxon,level)
#
# Args: $taxonomy_information_file
#       
#
#
# Returns: number of taxids

sub process_taxonomy_tree {

    my $sub_name = "process_taxonomy_tree()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_taxonomy_file) = @_;
    my $local_nextline; #one line of taxonomy information
    my @local_fields; #split of one line of taxonomy information
    my $local_num_taxids;

    open(TAXONOMY, "<$local_taxonomy_file") or die "Cannot open 1 $local_taxonomy_file\n"; 
    
    $local_num_taxids = 0;
    while(defined($local_nextline = <TAXONOMY>)) {
	chomp($local_nextline);
	@local_fields = split /\t/, $local_nextline;
        $taxids[$local_num_taxids] = $local_fields[$SELF_COLUMN];
	$taxonomy_parent{$local_fields[0]} = $local_fields[$PARENT_COLUMN];
	$taxonomy_rank{$local_fields[0]} = $local_fields[$RANK_COLUMN];
	$local_num_taxids++;
    }
    close(TAXONOMY);
    return($local_num_taxids);
}

# Subroutine: assign_levels()
# Synopsis: assign levels to the taxids, one level at a time
#
# Args: $num_taxids
#       
#
#
# Returns: nothing

sub assign_levels {

    my $sub_name = "assign_levels()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_num_taxids) = @_;
    my $local_i; #loop index
    my $local_new_assignment; #did we assign a new value in this round

    for($local_i = 0; $local_i < $local_num_taxids; $local_i++) {
	$taxonomy_level{$taxids[$local_i]} = 0;
    }

    $taxonomy_level{$TAXONOMY_ROOT} = 1;
    $local_new_assignment = 1;
    while ($local_new_assignment) {
	$local_new_assignment = 0;
	print "In while loop\n";
	for($local_i = 0; $local_i < $local_num_taxids; $local_i++) {
	    if (10239 == $taxids[$local_i]) {
		printf "Testing 10239 in for loop and its parent is $taxonomy_parent{$taxids[$local_i]} with level $taxonomy_level{$taxonomy_parent{$taxids[$local_i]}}\n"
	    }
	    if ((0 == $taxonomy_level{$taxids[$local_i]}) && (0 < $taxonomy_level{$taxonomy_parent{$taxids[$local_i]}})) {
		$taxonomy_level{$taxids[$local_i]} = $taxonomy_level{$taxonomy_parent{$taxids[$local_i]}} + 1;
		$local_new_assignment = 1;
	    }
	}
    }
}
