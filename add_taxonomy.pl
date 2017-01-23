#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer
# Code to add taxonomy information to parsed vecscreen output for vector screening
# Usage: add_taxonomy.pl --input_summary <input file> --input_taxa <taxonomy_file> (optional) --verbose --outfile <output file>

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";

my $input_summary_file; #input file of parsed vecscreen output
my $input_taxa_file;           #input file with taxonomy in a tab-delimited four-column format
                               #(taxid, parent taxid, rank, depth where the depth of the root is 1)
my $output_file; #output file that has the same rows and columns as input file with one or two taxonomy columns added
my $output_filehandle; #filehandle for output_file
my $temp_accession_file = "temp_accession_file.txt"; #temporary file holding accessions as input for call to srcchk
my $accession_filehandle; #filehandle for $temp_accession_file
my $nextline; #one line of input file
my $i; #loop index
my %accession_hash; #which accessions have we seen already

# hashes mapping accessions to taxonomy information

my %genus_hash; #maps an accession to its genus if known, or 1 otherwise
my %species_hash; #maps an accession to its species if known, or 1 otherwise

# hashes mapping taxon to other information
my %taxonomy_parent;              # hash, maps a taxon to its parent
my %taxonomy_level;               # hash, maps a taxon to its level (a.k.a. depth) in the tree, where the level of the root is 1
my %rank_hash;                    # hash, maps each taxon to its taxonomic rank (e.g., phylum)

my $verbose_mode; #did user specify verbose mode, in which case extra columns are printed
my $NUM_COLUMNS_VERBOSE = 9; #number of columns if input is in verbose mode
my $NUM_COLUMNS_TERSE = 4; #number of columns if input is in verbose mode
my $num_input_columns; #number of columns expected in input file, which will be either $NUM_COLUMNS_VERBOSE or $NUM_COLUMNS_TERSE

#variables for debugging
my $DEBUG = 0;
my $genus_to_test     = 325454;


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
opt_Add("--input_summary",           "string", undef,                        1,    undef, undef,     "input parsed vecscreen file",  "File name <s> with vecscreen output",     \%opt_HH, \@opt_order_A);
opt_Add("--input_taxa",                "string", undef,      1,    undef, undef,   "input file mapping vescreen matches to taxa",                                  "File name <s> mapping vecscreen matches to taxa",         \%opt_HH, \@opt_order_A);
opt_Add("--verbose",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                   "be verbose in output",  \%opt_HH, \@opt_order_A);
opt_Add("--outfile",           "string", undef,                        1,    undef, undef,     "output of sequences with terminal matches",  "File name <s> to hold output sequences, with terminal matches",     \%opt_HH,    \@opt_order_A);


# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "add_taxonomy.pl: add taxonomy information to parsed vecscreen output \n";
$usage    .= "add_taxonomy.pl --input_summary <input parsed vecscreen file> --input_taxa <input file with NCBI's taxonomy> (optional) --verbose --outfile <output> \n";
$usage   .= "\nFor example:\nadd_taxonomy.pl --input_summary parsed_vecscreen_output.txt --input_taxa taxonomy_wlevels.txt --outfile parsed_vecscreen_output_wtaxonomy.txt\n";

my $options_okay =
    &GetOptions('h'            => \$GetOptions_H{"-h"},
# basic options
                'input_summary=s'            => \$GetOptions_H{"--input_summary"},
                'input_taxa=s'                 => \$GetOptions_H{"--input_taxa"},
                'verbose'            => \$GetOptions_H{"--verbose"},
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

$input_summary_file = opt_Get("--input_summary", \%opt_HH);
$input_taxa_file           = opt_Get("--input_taxa", \%opt_HH);
$verbose_mode = opt_Get("--verbose", \%opt_HH);
$output_file = opt_Get("--outfile", \%opt_HH);

if (!defined($input_summary_file)) {
    die("--input_summary was not specified; exiting\n");
} 

if (!defined($output_file)) {
    die("--outfile was not specified; exiting\n");
} 



process_taxonomy_tree($input_taxa_file);
open(INPUT, "<$input_summary_file") or die "Cannot open 1 $input_summary_file\n"; 
if ($DEBUG) {
    print "Opened $input_summary_file\n";
}
open($accession_filehandle, ">$temp_accession_file") or die "Cannot open 2 $temp_accession_file\n"; 
open($output_filehandle, ">$output_file") or die "Cannot open 3 $output_file\n"; 

read_input($accession_filehandle);
close($accession_filehandle);
call_srcchk_and_parse();
seek INPUT, 0, 0;
add_taxonomy_columns($verbose_mode);

close(INPUT);
close($output_filehandle);

# Subroutine: read_input()
# Synopsis: read the input file of parsed vecscreen output and extract the accessions for taxonomy lookup;
#           accessions are printed once each in the file associated with $accession_filehandle
#
# Args: $accession_filehandle, $verbose_mode
#
#
#
# Returns: nothing
sub read_input {
    my $sub_name = "read_input()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_accession_filehandle) = @_;
    my $local_nextline; #one line with a single vecscreen match
    my @local_fields;   #split of $local_nextline
    my $local_accession; #one accession
    my $local_count = 0; #count of accessions

    if ($DEBUG) {
	print "In read_input\n";
    }
    while(defined($local_nextline = <INPUT>)) {
	chomp($local_nextline);
        @local_fields = split /\t/, $local_nextline;
	$local_accession = $local_fields[0];
	#check if accession has been seen before to avoid printing it twice to srcchk input file
	if ($DEBUG) {
	    print "Testing accession $local_accession\n";
	}
	if (!defined($accession_hash{$local_accession})) {
	    $accession_hash{$local_accession} = 1;
	    print $local_accession_filehandle "$local_accession\n";
	    $local_count++;
	}
    }
    if ($DEBUG) {
	print "Ended read_input with $local_count accessions\n";
    }
}


# Subroutine: call_srcchk_and_parse()
# Synopsis: call srcchk to find the taxonomy information for the accessions in $temp_accession_file
#
# Args: none
#
#
#
# Returns: nothing
sub call_srcchk_and_parse {
    my $sub_name = "call_srcck_and_parse()";
    my $nargs_exp = 0;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my $srcchk_output_file = "srcchk_output.txt"; #temporary file holding the results from the call to srcchk 
    my $local_ACCESSION_COLUMN   = 0;   
    my $local_TAXID_COLUMN       = 1;
    my $local_nextline; #one line with information about one taxid
    my @local_fields; #spli of $local_nextline 
    my $local_species; #species for one accession
    my $local_genus; #genus for one accession

    system("srcchk -i $temp_accession_file -f TaxID -o $srcchk_output_file");
    open(TAXINFO, "<", $srcchk_output_file) or die "Cannot open $srcchk_output_file for input in $sub_name\n";

    <TAXINFO>; #skip header line
    while(defined($local_nextline = <TAXINFO>)) {
        chomp($local_nextline);
        @local_fields = split /\t/, $local_nextline;    
	$local_species = find_ancestor($local_fields[$local_TAXID_COLUMN], "species");
	$local_genus = find_ancestor($local_species, "genus");
	$species_hash{$local_fields[$local_ACCESSION_COLUMN]} = $local_species;
	$genus_hash{$local_fields[$local_ACCESSION_COLUMN]} = $local_genus;
    }
    close(TAXINFO);
    system("rm -f $temp_accession_file");
}


# Subroutine: add_taxonomy_columns() 
# Synopsis: print the input matches with one or two taxonomy columns added after the first column
#
# Args: $local_verbose_mode
#             
#             
#
# Returns: nothing
sub add_taxonomy_columns {
    my $sub_name = "add_taxonomy_columns()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 

    my ($local_verbose_mode) = @_;
    my $local_nextline; #one line with a single vecscreen match
    my @local_fields;   #split of $local_nextline
    my $local_accession; #one accession



    while(defined($local_nextline = <INPUT>)) {
	chomp($local_nextline);
        @local_fields = split /\t/, $local_nextline;
	$local_accession = $local_fields[0];
	print  $output_filehandle "$local_accession\t";
	if (defined($genus_hash{$local_accession})) {
	    print  $output_filehandle "$genus_hash{$local_accession}\t";
	}
	if ($local_verbose_mode) {
	    if (defined($species_hash{$local_accession})) {
		print  $output_filehandle "$species_hash{$local_accession}\t";
	    }
	    for($i = 1; $i < $NUM_COLUMNS_VERBOSE; $i++) {
		print  $output_filehandle "$local_fields[$i]\t";
	    }
	}
	else {
	    for($i = 1; $i < $NUM_COLUMNS_TERSE; $i++) {
		print  $output_filehandle "$local_fields[$i]\t";
	    }
	}
	print $output_filehandle "\n";
    }
}

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

    open(TAXONOMY, "<", $local_taxonomy_file) or die "Cannot open $local_taxonomy_file for input in $sub_name\n";

    while(defined($local_nextline = <TAXONOMY>)) {
        chomp($local_nextline);
        @local_fields = split /\t/, $local_nextline;
        $taxonomy_parent{$local_fields[$local_TAXID_COLUMN]} = $local_fields[$local_PARENT_COLUMN];
        $taxonomy_level{$local_fields[$local_TAXID_COLUMN]}  = $local_fields[$local_LEVEL_COLUMN];
        if ($DEBUG && ($genus_to_test == $local_fields[0])) {
            print "Testing taxon $genus_to_test\n";
        }
        $rank_hash{$local_fields[$local_TAXID_COLUMN]} = $local_fields[$local_FORMAL_RANK_COLUMN];
    }
    close(TAXONOMY);
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

    $local_ancestor = $local_taxon;
    while ((1 != $local_ancestor) && (defined($taxonomy_parent{$local_ancestor})) && ((!(defined($rank_hash{$local_ancestor})))  || (!($rank_hash{$local_ancestor} eq $local_target_rank)))) {
        $local_ancestor = $taxonomy_parent{$local_ancestor};
    }
    if ((1 != $local_ancestor) && (!defined($taxonomy_parent{$local_ancestor}))) {
	return(1);
    }
    return($local_ancestor);
}
