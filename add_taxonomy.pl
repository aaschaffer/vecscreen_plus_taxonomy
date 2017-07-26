#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer and Eric Nawrocki
# Code to add taxonomy information to parsed vecscreen output for 
# vector screening
# Usage: add_taxonomy.pl --input_summary <input file> \ [REQUIRED]
#                        --input_taxa <taxonomy_file> \ [REQUIRED] 
#                        --outfile <output file>      \ [REQUIRED]
#                        --verbose                    \ [OPTIONAL]
#                        --keep                       \ [OPTIONAL]
#                        --debug                        [OPTIONAL]
use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";

# variables related to input/output files and command line options
my $input_summary_file;        # input file of parsed vecscreen output
my $input_taxa_file;           # input file with taxonomy in a tab-delimited four-column format
                               # (taxid, parent taxid, rank, depth where the depth of the root is 1)
my $output_file;               # output file that has the same rows and columns as input file with 
                               # one or two taxonomy columns added
my $verbose_mode;              # '1' if --verbose enabled, in which case extra columns are printed
my $keep_mode;                 # '1' if --keep enabled, in which case intermediate files are not removed
my $debug_mode;                # '1' if --debug enabled, in which case we output debugging print statements
my $temp_accession_file;       # file name for temporary accession file
my $temp_srcchk_file;          # file name for temporary srcchk output file
my $input_summary_FH;          # filehandle for reading $input_summary_file
my $output_FH;                 # filehandle for writing $output_file
my $output_accession_FH;       # filehandle for writing $temp_accession_file

my $nextline; # one line of input file
my $i;        #loop index
my %accession_H; # key: accession, value: 1 if we've seen this accession already

# hashes mapping accessions to taxonomy information
my %genus_H;     # key: accession, value: genus   (if known), or 1 otherwise
my %species_H;   # key: accession, value: species (if known), or 1 otherwise

# hashes mapping taxon to other information
my %taxonomy_parent_H;  # hash. key: taxon, value parent
my %taxonomy_rank_H;    # hash. key: taxon, value taxonomic rank (e.g., phylum)
my %taxonomy_level_H;   # hash. key: taxon, value level (a.k.a. depth) in the tree, where the level of the root is 1
                        # not currently used

# hard coded column counts 
my $NUM_COLUMNS_VERBOSE = 9; # number of columns if input is in verbose mode
my $NUM_COLUMNS_TERSE   = 4; # number of columns if input is not in verbose mode
my $num_input_columns;       # number of columns expected in input file, which will be 
                             # either $NUM_COLUMNS_VERBOSE or $NUM_COLUMNS_TERSE
# variables for debugging
my $genus_to_test = 325454;

# variables related to command line options, see epn-options.pm
my %opt_HH = ();
my @opt_order_A = ();
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
#       option                 type       default  group   requires incompat   preamble-outfile                                help-outfile
opt_Add("-h",                  "boolean", 0,           0,    undef, undef,     undef,                                          "display this help",                           \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"1"} = "required options";
opt_Add("--input_summary",     "string",  undef,       1,    undef, undef,     "input parsed vecscreen file",                  "REQUIRED: file name <s> of terminal summary", \%opt_HH, \@opt_order_A);
opt_Add("--input_taxa",        "string",  undef,       1,    undef, undef,     "input file mapping vecscreen matches to taxa", "REQUIRED: file name <s> of internal summary", \%opt_HH, \@opt_order_A);
opt_Add("--outfile",           "string",  undef,       1,    undef, undef,     "output file",                                  "REQUIRED: name <s> of output file to create", \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"2"} = "other options (not required)";
opt_Add("--verbose",           "boolean", 0,           2,    undef, undef,     "be verbose",                                   "be verbose in output",                        \%opt_HH, \@opt_order_A);
opt_Add("--keep",              "boolean", 0,           2,    undef, undef,     "keep intermediate files",                      "keep intermediate files",                     \%opt_HH, \@opt_order_A);
opt_Add("--debug",             "boolean", 0,           2,    undef, undef,     "debugging mode on",                            "turn debugging mode on",                      \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $all_options_recognized =
    &GetOptions('h'               => \$GetOptions_H{"-h"},
                'input_summary=s' => \$GetOptions_H{"--input_summary"},
                'input_taxa=s'    => \$GetOptions_H{"--input_taxa"},
                'outfile=s'       => \$GetOptions_H{"--outfile"},
                'verbose'         => \$GetOptions_H{"--verbose"},
                'keep'            => \$GetOptions_H{"--keep"},
                'debug'           => \$GetOptions_H{"--debug"});

my $synopsis = "add_taxonomy.pl :: add taxonomy information to parsed vecscreen output\n\n";
my $usage    = "Usage: perl add_taxonomy.pl ";

my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.01";
my $releasedate   = "Jan 2017";

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# define file names
$input_summary_file  = opt_Get("--input_summary",  \%opt_HH);
$input_taxa_file     = opt_Get("--input_taxa",     \%opt_HH);
$output_file         = opt_Get("--outfile",        \%opt_HH);
$verbose_mode        = opt_Get("--verbose",        \%opt_HH);
$keep_mode           = opt_Get("--keep",           \%opt_HH);
$debug_mode          = opt_Get("--debug",          \%opt_HH);

# We die if any of: 
# - non-existent option is used
# - any of the required options are not used. 
# - -h is used
my $reqopts_errmsg = "";
if(! defined $input_summary_file) { $reqopts_errmsg .= "ERROR, --input_summary not used. It is required.\n"; }
if(! defined $input_taxa_file)    { $reqopts_errmsg .= "ERROR, --input_taxa option not used. It is required.\n"; }
if(! defined $output_file)        { $reqopts_errmsg .= "ERROR, --outfile option not used. It is required.\n"; }

if($GetOptions_H{"-h"}) {
  opt_OutputHelp(*STDOUT, $synopsis . $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  exit 0;
}
if(($reqopts_errmsg ne "") || (! $all_options_recognized)) { 
  if   ($reqopts_errmsg ne "") { die $reqopts_errmsg; }
  else                         { die "ERROR, unrecognized option;"; }
}
# define the file names that depend on $output_file, now that we
# know that $output_file is valid
$temp_accession_file = $output_file . ".tmp_accession.txt";
$temp_srcchk_file    = $output_file . ".tmp_srcchk.txt";

# open the input and output files
open($input_summary_FH,    "<", "$input_summary_file")  or die "Cannot open $input_summary_file for reading\n"; 
open($output_FH,           ">", "$output_file")         or die "Cannot open $output_file for writing\n"; 
open($output_accession_FH, ">", "$temp_accession_file") or die "Cannot open $temp_accession_file for writing\n"; 

# process the taxonomy file
process_taxonomy_tree($input_taxa_file, \%taxonomy_parent_H, \%taxonomy_level_H, $debug_mode);

# read the input summary file and 
read_input($input_summary_FH, $output_accession_FH, \%accession_H, $verbose_mode, $debug_mode);
close($output_accession_FH); # important to close this here, call_srcchk_and_parse() may remove it

# call srcchk and parse its output
call_srcchk_and_parse($temp_accession_file, $temp_srcchk_file, $keep_mode);

# reset file handle to start of file
add_and_output_taxonomy_columns($input_summary_FH, $output_FH, \%genus_H, \%species_H, $verbose_mode, 
                                $NUM_COLUMNS_VERBOSE, $NUM_COLUMNS_TERSE);

# close file handles and exit
close($input_summary_FH);
close($output_FH);


#################################################################
# SUBROUTINES
#################################################################

#################################################################
# Subroutine: process_taxonomy_tree()
#
# Synopsis: Reads a file that includes NCBI's taxonomy information 
#           in four columns (taxon, parent taxon, rank, level).
#
# Args: $taxonomy_information_file: the NCBI taxonomy file
#       $taxonomy_parent_HR:        ref to hash, key is taxon, value is parent
#       $taxonomy_level_HR:         ref to hash, key is taxon, value is level

# Returns: nothing
#
# Dies: If file is in unexpected format: one or more lines does
#       not have exactly 4 columns.
# 
##################################################################
sub process_taxonomy_tree {
  my $sub_name = "process_taxonomy_tree()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($local_taxonomy_file, $taxonomy_parent_HR, $taxonomy_level_HR, $debug_mode) = @_;

  my $local_nextline;   # one line of taxonomy information
  my @local_fields_A;   # split of one line of taxonomy information

  # hard coded column meanings
  my $local_TAXID_COLUMN       = 0;
  my $local_PARENT_COLUMN      = 1;
  my $local_FORMAL_RANK_COLUMN = 2;
  my $local_LEVEL_COLUMN       = 3;

  open(TAXONOMY, "<", $local_taxonomy_file) or die "Cannot open $local_taxonomy_file for input in $sub_name\n";

  while(defined($local_nextline = <TAXONOMY>)) {
    chomp($local_nextline);
    @local_fields_A = split /\t/, $local_nextline;
    if(scalar(@local_fields_A) != 5) { die "ERROR in $sub_name, unexpected number of fields in $local_taxonomy_file in line: $local_nextline"; }
    $taxonomy_parent_H{$local_fields_A[$local_TAXID_COLUMN]} = $local_fields_A[$local_PARENT_COLUMN];
    $taxonomy_level_H{$local_fields_A[$local_TAXID_COLUMN]}  = $local_fields_A[$local_LEVEL_COLUMN];
    $taxonomy_rank_H{$local_fields_A[$local_TAXID_COLUMN]}   = $local_fields_A[$local_FORMAL_RANK_COLUMN];

    # debugging print statement
    if($debug_mode && ($genus_to_test == $local_fields_A[0])) { 
      print "Testing taxon $genus_to_test\n";
    }
  }
  close(TAXONOMY);

  return;
}

#################################################################
# Subroutine: read_input()
#
# Synopsis: Read the input file of parsed vecscreen output and 
#           extract the accessions for taxonomy lookup; Accessions 
#           are printed once each in the file associated with 
#           $accession_FH.
#
# Args: $input_summary_FH:    file handle to read from
#       $output_accession_FH: file handle to output to 
#       $accession_HR:        ref to hash of accession info (key: accession, value: 1)
#       $verbose_mode:        '1' if we're in verbose mode, else '0'
#       $debug_mode:          '1' if we're debugging (print extra info)
#
# Returns: nothing
# 
#################################################################
sub read_input {
  my $sub_name = "read_input()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($input_summary_FH, $output_accession_FH, $accession_HR, $verbose_mode, $debug_mode) = @_;
  my $local_nextline;  # one line with a single vecscreen match
  my @local_fields_A;  # split of $local_nextline
  my $local_accession; # one accession
  my $local_count = 0; # count of accessions
  
  if ($debug_mode) { 
    print "#DEBUG: In read_input\n";
  }
  while(defined($local_nextline = <$input_summary_FH>)) {
    # example line (columns are separated by tabs):
    #XM_010492096.1	36	258	uv|KF718979.1:2206-6729	2087	2309	Strong	Strong	Suspect[1,35];Suspect[259,265];
    chomp($local_nextline);
    @local_fields_A = split /\t/, $local_nextline;
    $local_accession = $local_fields_A[0];
    # check if accession has been seen before to avoid printing it 
    # twice to srcchk input file
    if ($debug_mode)  { 
      print "#DEBUG: Testing accession $local_accession\n";
    }
    if (!defined($accession_HR->{$local_accession})) {
      $accession_HR->{$local_accession} = 1;
      print $output_accession_FH "$local_accession\n";
      $local_count++;
    }
  }
  if ($debug_mode) { 
    print "#DEBUG: Ended read_input with $local_count accessions\n";
  }

  return;
}

###########################################################
# Subroutine: call_srcchk_and_parse()
#
# Synopsis: call srcchk to find the taxonomy information 
# for the accessions in $temp_accession_file.
#
# Args: $temp_accession_file: file name for temporary accession file
#                             (already created)
#       $temp_srcchk_file     file name for temporary srcchk output file 
#                             (created here)
#       $keep_mode:           '1' if we should keep intermediate files
#                             '0' if we should delete them
# Returns: nothing
#
###########################################################
sub call_srcchk_and_parse {
  my $sub_name = "call_srcck_and_parse()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($temp_accession_file, $temp_srcchk_file) = @_;

  # hard coded column info in srcchk output file
  my $local_ACCESSION_COLUMN = 0;   
  my $local_TAXID_COLUMN     = 1;
  my $local_nextline;       # one line with information about one taxid
  my @local_fields_A;       # split of $local_nextline 
  my $local_species;        # species for one accession
  my $local_genus;          # genus for one accession

  run_command("srcchk -i $temp_accession_file -f TaxID -o $temp_srcchk_file");

  # now parse the srcchk outfile that was just created
  open(TAXINFO, "<", $temp_srcchk_file) or die "Cannot open $temp_srcchk_file for input in $sub_name\n";

  <TAXINFO>; #skip header line
  # example of header line
  #accession	taxid	
  while(defined($local_nextline = <TAXINFO>)) {
    # example line: 
    #XM_010492096.1	90675	
    chomp($local_nextline);
    @local_fields_A = split /\t/, $local_nextline;    
    $local_species  = find_ancestor($local_fields_A[$local_TAXID_COLUMN], "species");
    $local_genus    = find_ancestor($local_species, "genus");
    $species_H{$local_fields_A[$local_ACCESSION_COLUMN]} = $local_species;
    $genus_H{$local_fields_A[$local_ACCESSION_COLUMN]}   = $local_genus;
  }
  close(TAXINFO);

  # remove accesion file, unless $keep_mode
  if(! $keep_mode) { 
    run_command("rm -f $temp_accession_file");
    run_command("rm -f $temp_srcchk_file");
  }
}

###########################################################
# Subroutine: add_and_output_taxonomy_columns() 
#
# Synopsis: Prints the input matches with one or two taxonomy 
# columns added after the first column. 
#
# Args: $input_summary_FH:    file handle to read from
#       $output_FH:           file handle to output to
#       $genus_HR:            ref to hash, key is accession, value is genus 
#                             (if known), or 1 otherwise
#       $species_HR:          ref to hash, key is accession, value is species
#                             (if known), or 1 otherwise
#       $verbose_mode:        '1' to output extra columns
#       $NUM_COLUMNS_VERBOSE: number of columns output if  $verbose_mode
#       $NUM_COLUMNS_TERSE:   number of columns output if !$verbose_mode
#
# Returns: nothing
#
###########################################################
sub add_and_output_taxonomy_columns {
  my $sub_name = "add_and_output_taxonomy_columns()";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  # rewind file so we can go back through it again
  seek $input_summary_FH, 0, 0; 
  
  my ($input_summary_FH, $output_FH, $genus_HR, $species_HR, $verbose_mode, 
      $NUM_COLUMNS_VERBOSE, $NUM_COLUMNS_TERSE) = @_;
  
  my $local_nextline;  # one line with a single vecscreen match
  my @local_fields_A;  # split of $local_nextline
  my $local_accession; # one accession
  
  while(defined($local_nextline = <$input_summary_FH>)) {
    chomp($local_nextline);
    @local_fields_A  = split /\t/, $local_nextline;
    $local_accession = $local_fields_A[0];

    print $output_FH "$local_accession\t";
    if (defined($genus_HR->{$local_accession})) {
      print $output_FH "$genus_HR->{$local_accession}\t";
    }

    # remainder of output differs depending on $verbose_mode
    if ($verbose_mode) {
      if (defined($species_HR->{$local_accession})) {
        print $output_FH "$species_HR->{$local_accession}\t";
      }
      for($i = 1; $i < $NUM_COLUMNS_VERBOSE; $i++) {
        print  $output_FH "$local_fields_A[$i]\t";
      }
    }
    else {
      for($i = 1; $i < $NUM_COLUMNS_TERSE; $i++) {
        print  $output_FH "$local_fields_A[$i]\t";
      }
    }
    print $output_FH "\n";
  }
  return;
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
    while ((1 != $local_ancestor) && (defined($taxonomy_parent_H{$local_ancestor})) && ((!(defined($taxonomy_rank_H{$local_ancestor})))  || (!($taxonomy_rank_H{$local_ancestor} eq $local_target_rank)))) {

        $local_ancestor = $taxonomy_parent_H{$local_ancestor};
    }
    if ((1 != $local_ancestor) && (!defined($taxonomy_parent_H{$local_ancestor}))) {
	return(1);
    }
    return($local_ancestor);
}

#################################################################
# Subroutine:  run_command()
# Incept:      EPN, Mon Dec 19 10:43:45 2016
#
# Purpose:     Runs a command using system() and exits in error 
#              if the command fails. 
#
# Arguments:
#   $cmd:         command to run, with a "system" command;
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if $cmd fails
#################################################################
sub run_command {
  my $sub_name = "run_command()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd) = @_;

  system($cmd);

  if($? != 0) { 
    die "ERROR in $sub_name, the following command failed:\n$cmd\n";
  }

  return;
}
