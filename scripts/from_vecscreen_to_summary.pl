#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Author: Alejandro Schaffer with help from Eric Nawrocki
#
# Code to go in several high-level steps from a set of accessions to a
# summary of parsed vecscreen output 
#
# Usage: from_vecscreen_to_summary.pl --input_fasta <input fasta file> \             [REQUIRED]
#                                --input_taxa <input file of taxonomy with levels> \ [REQUIRED]
#                                --output_root <root for naming output files> \      [REQUIRED]
#                                --verbose \                                         [OPTIONAL]
#                                --combine_output \                                  [OPTIONAL]
#                                --keep                                              [OPTIONAL]

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday); # for timings

require "epn-options.pm";

# make sure the VECPLUSDIR environment variable is set
my $vecplusdir = $ENV{'VECPLUSDIR'};
if(! exists($ENV{'VECPLUSDIR'})) { 
  printf STDERR ("\nERROR, the environment variable VECPLUSDIR is not set, please set it to the directory that was created when you unpackaged the vecscreen_plus_taxonomy package.\n"); 
  exit(1); 
}
if(! (-d $vecplusdir)) { 
  printf STDERR ("\nERROR, the vecscreen_plus_taxonomy directory specified by your environment variable VECPLUSDIR does not exist.\n"); 
  exit(1); 
}    
my $vecplus_exec_dir = $vecplusdir . "/scripts/"; 

my $output_root      = undef; # root for naming output files
my $input_fasta_file = undef; # input fasta file
my $input_taxa_file  = undef; # input file with taxonomy in a tab-delimited four-column format
                              # (taxid, parent taxid, rank, depth where the depth of the root is 1)
my $verbose_mode;             # did user specify verbose mode, in which case extra columns are printed
my $verbose_string;           # string to add as argument to called programs depending on $verbose_mode 
my $keep_mode;                # string to add as argument to called programs depending on $keep_mode 
my $keep_string;              # string to add as argument to called programs depending on $verbose_mode 
my $combine_summaries_mode;   # should the internal and terminal matches be combined into one file

my %opt_HH = ();
my @opt_order_A = ();
my %opt_group_desc_H = ();


# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
#       option              type       default group requires incompat preamble-outfile                                       help-outfile
opt_Add("-h",               "boolean", 0,         0,    undef, undef,  undef,                                                 "display this help",                                         \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"1"} = "required options";
opt_Add("--input_fasta",    "string",  undef,     1,    undef, undef,  "input fasta file",                                    "REQUIRED: file name <s> with sequences in fasta format",    \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"2"} = "other options (not required)";
opt_Add("--input_taxa",     "string",  undef,     1,    undef, undef,  "input file mapping vecscreen matches to taxa",        "REQUIRED: file name <s> mapping vecscreen matches to taxa", \%opt_HH, \@opt_order_A);
opt_Add("--output_root",    "string",  undef,     1,    undef, undef,  "output files will be named starting with",            "REQUIRED: output files will be named starting with <s>",    \%opt_HH, \@opt_order_A);
opt_Add("--verbose",        "boolean", 0,         2,    undef, undef,  "output 11 columns instead of 5",                      "output 11 columns instead of 5",                            \%opt_HH, \@opt_order_A);
opt_Add("--combine_output", "boolean", 0,         2,    undef, undef,  "combine internal and terminal matches",               "combine internal and terminal matches",                     \%opt_HH, \@opt_order_A);
opt_Add("--keep",           "boolean", 0,         2,    undef, undef,  "keep all intermediate files (e.g. vecscreen output)", "keep all intermediate files (e.g. vecscreen output)",       \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $all_options_recognized =
    &GetOptions('h'                  => \$GetOptions_H{"-h"},
                'output_root=s'      => \$GetOptions_H{"--output_root"},
                'input_fasta=s'      => \$GetOptions_H{"--input_fasta"},
                'input_taxa=s'       => \$GetOptions_H{"--input_taxa"},
                'verbose'            => \$GetOptions_H{"--verbose"},
                'combine_output'     => \$GetOptions_H{"--combine_output"},
                'keep'               => \$GetOptions_H{"--keep"});

my $synopsis = "from_vecscreen_to_summary.pl :: run vecscreen, parse its output and add taxonomy information";
my $usage    = "Usage: from_vecscreen_to_summary.pl ";

my $total_seconds = -1 * seconds_since_epoch(); # by multiplying by -1, we can just add another seconds_since_epoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.06";
my $releasedate   = "Aug 2017";

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# define file names
$output_root            = opt_Get("--output_root", \%opt_HH);
$input_fasta_file       = opt_Get("--input_fasta", \%opt_HH);
$input_taxa_file        = opt_Get("--input_taxa", \%opt_HH);
$verbose_mode           = opt_Get("--verbose", \%opt_HH);
$combine_summaries_mode = opt_Get("--combine_output", \%opt_HH);
$keep_mode              = opt_Get("--keep", \%opt_HH);

# exit if necessary (if options were messed up)
# first, determine if all required options were actually used
my $reqopts_errmsg = "";
if(! defined $output_root)        { $reqopts_errmsg .= "ERROR, --output_root option not used. It is required.\n"; }
if(! defined $input_fasta_file)   { $reqopts_errmsg .= "ERROR, --input_fasta option not used. It is required.\n"; }
if(! defined $input_taxa_file)    { $reqopts_errmsg .= "ERROR, --input_taxa option not used. It is required.\n"; }

if(($reqopts_errmsg ne "") || (! $all_options_recognized) || ($GetOptions_H{"-h"})) {
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if   ($reqopts_errmsg ne "")     { die $reqopts_errmsg; }
  elsif(! $all_options_recognized) { die "ERROR, unrecognized option;"; }
  else                             { exit 0; } # -h, exit with 0 status
}

$verbose_string = ($verbose_mode) ? "--verbose" : "";
$keep_string    = ($keep_mode)    ? "--keep"    : "";

# make sure the executable files we need exist in $VECPLUSDIR
my %execs_H = (); # hash with paths to all required executables
$execs_H{"vecscreen"}         = $vecplus_exec_dir . "vecscreen";
$execs_H{"srcchk"}            = $vecplus_exec_dir . "srcchk";
$execs_H{"parse_vecscreen"}   = $vecplus_exec_dir . "parse_vecscreen.pl";
$execs_H{"combine_summaries"} = $vecplus_exec_dir . "combine_summaries.pl";
$execs_H{"add_taxonomy"}      = $vecplus_exec_dir . "add_taxonomy.pl";
validate_executable_hash(\%execs_H);

# set output file names
my $temp_vecscreen_output_file     = $output_root . ".vecscreen_output.txt";
my $temp_internal_output_file      = $output_root . ".output_internal.txt";
my $temp_terminal_output_file      = $output_root . ".output_terminal.txt";
my $temp_combined_output_file      = $output_root . ".output_combined.txt";
my $combined_wtaxonomy_output_file = $output_root . ".output_combined_wtaxonomy.txt";
my $terminal_wtaxonomy_output_file = $output_root . ".output_terminal_wtaxonomy.txt";
my $internal_wtaxonomy_output_file = $output_root . ".output_internal_wtaxonomy.txt";

# output banner and preamble
my @arg_desc_A = (); # necessary to pass into opt_OutputPreamble()
my @arg_A      = (); # necessary to pass into opt_OutputPreamble()
output_banner(*STDOUT, $version, $releasedate, $synopsis, $date);
opt_OutputPreamble(*STDOUT, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

##########################
# Step 1. Call vecscreen #
##########################
my $progress_w = 51; # the width of the left hand column in our progress output, hard-coded
my $start_secs = output_progress_prior("Running vecscreen", $progress_w, undef, *STDOUT);
run_command($execs_H{"vecscreen"} . " -query $input_fasta_file -text_output > $temp_vecscreen_output_file", 0); # 0: don't echo command to STDOUT
my $desc_str = sprintf("output saved as $temp_vecscreen_output_file%s", $keep_mode ? "]" : " (temporarily)"); 
output_progress_complete($start_secs, $desc_str, undef, *STDOUT);

##################################
# Step 2. Parse vecscreen output #
##################################
$start_secs = output_progress_prior("Parsing vecscreen output", $progress_w, undef, *STDOUT);
run_command($execs_H{"parse_vecscreen"} . " --input $temp_vecscreen_output_file $verbose_string --outfile_internal $temp_internal_output_file --outfile_terminal $temp_terminal_output_file", 0); # 0: don't echo command to STDOUT
$desc_str = sprintf("output saved as $temp_internal_output_file and $temp_terminal_output_file%s", $keep_mode ? "]" : " (temporarily)"); 
output_progress_complete($start_secs, $desc_str, undef, *STDOUT);

###############################################
# Step 3. [OPTIONAL] Combine vecscreen output #
###############################################
if ($combine_summaries_mode) {
  $start_secs = output_progress_prior("Combining output [due to --combine_output]", $progress_w, undef, *STDOUT);
  run_command($execs_H{"combine_summaries"} . " --input_terminal $temp_terminal_output_file --input_internal $temp_internal_output_file $verbose_string --outfile $temp_combined_output_file", 0); # 0: don't echo command to STDOUT
  $desc_str = sprintf("output saved as $temp_combined_output_file%s", $keep_mode ? "]" : " (temporarily)"); 
  output_progress_complete($start_secs, $desc_str, undef, *STDOUT);
}

####################################
# Step 4. Add taxonomy information #
####################################
if($combine_summaries_mode) { 
  $start_secs = output_progress_prior("Adding taxonomy information to output", $progress_w, undef, *STDOUT);
  run_command($execs_H{"add_taxonomy"} . " --input_summary $temp_combined_output_file --input_taxa $input_taxa_file $verbose_string $keep_string --outfile $combined_wtaxonomy_output_file", 0); # 0: don't echo command to STDOUT
  $desc_str = "output saved as $combined_wtaxonomy_output_file";
  output_progress_complete($start_secs, $desc_str, undef, *STDOUT);
}
else {
  $start_secs = output_progress_prior("Adding taxonomy information to output in two stages", $progress_w, undef, *STDOUT);
  run_command($execs_H{"add_taxonomy"} . " --input_summary $temp_internal_output_file --input_taxa $input_taxa_file $verbose_string --outfile $internal_wtaxonomy_output_file", 0); # 0: don't echo command to STDOUT
  run_command($execs_H{"add_taxonomy"} . " --input_summary $temp_terminal_output_file --input_taxa $input_taxa_file $verbose_string --outfile $terminal_wtaxonomy_output_file", 0); # 0: don't echo command to STDOUT
  $desc_str = "output saved as $internal_wtaxonomy_output_file and $combined_wtaxonomy_output_file";
  output_progress_complete($start_secs, $desc_str, undef, *STDOUT);
}

####################
# Step 5. Clean up #
####################
if(! $keep_mode) { 
  $start_secs = output_progress_prior("Cleaning up temporary files", $progress_w, undef, *STDOUT);
  run_command("rm -f $temp_vecscreen_output_file $temp_internal_output_file $temp_terminal_output_file", 0); # 0: don't echo command to STDOUT
  $desc_str = "deleted $temp_vecscreen_output_file, $temp_internal_output_file";
  if($combine_summaries_mode) { 
    run_command("rm -f $temp_combined_output_file", 0); # 0: don't echo command to STDOUT
    $desc_str .= ", $temp_terminal_output_file and $temp_combined_output_file";
  }
  else { 
    $desc_str .= " and $temp_terminal_output_file";
  }
  output_progress_complete($start_secs, $desc_str, undef, *STDOUT);
}

#############################
# Step 6. Conclude and exit #
#############################
$total_seconds += seconds_since_epoch();
printf STDOUT ("#\n");
if($combine_summaries_mode) { 
  printf STDOUT ("# Final combined output file in %d column format saved to:\n#\t%s\n#\n", ($verbose_mode ? 11 : 5), $combined_wtaxonomy_output_file);
}
else { 
  printf STDOUT ("# Final internal match output file in %d column format saved to:\n#\t%s\n", ($verbose_mode ? 11 : 5), $internal_wtaxonomy_output_file);
  printf STDOUT ("# Final terminal match output file in %d column format saved to:\n#\t%s\n#\n", ($verbose_mode ? 11 : 5), $terminal_wtaxonomy_output_file);
}
output_column_descriptions(*STDOUT, $verbose_mode);
printf STDOUT ("# Total time: %.1f seconds\n", $total_seconds);
printf STDOUT ("# \n");
printf STDOUT ("# SUCCESS\n");

exit 0;


#####################################################################
# SUBROUTINES 
#####################################################################
# List of subroutines:
#
# Functions for output: 
# output_column_descriptions: output descriptions of columns in final output file
# output_banner:              output the banner with info on the script and options used
# output_progress_prior:      output routine for a step, prior to running the step
# output_progress_complete:   output routine for a step, after the running the step
#
# Miscellaneous functions:
# run_command:              run a command using system()
# seconds_since_epoch:      number of seconds since the epoch, for timings
#
#################################################################
# Subroutine:  output_column_descriptions()
# Incept:      EPN, Thu Jan  5 10:17:15 2017
#
# Purpose:     Output descriptions of all output columns 
#              to $FH. One line per column. If $verbose_mode
#              is '1', 11 column descriptions are output,
#              else 5 column descriptions are output.
#
# Arguments:
#   $FH:           file handle to output to
#   $verbose_mode: '1' if in verbose mode, '0' if default mode
# 
# Returns:    nothing
#
# Dies:       neve
#################################################################
sub output_column_descriptions {
  my $sub_name = "output_column_descriptions()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($FH, $verbose_mode) = @_;

  my $accn_desc   = "Accession of query";
  my $genus_desc  = "Genus of query if known, or 1 otherwise";
  my $vector_desc = "Matching vector, starting with uv|";
  my $end1_desc   = "One end of the alignment in the vector";
  my $end2_desc   = "The other end of the alignment in the vector";

  if($verbose_mode) { 
    print $FH "# Descriptions of the 11 columns:\n#\n";
    print $FH "# Column  1: $accn_desc\n";
    print $FH "# Column  2: $genus_desc\n";
    print $FH "# Column  3: Species of query if known, or 1 otherwise\n";
    print $FH "# Column  4: Lower end of the alignment in the vector\n";
    print $FH "# Column  5: Upper end of the alignment in the vector\n";
    print $FH "# Column  6: $vector_desc\n";
    print $FH "# Column  7: $end1_desc\n";
    print $FH "# Column  8: $end2_desc\n";
    print $FH "# Column  9: The strength of this vecscreen match \n";
    print $FH "# Column 10: The strength of the strongest vecscreen match for this query\n";
    print $FH "# Column 11: Whether there is any dangling part (called \"Suspect\" by vecscreen) at either end of the query\n#\n";
    print $FH "# A dangling part is an unmatched segment of <= 25 nucleotides.\n#\n"  
  }
  else { # default mode: 5 columns
    print $FH "# Descriptions of the 5 columns:\n#\n";
    print $FH "# Column 1: $accn_desc\n";
    print $FH "# Column 2: $genus_desc\n";
    print $FH "# Column 3: $vector_desc\n";
    print $FH "# Column 4: $end1_desc\n";
    print $FH "# Column 5: $end2_desc\n#\n";
  }

  return;
}

#################################################################
# Subroutine : output_progress_prior()
# Incept:      EPN, Fri Feb 12 17:22:24 2016 [dnaorg.pm]
#
# Purpose:      Output to $FH1 (and possibly $FH2) a message indicating
#               that we're about to do 'something' as explained in
#               $outstr.  
#
#               Caller should call *this* function, then do
#               the 'something', then call output_progress_complete().
#
#               We return the number of seconds since the epoch, which
#               should be passed into the downstream
#               output_progress_complete() call if caller wants to
#               output running time.
#
# Arguments: 
#   $outstr:     string to print to $FH
#   $progress_w: width of progress messages
#   $FH1:        file handle to print to, can be undef
#   $FH2:        another file handle to print to, can be undef
# 
# Returns:     Number of seconds and microseconds since the epoch.
#
################################################################# 
sub output_progress_prior { 
  my $nargs_expected = 4;
  my $sub_name = "output_progress_prior()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($outstr, $progress_w, $FH1, $FH2) = @_;

  if(defined $FH1) { printf $FH1 ("# %-*s ... ", $progress_w, $outstr); }
  if(defined $FH2) { printf $FH2 ("# %-*s ... ", $progress_w, $outstr); }

  return seconds_since_epoch();
}

#################################################################
# Subroutine : output_progress_complete()
# Incept:      EPN, Fri Feb 12 17:28:19 2016 [dnaorg.pm]
#
# Purpose:     Output to $FH1 (and possibly $FH2) a 
#              message indicating that we've completed 
#              'something'.
#
#              Caller should call *this* function,
#              after both a call to output_progress_prior()
#              and doing the 'something'.
#
#              If $start_secs is defined, we determine the number
#              of seconds the step took, output it, and 
#              return it.
#
# Arguments: 
#   $start_secs:    number of seconds either the step took
#                   (if $secs_is_total) or since the epoch
#                   (if !$secs_is_total)
#   $extra_desc:    extra description text to put after timing
#   $FH1:           file handle to print to, can be undef
#   $FH2:           another file handle to print to, can be undef
# 
# Returns:     Number of seconds the step took (if $secs is defined,
#              else 0)
#
################################################################# 
sub output_progress_complete { 
  my $nargs_expected = 4;
  my $sub_name = "output_progress_complete()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($start_secs, $extra_desc, $FH1, $FH2) = @_;

  my $total_secs = undef;
  if(defined $start_secs) { 
    $total_secs = seconds_since_epoch() - $start_secs;
  }

  if(defined $FH1) { printf $FH1 ("done."); }
  if(defined $FH2) { printf $FH2 ("done."); }

  if(defined $total_secs || defined $extra_desc) { 
    if(defined $FH1) { printf $FH1 (" ["); }
    if(defined $FH2) { printf $FH2 (" ["); }
  }
  if(defined $total_secs) { 
    if(defined $FH1) { printf $FH1 (sprintf("%.1f seconds%s", $total_secs, (defined $extra_desc) ? ", " : "")); }
    if(defined $FH2) { printf $FH2 (sprintf("%.1f seconds%s", $total_secs, (defined $extra_desc) ? ", " : "")); }
  }
  if(defined $extra_desc) { 
    if(defined $FH1) { printf $FH1 $extra_desc };
    if(defined $FH2) { printf $FH2 $extra_desc };
  }
  if(defined $total_secs || defined $extra_desc) { 
    if(defined $FH1) { printf $FH1 ("]"); }
    if(defined $FH2) { printf $FH2 ("]"); }
  }

  if(defined $FH1) { printf $FH1 ("\n"); }
  if(defined $FH2) { printf $FH2 ("\n"); }
  
  return (defined $total_secs) ? $total_secs : 0.;
}

#####################################################################
# Subroutine: output_banner()
# Incept:     EPN, Thu Oct 30 09:43:56 2014 (rnavore)
# 
# Purpose:    Output the banner with info on the script, input arguments
#             and options used.
#
# Arguments: 
#    $FH:                file handle to print to
#    $version:           version
#    $releasedate:       month/year of version (e.g. "Feb 2016")
#    $synopsis:          string reporting the date
#    $date:              date information to print
#
# Returns:    Nothing, if it returns, everything is valid.
# 
# Dies: never
####################################################################
sub output_banner {
  my $nargs_expected = 5;
  my $sub_name = "output_banner()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($FH, $version, $releasedate, $synopsis, $date) = @_;

  print $FH ("\# $synopsis\n");
  print $FH ("\# version: $version ($releasedate)\n");
  print $FH ("\# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  if(defined $date)    { print $FH ("# date:    $date\n"); }
  printf $FH ("#\n");

  return;
}

#################################################################
# Subroutine:  run_command()
# Incept:      EPN, Mon Dec 19 10:43:45 2016
#
# Purpose:     Runs a command using system() and exits in error 
#              if the command fails. If $be_verbose, outputs
#              the command to stdout. If $FH_HR->{"cmd"} is
#              defined, outputs command to that file handle.
#
# Arguments:
#   $cmd:         command to run, with a "system" command;
#   $be_verbose:  '1' to output command to stdout before we run it, '0' not to
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if $cmd fails
#################################################################
sub run_command {
  my $sub_name = "run_command()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $be_verbose) = @_;
  
  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));

  system($cmd);

  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));

  if($? != 0) { 
    die "ERROR in $sub_name, the following command failed:\n$cmd\n";
  }

  return ($stop_time - $start_time);
}

#################################################################
# Subroutine : seconds_since_epoch()
# Incept:      EPN, Sat Feb 13 06:17:03 2016
#
# Purpose:     Return the seconds and microseconds since the 
#              Unix epoch (Jan 1, 1970) using 
#              Time::HiRes::gettimeofday().
#
# Arguments:   NONE
# 
# Returns:     Number of seconds and microseconds
#              since the epoch.
#
################################################################# 
sub seconds_since_epoch { 
  my $nargs_expected = 0;
  my $sub_name = "seconds_since_epoch()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seconds, $microseconds) = gettimeofday();
  return ($seconds + ($microseconds / 1000000.));
}

#################################################################
# Subroutine : validate_executable_hash()
#
# Purpose:     Given a reference to a hash in which the 
#              values are paths to executables, validate
#              that those files are executable.
#
# Arguments: 
#   $execs_HR: REF to hash, keys are short names to executable
#              e.g. "vecscreen", values are full paths to that
#              executable, e.g. "/usr/local/infernal/1.1.1/bin/vecscreen"
# 
# Returns:     void
#
# Dies:        if one or more executables does not exist
#
################################################################# 
sub validate_executable_hash { 
  my $nargs_expected = 1;
  my $sub_name = "validate_executable_hash()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($execs_HR) = (@_);

  my $fail_str = undef;
  foreach my $key (sort keys %{$execs_HR}) { 
    if(! -e $execs_HR->{$key}) { 
      $fail_str .= "\t$execs_HR->{$key} does not exist.\n"; 
    }
    elsif(! -x $execs_HR->{$key}) { 
      $fail_str .= "\t$execs_HR->{$key} exists but is not an executable file.\n"; 
    }
  }
  
  if(defined $fail_str) { 
    die "ERROR in $sub_name(),\n$fail_str"; 
  }

  return;
}
