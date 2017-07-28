#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# Authors: Alejandro Schaffer and Eric Nawrocki
# Code to combine summaries of vecscreen matches from terminal and internal specific summaries; 
# if strength is available, then Strong goes above Moderate goes above Weak
# Usage: combine_summaries.pl --input_terminal <input terminal file> \   [REQUIRED]
#                             --input_internal <input internal file> \   [REQUIRED]
#                             --outfile <combined output summary file> \ [REQUIRED]
#                             --verbose                                  [OPTIONAL]
#                             --debug                                    [OPTIONAL]

use strict;
use warnings;
use Getopt::Long;

require "epn-options.pm";

# variables related to command line options
my $input_terminal_file = undef; # input fasta file 
my $input_internal_file = undef; # input fasta file 
my $output_file         = undef; # output fasta file for sequences that do not have the forbidden terms
my $verbose_mode;                # did user specify verbose mode, in which case match strengths are kept track of 
                                 # and extra columns are printed
my $debug_mode;                  # did user specify debug mode, in which case additional diagnostics are printed

# hard-coded information on meaning of columns in the input
my $ACC_COLUMN = 0;                # the column with the identifier in the input files
my $STR_COLUMN = 6;                # the column with the strength of the match in the input files
my $OVERALL_STR_COLUMN = 7;        # the column with the overal highest strength for this match in the input files

# variables used to store info from terminal and internal match input files
my $terminal_FH;                   # file handle for $input_terminal_file
my $internal_FH;                   # file handle for $input_internal_file
my @terminal_lines_A;              # array of forbidden terms
my @internal_lines_A;              # number of forbidden terms
my @terminal_accessions_A;         # accessions for each terminal line
my @internal_accessions_A;         # accessions for each internal line
my @terminal_strengths_A;          # strengths for each terminal line
my @internal_strengths_A;          # strengths for each internal line
my @terminal_overall_strengths_A;  # overall strengths for each match
my @internal_overall_strengths_A;  # strengths for each match
my @terminal_printed_A;            # has each terminal line been printed
my @internal_printed_A;            # has each internal line been printed
my %internal_accession_H;          # hash. key: accession, value first line of that accession in the internal set
my %crossfile_overall_strengths_H; # hash. key: accession, value is strongest strength observed for this accession in either file
my $num_internal_lines;            # number of lines in internal input file
my $num_terminal_lines;            # number of lines in terminal input file

# variables related to output
my $current_strength;      # strength of strongest matches in current round of output ("Strong", "Moderate", or "Weak")
my $terminal_idx;          # index of terminal line to be considered for printing
my $internal_idx;          # index of internal line to be considered for printing
my $terminal_accession;    # terminal accession $terminal_accession_A[$terminal_idx]
my $internal_accession;    # internal accession $internal_accession_A[$internal_idx]
my $repeated_accession;    # accession that exists in both terminal and internal files
my $special_internal_idx;  # index of internal line to be considered for printing when accession matches terminal
my $num_printed;           # number of lines printed so far

# variables related to command line options, see epn-options.pm
my %opt_HH = ();
my @opt_order_A = ();
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
#       option                 type       default  group   requires incompat   preamble-outfile                help-outfile
opt_Add("-h",                  "boolean", 0,           0,    undef, undef,     undef,                          "display this help",                           \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"1"} = "required options";
opt_Add("--input_terminal",    "string",  undef,       1,    undef, undef,     "input terminal summary file",  "REQUIRED: file name <s> of terminal summary", \%opt_HH, \@opt_order_A);
opt_Add("--input_internal",    "string",  undef,       1,    undef, undef,     "input internal summary file",  "REQUIRED: file name <s> of internal summary", \%opt_HH, \@opt_order_A);
opt_Add("--outfile",           "string",  undef,       1,    undef, undef,     "output file",                  "REQUIRED: name <s> of output file to create", \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"2"} = "other options (not required)";
opt_Add("--verbose",           "boolean", 0,           2,    undef, undef,     "be verbose",                   "be verbose in output",                        \%opt_HH, \@opt_order_A);
opt_Add("--debug",             "boolean", 0,           2,    undef, undef,     "debugging mode on",            "turn debugging mode on",                      \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $all_options_recognized =
    &GetOptions('h'                => \$GetOptions_H{"-h"},
                'input_terminal=s' => \$GetOptions_H{"--input_terminal"},
                'input_internal=s' => \$GetOptions_H{"--input_internal"},
                'outfile=s'        => \$GetOptions_H{"--outfile"},
                'verbose'          => \$GetOptions_H{"--verbose"},
                'debug'            => \$GetOptions_H{"--debug"});

my $synopsis = "combine_summaries.pl :: combine summaries of terminal and internal vecscreen matches\n\n";
my $usage    = "Usage: combine_summaries.pl ";

my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.01";
my $releasedate   = "Jan 2017";

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# define file names
$input_terminal_file  = opt_Get("--input_terminal", \%opt_HH);
$input_internal_file  = opt_Get("--input_internal", \%opt_HH);
$output_file          = opt_Get("--outfile", \%opt_HH);
$verbose_mode         = opt_Get("--verbose", \%opt_HH);
$debug_mode           = opt_Get("--debug", \%opt_HH);

# Die if any of:
# - non-existent option is used
# - any of the required options are not used. 
# - -h is used
my $reqopts_errmsg = "";
if(! defined $input_terminal_file) { $reqopts_errmsg .= "ERROR, --input_terminal not used. It is required.\n"; }
if(! defined $input_internal_file) { $reqopts_errmsg .= "ERROR, --input_internal option not used. It is required.\n"; }
if(! defined $output_file)         { $reqopts_errmsg .= "ERROR, --outfile option not used. It is required.\n"; }

if($GetOptions_H{"-h"}) {
  opt_OutputHelp(*STDOUT, $synopsis . $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  exit 0;
}
if(($reqopts_errmsg ne "") || (! $all_options_recognized)) { 
  if   ($reqopts_errmsg ne "") { die $reqopts_errmsg; }
  else                         { die "ERROR, unrecognized option;"; }
}


# open the input and output files
open($terminal_FH, "<", "$input_terminal_file") or die "Cannot open $input_terminal_file for reading\n"; 
open($internal_FH, "<", "$input_internal_file") or die "Cannot open $input_internal_file for reading\n"; 
open(OUTPUT,       ">", "$output_file")         or die "Cannot open $output_file for writing\n"; 

########################
# Parse the input files 
########################
$num_terminal_lines = parseInputFile($terminal_FH, \@terminal_lines_A, \@terminal_accessions_A, \@terminal_strengths_A, 
                                     \@terminal_overall_strengths_A, undef, \%crossfile_overall_strengths_H, $verbose_mode, 
                                     $ACC_COLUMN, $STR_COLUMN, $OVERALL_STR_COLUMN);
$num_internal_lines = parseInputFile($internal_FH, \@internal_lines_A, \@internal_accessions_A, \@internal_strengths_A, 
                                     \@internal_overall_strengths_A, \%internal_accession_H, \%crossfile_overall_strengths_H, $verbose_mode, 
                                     $ACC_COLUMN, $STR_COLUMN, $OVERALL_STR_COLUMN);

##########################################################################
# Create and print output 
#
# We loop through 3 times, once for each overall strength: "Strong",
# "Moderate", "Weak" (if --verbose was NOT used then all strengths
# have been stored as "Strong" since strength will not be part of the
# output)
#
# For each overall strength, we print the accessions in the order they
# existed in the terminal file followed by accessions in the order
# they existed in the internal file. For any accessions that existed
# in both files (identified and stored as $repeated_accession) below,
# all lines for that accession are printed consecutively, in decreasing order of 
# strength.
##########################################################################

# initialize *_printed_A arrays
for($terminal_idx = 0; $terminal_idx < $num_terminal_lines; $terminal_idx++) { 
  $terminal_printed_A[$terminal_idx] = 0;
}
for($internal_idx = 0; $internal_idx < $num_internal_lines; $internal_idx++) { 
  $internal_printed_A[$internal_idx] = 0;
}

$num_printed = 0;
foreach $current_strength ("Strong", "Moderate", "Weak") { 
  $terminal_idx = 0;
  $internal_idx = 0;

  while (($terminal_idx < $num_terminal_lines) && ($internal_idx < $num_internal_lines)) {
     if($debug_mode) { 
       print "#DEBUG: terminal $terminal_idx  $terminal_accessions_A[$terminal_idx]\n";
       print "#DEBUG: strength is $crossfile_overall_strengths_H{$terminal_accessions_A[$terminal_idx]}\n";
     }
    if ($terminal_idx < $num_terminal_lines) {
      $terminal_accession = $terminal_accessions_A[$terminal_idx];
      if ($crossfile_overall_strengths_H{$terminal_accession} eq $current_strength) { 

        # do we also have internal hits for this accession?
        if (defined($internal_accession_H{$terminal_accession})) { 
          $special_internal_idx = $internal_accession_H{$terminal_accessions_A[$terminal_idx]};
          $repeated_accession = $terminal_accession;
          if($debug_mode) { 
            print "#DEBUG: Repeated accession is $terminal_accession $terminal_idx  $num_terminal_lines $special_internal_idx $num_internal_lines\n";
          }

          # $repeated_accession exists in both internal and terminal data structures
          # output all lines from either internal or terminal for $repeated_accession in order of strength
          # while there are >= 1 lines left in BOTH terminal and internal data structures
          while (($terminal_idx         < $num_terminal_lines) && ($terminal_accessions_A[$terminal_idx]         eq $repeated_accession) &&
                 ($special_internal_idx < $num_internal_lines) && ($internal_accessions_A[$special_internal_idx] eq $repeated_accession)){
            if (compareStrengths($terminal_strengths_A[$terminal_idx], $internal_strengths_A[$special_internal_idx])) {
              print OUTPUT "$terminal_lines_A[$terminal_idx]";
              print OUTPUT "\n";
              $terminal_printed_A[$terminal_idx] = 1;
              $num_printed++; 
              $terminal_idx++;
            } 
            else {
              print OUTPUT "$internal_lines_A[$special_internal_idx]";
              print OUTPUT "\n";
              $internal_printed_A[$special_internal_idx] = 1;
              $num_printed++; 
              $special_internal_idx++;
            }
          }

          # no longer >= 1 lines left for $repeated_accession in BOTH terminal and internal data structures, finish off any left in terminal
          while (($terminal_idx < $num_terminal_lines) && ($terminal_accessions_A[$terminal_idx] eq $repeated_accession)) {
            print OUTPUT "$terminal_lines_A[$terminal_idx]";
            print OUTPUT "\n";
            $terminal_printed_A[$terminal_idx] = 1;
            $num_printed++; 
            $terminal_idx++;
          }
          # no longer >= 1 lines for $repeated_accession left in terminal, finish off any left in internal
          while (($special_internal_idx < $num_internal_lines) && ($internal_accessions_A[$special_internal_idx] eq $repeated_accession)){
            print OUTPUT "$internal_lines_A[$special_internal_idx]";
            print OUTPUT "\n";
            $internal_printed_A[$special_internal_idx] = 1;
            $num_printed++; 
            $special_internal_idx++;
          }
        } # end of 'if (defined($internal_accession_H{$terminal_accession}))' 
          # which was entered if $terminal_accession is in both internal and terminal files
        else {
          # $terminal_accession is ONLY in terminal data structure, output it (if we didn't previously)
          if ((! ($terminal_printed_A[$terminal_idx])) && ($crossfile_overall_strengths_H{$terminal_accessions_A[$terminal_idx]} eq $current_strength)) { 
            print OUTPUT "$terminal_lines_A[$terminal_idx]";
            print OUTPUT "\n";
            $terminal_printed_A[$terminal_idx] = 1;
            $num_printed++; 
          }
          $terminal_idx++;
        }
      } # end of 'if ($crossfile_overall_strengths_H{$terminal_accession} eq $current_strength)'
      else {
        $terminal_idx++;
      }
    }
    # done with all accessions that had a terminal match, move onto remaining internal accessions

    while (($terminal_idx == $num_terminal_lines) && ($internal_idx < $num_internal_lines)) {
      if($debug_mode) { 
        print "#DEBUG: Internal $internal_idx\n";
        print "#DEBUG: $internal_accessions_A[$internal_idx]\n";
        print "#DEBUG: Internal strength is $crossfile_overall_strengths_H{$internal_accessions_A[$internal_idx]}\n";
      }
      $internal_accession = $internal_accessions_A[$internal_idx];
      if ((! ($internal_printed_A[$internal_idx])) && ($crossfile_overall_strengths_H{$internal_accession} eq $current_strength)) { 
        print OUTPUT "$internal_lines_A[$internal_idx]";
        print OUTPUT "\n";
        $internal_printed_A[$internal_idx] = 1;
        $num_printed++; 
      }
      $internal_idx++;
    }
  }
}
close (OUTPUT);

#################################################################
# Subroutine: parseInputFile()
# Synopsis:   Parses either the terminal or internal input file.
# 
# Args:       $FH:                             file handle to read from
#             $lines_AR:                       ref to array of lines, filled here
#             $accessions_AR:                  ref to array of accessions, filled here
#             $strengths_AR:                   ref to array of strengths, filled here
#             $overall_strengths_AR:           ref to array of overall strengths, filled here
#             $accessions_line_HR:             ref to hash of accessions, keys is accession
#                                              value is line number, filled here but can be undef 
#             $crossfile_overall_strengths_HR: ref to hash of overall strengths, best from either
#                                              file, updated here
#             $verbose_mode:                   '1' if we're in verbose mode, '0' if not
#             $ACC_COLUMN:                     index of accession column in file
#             $STR_COLUMN:                     index of strength column in file
#             $OVERALL_STR_COLUMN:             index of overall strength column in file
#
# Returns: Number of lines read
#
#################################################################
sub parseInputFile { 
  my $sub_name = "parseInputFile";
  my $nargs_exp = 11;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($FH, $lines_AR, $accessions_AR, $strengths_AR, $overall_strengths_AR, $accessions_line_HR, 
      $crossfile_overall_strengths_HR, $verbose_mode, $ACC_COLUMN, $STR_COLUMN, $OVERALL_STR_COLUMN) = @_;
  
  my $num_lines = 0;
  my $prev_accession = undef;
  while(defined(my $nextline = <$FH>)) {
    # Example terminal line (columns are separated by tabs)
    #XM_010492096.1	266	294	uv|KF718979.1:2206-6729	2404	2432	Moderate	Strong	Suspect[1,35];Suspect[259,265];
    # Example internal line (columns are separated by tabs)
    #XM_715279.1	1408	1436	uv|U29899.1:1847-4463	2140	2112	Moderate	Moderate	Suspect[1437,1467];
    chomp($nextline);
    $lines_AR->[$num_lines] = $nextline;
    my @fields = split /\t/, $nextline;
    my $accession = $fields[$ACC_COLUMN];

    # update array elements, unique for each file
    $accessions_AR->[$num_lines]        = $accession;
    $strengths_AR->[$num_lines]         = ($verbose_mode) ? $fields[$STR_COLUMN] : "Strong";
    $overall_strengths_AR->[$num_lines] = $fields[$OVERALL_STR_COLUMN];

    # update the one hash that is cross-file
    if (!defined($crossfile_overall_strengths_HR->{$accession})) {
      $crossfile_overall_strengths_HR->{$accession} = ($verbose_mode) ? $fields[$OVERALL_STR_COLUMN] : "Strong";
    }
    else {
      if ($verbose_mode) {
        if(compareStrengths($fields[$OVERALL_STR_COLUMN], $crossfile_overall_strengths_HR->{$accession})) {
          $crossfile_overall_strengths_HR->{$accession} = $fields[$OVERALL_STR_COLUMN];
        }
      }
    }

    # update the accessions_line hash which keeps track of the first line of each accession, if it's defined
    if(defined $accessions_line_HR) { 
      if ((! defined $prev_accession) || ($prev_accession ne $accession)) {
        $accessions_line_HR->{$accession} = $num_lines; 
      }
    }

    $num_lines++;
    $prev_accession = $accession;
  }
  close $FH;
  return $num_lines;  
}

#################################################################
# Subroutine: compareStrengths()
# Synopsis:   return 1 if first strength is greater than or 
#             equal to second strength
#
# Args: $local_strength1: first strength
#       $local_strength2: second strength
#
# Returns: 1 if $local_strength >= $local_strength2, else 0
#
# Dies: if either $local_strength1 or $local_strength2 is none of:
#       "Strong", "Moderate", or "Weak" (actual 'die' is in 
#       strengthWordToNumber())
#
#################################################################
sub compareStrengths {
    my $sub_name = "compareStrengths()";
    my $nargs_exp = 2;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_strength1, $local_strength2) = @_;

    return((strengthWordToNumber($local_strength1)) >= 
           (strengthWordToNumber($local_strength2)));
}

#################################################################
# Subroutine: strengthWordToNumber()
# Synopsis:   Converts a strength 'word' to a number as follows
#             "Strong"   => 3
#             "Moderate" => 2
#             "Weak"     => 1
#
# Args: $strength_word: strength word 
#
# Returns: strength number as shown in synopsis.
#
# Dies: if $strength_word is not "Strong", "Moderate" or "Weak"
#
#################################################################
sub strengthWordToNumber { 
    my $sub_name = "strengthWordToNumber";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($strength_word) = @_;

    if ($strength_word eq "Strong") { 
      return 3;
    }
    elsif($strength_word eq "Moderate") { 
      return 2;
    }
    elsif ($strength_word eq "Weak") { 
      return 1;
    }
    elsif ($strength_word eq "None") { 
      return 0;
    }

    die "ERROR in $sub_name, unexpected <strength_word>: $strength_word";

    return 0; # NEVER REACHED
}
