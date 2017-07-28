#!/usr/bin/perl -w
# the first line of perl code has to be above
#
# test_vecscreen_plus_taxonomy_scriptsscripts.pl: 
# Testing script that runs example commands and compares output files 
# with expected files to make sure they are identical.
#
# EPN, Mon Jan 30 09:39:20 2017
# 
use strict;

# make *STDOUT file handle 'hot' so it automatically flushes whenever we print to it
select *STDOUT;
$| = 1;

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

# make sure the executable files we need exist in $VECPLUSDIR
my %execs_H = (); # hash with paths to all required executables
$execs_H{"vecscreen"}                 = $vecplus_exec_dir . "vecscreen";
$execs_H{"srcchk"}                    = $vecplus_exec_dir . "srcchk";
$execs_H{"parse_vecscreen"}           = $vecplus_exec_dir . "parse_vecscreen.pl";
$execs_H{"combine_summaries"}         = $vecplus_exec_dir . "combine_summaries.pl";
$execs_H{"add_taxonomy"}              = $vecplus_exec_dir . "add_taxonomy.pl";
$execs_H{"from_vecscreen_to_summary"} = $vecplus_exec_dir . "from_vecscreen_to_summary.pl";
$execs_H{"compare_vector_matches"}    = $vecplus_exec_dir . "compare_vector_matches_wtaxa.pl";
validate_executable_hash(\%execs_H);

########################################################################
# Ensure we have all input files and expected output files (to compare
# to actual output files we create) that we need.
########################################################################
# This needs to be manually updated. The code does not ensure that all 
# input files used in the 'Test *.pl' sections below are checked for
# in this block. So, if you add an input or expected file to a 'Test *.pl'
# block, add it here too.
# 

my $command_width = 40;
printf("%-*s ... ", $command_width, "Checking that required input files exist");
my $input_dir             = "$vecplusdir/test-files";
my $input_root            = "test";
my $input_dir_and_root    = $input_dir . "/" . $input_root;
my @input_suffixes_A      = ("input_sequence_file.fa", "sample_input_final_step.txt", "sample_candidates.fa");
# and a second input dir, with the big taxonomy_tree file and the exclusions files required by compare_vector_matches_wtaxa.pl
my $input_dir2            = "$vecplusdir/info-files";
my $input_dir_and_root2   = $input_dir2 . "/";
my @input_suffixes2_A     = ("taxonomy_tree_wlevels.txt", # this file is too big to require a copy to exist in test-files
                             "artificial_whole_sequences.txt",
                             "artificial_intervals.txt",
                             "biological_exclusions.txt",
                             "superkingdom_level_exclusions.txt",
                             "kingdom_level_exclusions.txt",
                             "phylum_level_exclusions.txt",
                             "class_level_exclusions.txt",
                             "order_level_exclusions.txt",
                             "family_level_exclusions.txt",
                             "tribe_level_exclusions.txt", 
                             "genus_level_exclusions.txt", 
                             "UniVec10_vs_amr_distinct_intervals.txt");

my $expected_dir          = "$vecplusdir/test-files";
my $expected_root         = "expected";
my $expected_dir_and_root = $expected_dir . "/" . $expected_root;
my @expected_suffixes_A   = ("output_combined.txt", "output_combined.verbose.txt", "add_taxonomy.out", "add_taxonomy.verbose.out", "output_final_step.txt");

my $output_dir          = "$vecplusdir/test-files";
my $output_root         = "tmp";
my $output_dir_and_root = $output_dir . "/" . $output_root;

my $output_dir2          = "./";
my $output_root2         = "tmp";
my $output_dir_and_root2 = $output_dir2 . $output_root2;

my $errmsg = "";

# check input files exist
foreach my $suffix (@input_suffixes_A) { 
  $errmsg .= check_file_exists_and_is_nonempty($input_dir . "/" . $input_root . "." . $suffix, "input");
}
foreach my $suffix (@input_suffixes2_A) { 
  $errmsg .= check_file_exists_and_is_nonempty($input_dir2 . "/" . $suffix, "input");
}
# check expected output files exist
foreach my $suffix (@expected_suffixes_A) { 
  $errmsg .= check_file_exists_and_is_nonempty($expected_dir . "/" . $expected_root . "." . $suffix, "expected output");
}
# die if any were missing
if($errmsg ne "") { 
  die "\n" . $errmsg;
}
printf("done.\n");

###########################
# Test combine_summaries.pl
###########################
printf("%-*s ... ", $command_width, "Testing combine_summaries.pl");
######################
my $cmd = $execs_H{"combine_summaries"} . " --input_terminal $input_dir_and_root.output_internal.txt --input_internal $input_dir_and_root.output_terminal.txt --outfile $output_dir_and_root.output_combined.txt > /dev/null"; # standard output is expected to be empty
run_command($cmd, 0);
my @out_A = ("$output_dir_and_root.output_combined.txt");
check_many_files_exist_and_are_nonempty(\@out_A, "output");
diff_files("$output_dir_and_root.output_combined.txt", "$expected_dir_and_root.output_combined.txt");
rm_files(\@out_A);

#####################
# same command, but with --verbose added
$cmd = $execs_H{"combine_summaries"} . " --input_terminal $input_dir_and_root.output_internal.txt --input_internal $input_dir_and_root.output_terminal.txt --outfile $output_dir_and_root.output_combined.verbose.txt --verbose > /dev/null"; # standard output is expected to be empty
run_command($cmd, 0);
@out_A = ("$output_dir_and_root.output_combined.verbose.txt");
check_many_files_exist_and_are_nonempty(\@out_A, "output");
diff_files("$output_dir_and_root.output_combined.verbose.txt", "$expected_dir_and_root.output_combined.verbose.txt");
rm_files(\@out_A);

#####################
printf("done.\n");

######################
# Test add_taxonomy.pl
######################
printf("%-*s ... ", $command_width, "Testing add_taxonomy.pl");
######################
$cmd = $execs_H{"add_taxonomy"} . " --input_summary $input_dir_and_root.combine_summaries --input_taxa " . $input_dir_and_root2 . "taxonomy_tree_wlevels.txt --outfile $output_dir_and_root.add_taxonomy.out"; 
run_command($cmd, 0);
@out_A = ("$output_dir_and_root.add_taxonomy.out");
check_many_files_exist_and_are_nonempty(\@out_A, "output");
diff_files("$output_dir_and_root.add_taxonomy.out", "$expected_dir_and_root.add_taxonomy.out");
rm_files(\@out_A);

#####################
# same command, but with --verbose added
$cmd = $execs_H{"add_taxonomy"} . " --input_summary $input_dir_and_root.combine_summaries --input_taxa " . $input_dir_and_root2 . "taxonomy_tree_wlevels.txt --outfile $output_dir_and_root.add_taxonomy.verbose.out --verbose"; 
run_command($cmd, 0);
@out_A = ("$output_dir_and_root.add_taxonomy.verbose.out");
check_many_files_exist_and_are_nonempty(\@out_A, "output");
diff_files("$output_dir_and_root.add_taxonomy.verbose.out", "$expected_dir_and_root.add_taxonomy.verbose.out");
rm_files(\@out_A);

#####################
printf("done.\n");

###################################
# Test from_vecscreen_to_summary.pl
###################################
printf("%-*s ... ", $command_width, "Testing from_vecscreen_to_summary.pl");
######################
$cmd = $execs_H{"from_vecscreen_to_summary"} . " --keep --output_root $output_dir_and_root2 --input_fasta $input_dir_and_root.input_sequence_file.fa --input_taxa " . $input_dir_and_root2 . "taxonomy_tree_wlevels.txt --verbose --combine_output > /dev/null"; # standard output is expected to be empty
run_command($cmd, 0);
@out_A = ("$output_dir_and_root2.output_combined_wtaxonomy.txt");
check_many_files_exist_and_are_nonempty(\@out_A, "output");
diff_files("$output_dir_and_root2.output_combined_wtaxonomy.txt", "$expected_dir_and_root.output_combined_wtaxonomy.txt");
rm_files(\@out_A);

######################################
# Test compare_vector_matches_wtaxa.pl
######################################
printf("%-*s ... ", $command_width, "Testing compare_vector_matches_wtaxa.pl");
######################
$cmd = $execs_H{"compare_vector_matches"} . " --input_summary $vecplusdir/test-files/test.sample_input_final_step.txt  --input_taxa $vecplusdir/info-files/taxonomy_tree_wlevels.txt --input_artificial_vectors $vecplusdir/info-files/artificial_whole_sequences.txt --input_artificial_segments $vecplusdir/info-files/artificial_intervals.txt --input_univec_sources $vecplusdir/info-files/biological_exclusions.txt  --input_amr $vecplusdir/info-files/UniVec10_vs_amr_distinct_intervals.txt --input_sequences $vecplusdir/test-files/test.sample_candidates.fa --outfile $output_dir_and_root2.output_final_step.txt"; 
run_command($cmd, 0);
@out_A = ("$output_dir_and_root2.output_final_step.txt");
check_many_files_exist_and_are_nonempty(\@out_A, "output");
diff_files("$output_dir_and_root2.output_final_step.txt", "$expected_dir_and_root.output_final_step.txt");
rm_files(\@out_A);

#####################
printf("done.\n");

printf("# All tests passed.\n");
printf("# SUCCESS\n");

#################################################################
# Subroutine:  check_many_files_exists_and_are_nonempty()
# Incept:      EPN, Mon Jan 30 10:26:28 2017
#
# Purpose:     Checks if several expected output files exist and 
#              are non-empty.
#
# Arguments:
#   $file_AR:  ref to array of files to check for
#   $type:     type of file (e.g. 'input' or 'output')
#
# Returns:    nothing
#
# Dies:       If >= 1 of the files does not exist. 
#################################################################
sub check_many_files_exist_and_are_nonempty { 
  my $sub_name = "check_many_files_exists_and_are_nonempty()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($file_AR, $type) = @_;

  foreach my $file (@{$file_AR}) { 
    $errmsg .= check_file_exists_and_is_nonempty($file, $type);
  }
  if($errmsg ne "") { 
    die "\n" . $errmsg;
  }

  return;
}

#################################################################
# Subroutine:  check_file_exists_and_is_nonempty()
# Incept:      EPN, Mon Jan 30 09:15:06 2017
#
# Purpose:     Checks if an expected output file exists and 
#              is non-empty.
#
# Arguments:
#   $file:     file to check
#   $type:     type of file (e.g. 'input' or 'output')
#
# Returns:    error message if file does not exist or is empty,
#             if "", then file exists and is non-empty
#
# Dies:       Never
#################################################################
sub check_file_exists_and_is_nonempty { 
  my $sub_name = "check_file_exists_and_is_nonempty()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($file, $type) = @_;

  if   (! -e $file) { return "ERROR, $type file $file does not exist\n"; }
  elsif(! -s $file) { return "ERROR, $type file $file exists but is empty\n"; }

  return "";
}

#################################################################
# Subroutine:  diff_files()
# Incept:      EPN, Mon Jan 30 09:03:43 2017
#
# Purpose:     Diffs two files and reports on any differences.
#
# Arguments:
#   $file1:    first  file to compare
#   $file2:    second file to compare
#
# Returns:    output from diff, if "", then both files are equal
#
# Dies:       if either $file1 or $file2 does not exist
#################################################################
sub diff_files { 
  my $sub_name = "diff_files()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($file1, $file2) = @_;

  if(! -e $file1) { die "ERROR, file $file1 in $sub_name does not exist"; }
  if(! -e $file2) { die "ERROR, file $file2 in $sub_name does not exist"; }

#  printf("diff'ing $file1 and $file2\n");

  return `diff $file1 $file2`;
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
#   $echo_cmd:    '1' to echo command to screen, '0' not to
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if $cmd fails
#################################################################
sub run_command {
  my $sub_name = "run_command()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $echo_cmd) = @_;

  if($echo_cmd) { 
    printf("\tExecuting command: $cmd\n");
  }
  system($cmd);

  if($? != 0) { 
    die "ERROR in $sub_name, the following command failed:\n$cmd\n";
  }

  return;
}

#################################################################
# Subroutine:  rm_files
# Incept:      EPN, Mon Jan 30 11:10:58 2017
#
# Purpose:     Remove the files listed in @{$file_AR}
#              using unlink.
#
# Arguments:
#   $file_AR:  ref to array of files to check for
#
# Returns:    nothing
#
# Dies:       Never.
#################################################################
sub rm_files { 
  my $sub_name = "rm_files()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($file_AR) = @_;

  foreach my $file (@{$file_AR}) { 
    unlink $file; 
  }

  return;
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
