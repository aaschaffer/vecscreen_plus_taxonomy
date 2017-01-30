# test-scripts.pl: Testing script that runs example commands and compares output files 
#                  with expected files to make sure they are identical.
#
# EPN, Mon Jan 30 09:39:20 2017
# 
use strict;

# make *STDOUT file handle 'hot' so it automatically flushes whenever we print to it
select *STDOUT;
$| = 1;

##############################################################
# Step 1: Ensure we have all input files and expected output
#         files (to compare to actual output files we create)
#         that we need.
##############################################################
printf("Checking that required input files exist ... ");
my $input_dir             = "./test-files";
my $input_root            = "test";
my $input_dir_and_root    = $input_dir . "/" . $input_root;
my @input_suffixes_A      = ("output_internal.txt", "output_terminal.txt");

my $expected_dir          = "./test-files";
my $expected_root         = "expected";
my $expected_dir_and_root = $expected_dir . "/" . $expected_root;
my @expected_suffixes_A   = ("output_combined.txt");

my $output_dir          = "./test-files";
my $output_root         = "tmp";
my $output_dir_and_root = $output_dir . "/" . $output_root;
my @output_suffixes_A   = ("output_combined.txt");

my $errmsg = "";

# check input files exist
foreach my $suffix (@input_suffixes_A) { 
  $errmsg .= check_file_exists_and_is_nonempty($input_dir . "/" . $input_root . "." . $suffix, "input");
}
# check expected output files exist
foreach my $suffix (@expected_suffixes_A) { 
  $errmsg .= check_file_exists_and_is_nonempty($expected_dir . "/" . $expected_root . "." . $suffix, "expected output");
}
# die if any were missing
if($errmsg ne "") { 
  die "\n" . $errmsg;
}
printf("done\n");

#####################################################################
# Step 2: Run commands and check that output files are as expected. #
#####################################################################
printf("Testing combine_summaries.pl ... \n");

######################
my $cmd = "./combine_summaries.pl --input_terminal $input_dir_and_root.output_internal.txt --input_internal $input_dir_and_root.output_terminal.txt --outfile $output_dir_and_root.output_combined.txt > /dev/null"; # standard output is expected to be empty
run_command($cmd, 0);
my @out_A = ("$output_dir_and_root.output_combined.txt");
check_many_files_exist_and_are_nonempty(\@out_A, "output");
diff_files("$output_dir_and_root.output_combined.txt", "$expected_dir_and_root.output_combined.txt");
#####################
# same command, but with --verbose added
my $cmd = "./combine_summaries.pl --input_terminal $input_dir_and_root.output_internal.txt --input_internal $input_dir_and_root.output_terminal.txt --outfile $output_dir_and_root.output_combined.verbose.txt --verbose > /dev/null"; # standard output is expected to be empty
run_command($cmd, 0);
my @out_A = ("$output_dir_and_root.output_combined.verbose.txt");
check_many_files_exist_and_are_nonempty(\@out_A, "output");
diff_files("$output_dir_and_root.output_combined.verbose.txt", "$expected_dir_and_root.output_combined.verbose.txt");
#####################

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
