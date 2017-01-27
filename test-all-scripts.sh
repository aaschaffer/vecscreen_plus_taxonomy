echo -n 'Testing combine_summaries.pl ...'


echo -n 'Testing combine_summaries.pl ...'
perl combine_summaries.pl --input_terminal test.output_internal.txt --input_internal test.output_terminal.txt --outfile tmp.output_combined.txt
diff tmp.output_combined.txt test-files/expected.output_combined.txt

echo -n 'Testing add_taxonomy.pl ...'
perl add_taxonomy.pl --keep --input_summary test-files/test.combine_summaries --input_taxa test-files/taxonomy_tree_wlevels.txt --outfile tmp.add_taxonomy.out
diff tmp.add_taxonomy.out.tmp_accession.txt test-files/expected.add_taxonomy.out.tmp_accession.txt
diff tmp.add_taxonomy.out.tmp_srcchk.txt    test-files/expected.add_taxonomy.out.tmp_srcchk.txt
diff tmp.add_taxonomy.out                   test-files/expected.add_taxonomy.out
echo ' done.'



# TODO: test verbose and non-verbose
# MAKE SURE TESTS PASS WITH ORIGINAL CODE
# review parse_vecscreen.pl, and add it here.
# add from_vecscreen_to_summary.pl tests here

