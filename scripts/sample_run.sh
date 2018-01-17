# next line (set -e) makes it so that this script exits with non-zero status if any of the commands within it
# (e.g. the compare_vector_matches_wtaxa.pl command) fail with non-zero status.
set -e 
$VECPLUSDIR/scripts/from_vecscreen_to_summary.pl --output_root mytest --input_fasta $VECPLUSDIR/test-files/test.input_sequence_file.fa --input_taxa $VECPLUSDIR/info-files/taxonomy_tree_wlevels.txt --combine_output --verbose > mytest.out
echo 'comparing expected output test-files/expected.output_combined_wtaxonomy.txt to observed output mytest.output_combined_wtaxonomy.txt'
diff -w $VECPLUSDIR/test-files/expected.output_combined_wtaxonomy.txt mytest.output_combined_wtaxonomy.txt
echo 'SUCCESS: Files are identical'
