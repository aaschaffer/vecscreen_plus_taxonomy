# next line (set -e) makes it so that this script exits with non-zero status if any of the commands within it
# (e.g. the compare_vector_matches_wtaxa.pl command) fail with non-zero status.
set -e 
$VECPLUSDIR/scripts/compare_vector_matches_wtaxa.pl --input_summary $VECPLUSDIR/test-files/sample_input_final_step.txt  --input_taxa $VECPLUSDIR/info-files/taxonomy_tree_wlevels.txt --input_artificial_vectors $VECPLUSDIR/info-files/artificial_whole_sequences.txt --input_artificial_segments $VECPLUSDIR/info-files/artificial_intervals.txt --input_univec_sources $VECPLUSDIR/info-files/biological_exclusions.txt  --input_amr $VECPLUSDIR/info-files/UniVec10_vs_amr_distinct_intervals.txt --input_sequences $VECPLUSDIR/test-files/sample_candidates.fa --outfile my_output_final_step.txt
echo 'comparing expected output test-files/sample_output_final_step.txt to observed output my_outputfinal_step.txt'
diff -w $VECPLUSDIR/test-files/sample_output_final_step.txt my_output_final_step.txt
echo 'SUCCESS: Files are identical'
