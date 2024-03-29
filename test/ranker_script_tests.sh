#!/usr/bin/env bash

source ~soft_bio_267/initializes/init_ruby
source ~soft_bio_267/initializes/init_netanalyzer
export PATH=/mnt/home/users/bio_267_uma/cvargasf/dev_gems/NetAnalyzer/bin:$PATH
test_data=ranker_validation/
out=output_test_scripts 
data_to_test=data_test_scripts
mkdir $out


# ranker by seeds -----------------------------------------------------------------------------------------------------------------------------------------------
##string
ranker_gene.rb -k $test_data/kernel_for_validating -n $test_data/kernel_for_validating.lst -s 'A,B' -o $out/ranker_by_seed_string_results
## file
ranker_gene.rb -k $test_data/kernel_for_validating -n $test_data/kernel_for_validating.lst -s $test_data/seed_genes_for_validating -o $out/ranker_by_seed_file_results

# ranker leaving one out ----------------------------------------------------------------------------------------------------------------------------------------
ranker_gene.rb -k $test_data/kernel_for_validating -n $test_data/kernel_for_validating.lst -s $test_data/seed_genes_for_validating -l true -o $out/ranker_leave_one_out_by_seed_results

# ranker with filter --------------------------------------------------------------------------------------------------------------------------------------------
ranker_gene.rb -k $test_data/kernel_for_validating -n $test_data/kernel_for_validating.lst -s $test_data/seed_genes_for_validating -f $test_data/genes2filter_for_validating -o $out/ranker_filter_results

# ranker with top=2 ---------------------------------------------------------------------------------------------------------------------------------------------
ranker_gene.rb -k $test_data/kernel_for_validating -n $test_data/kernel_for_validating.lst -s $test_data/seed_genes_for_validating -t '2' --output_top $out/ranker_top_results


for file_to_test in `ls $out`; do
	echo $file_to_test
	diff $out/$file_to_test $data_to_test/$file_to_test"_to_test"
done

