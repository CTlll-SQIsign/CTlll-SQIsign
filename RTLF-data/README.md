# RTLF tests and data

This folder contains the in- and outputs of the RTLF tests, separated in three main folders:

## datafiles

This (compressed) folder contains data files with 144.000 lattices generated via `create_inputs.sh` run in 10 threads in parallel in 10 different folders, each containing a copy of this repo. This creates 12 data files with 1200 lattices each per thread. In the analysis, the 1st and the 12th file are ignored, since the 1st served for warming up the processor, and the last one was potentially less impacted by parallelism as other threads might alredy stopped while it was finished. The remaining files were merged first by thread, then all together using the `merger.py` script (with slight adaptations for the merge of all threads). 

`create_inputs.sh` mmust be run from the root folder, and the datafiles will be created there too. 

Once benchmarks were also made, the `post_processor.py` script can be used to separate the test data sets into two groups x and y depending on the size of the first coefficiebt of the basis matrix.

## results_bkz_10_cores

This contains, again ordered by thread and merged, the benchmarks of our run on all datafiles from the datafiles folder, in 10 threads in parallel using our BKZ2 implementation. This folder also includes plots of the output and the RTLF results.

This can be executed by running `runner.sh` in the root directory, in which must be in the state where `create_inputs.sh` left (or would leave) it. This means the constants in the headerfiles must be appropriately set, and the data files (with the exact names `create_inputs.sh` gave them) are also in that folder.


## results_lll_10_cores

This folder is identical to folder `results_bkz_10_cores`, except that the data it contains comes from our measures of LLL (in the implementation from the SQIsign Round 1 NIST submission, but with our integers from `src/ct_intbig` instaed of the original and non-constant-time module) on the same inputs and in the same condition.

The data can be optained by running `runner_lll.sh` in the same conditions described above for `runner.sh`.