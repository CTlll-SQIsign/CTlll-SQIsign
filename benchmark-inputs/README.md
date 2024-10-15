# Inputs for the reported benchmarks

The benchmarks in the paper were run for the inputs in this folder. The lattices in data files for the same level are identical, just represented using larger or smaller integers. The remaining files are supposed to be copied into `src/` and `src/lll` in order to set the constants as needed for the datafiles. Be careful to adjust the scratch space in `limbnum.h` to what your machine and GMP version need. 

`bkz_constants.h` can be adjusted manually to the tours one wants to benchmark.

For benchmarking LLL (only  possible on LVL5_72, LVL3_55 and LVL1_37), inputs.h must be edited by adding the lines `#define LLL_BENCH 1` and `#define BKZ_BENCH 0`. 

Once these preparations done, one can directly run `make bench` from the root directory.

The `sage_command` in each folder explains how to generate a data file with the same parameters using `create_ideals.sage`

When running the benchmarks, the output should be saved to a file (except the maybe header lines). An average over each of these files can then be computed using `compute_average.py` with the corresponding file as argument.


### The benchmarks for LLL using non-constant-time integers

For comparison with the current LLL implementation in SQIsign with the variable-size and variable-time integers SQIsign currently uses, we also provide support for running benchmarks of LLL on these integers. Set the header files as described above for the desired level (the integer size does not matter), then run `make bench_nct_lll` from the root directory. Averages can be computed as above.

## The results folder

The results folder contains (examples of) results obtained by these methods.