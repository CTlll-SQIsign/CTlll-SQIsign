# Inputs for the reported benchmarks

The benchmarks in the paper were run for the inputs in this folder. The lattices in data files for the same level are identical, just represented using larger or smaller integers. The remaining files are supposed to be copied into `src/` and `src/lll` in order to set the constants as needed for the datafiles. Be careful to adjust the scratch space in `limbnum.h` to what your machine and GMP version need. 

`bkz_constants.h` can be adjusted manually to the tours one wants to benchmark.

For benchmarking LLL (only  possible on LVL5_72, LVL3_55 and LVL1_37), inputs.h must be edited by adding the lines `#define LLL_BENCH 1` and `#define BKZ_BENCH 0`. 

Once these preparations done, one can directly run `make bench` from the root directory.

The `sage_command` fdil in each folder explains how to generate a data file with the same parameters using `create_ideals.sage`

When running the benchmarks, the output should be saved to a file (except the maybe header lines). An average over each of these files can then be computed using `compute_average.py` with the corresponding file as argument.