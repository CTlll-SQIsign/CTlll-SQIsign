# CTLLL-SQIsign

Constant-time lattice reduction for SQIsign

We will first explain how to use the present code, then we will detail the structure of the codebase. 

------

##Origins of the code


Many parts of this code are adapted from the NIST submission of [SQIsign](https://sqisign.org/). 

In particular the file `src/lll/lll.c` is directly taken from that implementation. Furthermore, almost all quaternion-related functions in `src/ct_helpers/ct_helpers.c`, `src/test_helpers/nct_helpers.c` and `src/test_helpers/test_helpers.c`  are taken over from SQIsign, and the `src/nct_intbig` module is obtained from SQIsign's intbig module by just renaming the functions to add the `nct_` prefix. Similarly, `src/nct_lll/lll.c` is an adaptation of SQIsign's `lll.c` file to use the `nct_ibz` functions.

The code in `src/bkz`, as well as most test code at the root of the `src` folder and all code outside of the `src` folder was newly written for this project.

------

## Prerequisites and first step

### Prerequisites

- Strictly required: A 64-bit little-endian architecture, gcc, the GNU multiprecision library GMP (version at least 6.2.0). 
- Required for generation scripts: Python and a recent version of sage
- Required for ctgrind/valgrind: Not an ARM processor (intel/AMD work)

### First steps

Type `make test`. Wait a few seconds, and check if the tests pass.

-----

## Generate inputs

Most executables described below require the three header files `src/inputs.h`, `bkz_constants.h` and `src/limbnum.h` to be set correctly and consistently. Some of them can also be run on inputs from a data file `data_file` in the root directory instead of `src/inputs.h`.

- `limbnum.h` should declare one constant called `LIMBNUM` containing the length of the integers (in 64-bit words). This file should also define the amount of scratch space (in words) needed by the specific machine and GMP version for performing multiplication (`MUL_SCRATCH`) and division (`DIV_SCRATCH`) of two integers of `LIMBNUM` limbs. This can be computed using GMP's `mpn_sec_mul_itch` and `mpn_sec_div_qr_itch` functions on inputs `LIMBNUM,LIMBNUM`, or by using the `get_scratch_size.c` file. 
- `bkz_constants.h` should declare three constants: `BKZ_TOURS` for the number of outer tours in each call to BKZ, `LAGRANGE_TOURS` for number of loop iteration in each call to the constant-time version of Lagrange's 2-dimensional lattice reduction, and `BENCHMARK_ITERATIONS` for the number of iterations over the list of lattices in `inputs.h` used for benchmarking. The last constant is only used by bechmarking via `bench`, and can be omitted otherwise. The whole file is not required for running `small_test`
- `inputs.h` can be used in two variations. If no data file will be used, it should declare an array of lattices `LATTICE_LIST` and an array of big integers `NORM_LIST` containing their norms (unused but useful for debugging). In this case it should also define `LATTICE_NUMBER` to be the number of lattices in `LATTICE_LIST`.  If a data file will be used, it should define a number `LATTICE_DATA_NUMBER` set to the number of lattices to be read from the file, and `DATA_FILE` which should contain the name of the data file used (usually `data_file`, but it can be changed). In both cases, it must declare and set a static big integer `QUATALG_P` (should be prime congruent to 3 mod 4) defining the quaternion algebra to be used. All integers must be given by their little-endian representation with as many limbs as `LIMBNUM` specifies. If this file defines constants called `LLL_TEST` (or `LLL_BENCH`) to 1 instad of 0, tests (or benchmarks) are also run on an implementation of LLL taken from SQIsign's NIST submission and slightly adapted to work with our integers. Similarly, defining constants called `LLL_TEST` (or `LLL_BENCH`) to 0 can disable tests (or benchmarks) for our BKZ-2 variant. If it defines `NO-DIVISION`, then the division by norm optimization (from the section on benchmarking) is disabled. It is enabled by default.
- `inputs.c` can be just a blank file, if for example a datafile is used. Its purpose is to contain the definitions of the lists `LATTICE_LIST` (and optionally `NORM_LIST`) if these are declared as external in `inputs.h`.
- (optional) `data_file` (or a file with another name, equal to the name in the `DATA_FILE` constant in `inputs.h`) contains the binary representation of at least `LATTICE_DATA_NUMBER` input lattices, in which each integer is represented using `LIMBNUM` 64-bit limbs. 
- (for LLL) `lll/lll_constants.c` must define the `ibz_t` constants for 0,1,2 and 3 with the number of limbs defined in `limbnum.h`.

Given the interdependencies between these files (in particular the dependency of most others on `limbnum.h`), manipulating these files by hand is likely to create non-compiling code. Therefore, some tools are given below.

### Automatic input generation using sage

Sage can be used to generate input ideals and their norms. The file `create_ideals.sage` gives an example of how to do that. So far it generates random ideals of a specific order (called O0 in all papers on SQIsign), with norm in a given range. The file must be exeuted with some inputs. Calling it with `-h` as argument gives a description of these. To generate a data file instead of a `inputs.c` file, run the script with an additional argument `-data`. The name of the datafile can be set by using `datafile=<name>` as last argument, it is `data_file` by default.

The toy example in the datafile `p335_nl11_1_small` given inputs.h was obtained by running `sage create_ideals.sage 33506587778976371636384290201141387 12890707586852862 187181532229268607942 9 11 10 6 -data datafile=p335_nl11_1_small`. 


### Automatic input from SQIsign prints

As explained at the beginning of the `parse_output.py` file, it is possible to slightly modify the SQIsign NIST submission in order to print out all lattices encountered in a run. The name of a file containing this output and an appropriate choice for `LIMBNUM`, `BKZ_TOURS`, `LAGRANGE_TOURS` and `BENCHMARK_ITERATIONS` can then be given as input to `parse_output.py` which will set `input.h`, `limbnum.h` and `bkz_constants.h` accordingly. 

An example of such a trace is given with the file `sqisign_lvl1_lll_calls_trace`, and a good limbsize is for example 100. The default are BKZ and Lagrange tour counts of 30 and 10, these can however be adjusted in the script.

To generate a data file, run the script with an additional argument `-data`. The name of the datafile can be set by using `datafile=<name>` as last argument, it is `data_file` by default.

### Custom input generation

To generate input ideals in another way, the functions in `serializer.py` can be used. Some notes on this: 

- In order to use `combine_lattices`, the lattices must be passed as pairs *(d,M)* where *d* is the common denominator and *M* is a matrix whose columns (divided by *d*) are algebra elements. It outputs the content of the `inputs.h` and the `inputs.c` as a tuple.
- The `standard_file_write` function only writes `input.h` and `inputs.c`, the files `limbnum.h` and `lll_constants.c` must be generated by a call to `adapt_limbnum` with the correct number as argument.
- The `adapt_limbnum` function sets not only `limbnum.h` but also adapts `lll_constants.c`.
- The `set_bkz_constants` function writes the given parameters directly to `bkz_constants.h` and defaults `BENCHMARK_ITERATIONS` to 1 if it is not given. 
- The file `data_serializer.py` provides additional functions for generating data files.
- The file `compute_scratch.py` and its dependency `get_scratch_size.c` compute the required scratch space to be declared in `limbnum.h`.

------

## Build and execute

First, create a directory called build in this folder. 

There are 4 main executables which can be built and run: 

### Small_tests

Runs some fixed and not randomized tests on submodules

Build and run it with `make small_test`. 

This is the only executable which can be run without using `inputs.h`, but it does require setting `LIMBNUM` and the scratch space in `limbnum.h`, and the constants in `lll_constants.c`.

### Test

Runs tests ensuring that BKZ reduces all lattices in the `inputs.h` file. 

The tests verify that the bound from Minkowski's second theorem and the bound on the shortest vector are respected. Furthermore, it is veriied that the output is a basis of the same lattice then the input.

Compile and run it with `make test`

Requires `inputs.h` (with corresponding `inputs.c`or datafile), `bkz_constants.h`, `lll_constants.c` and `limbnum.h`to be generated correctly and consistently. 

If `DATA_FILE` and `LATTICE_DATA_NUMBER` are defined in `inputs.h`, it the script will be compiled so that it uses these. The lattices in the data file are assumed to be integral ideals of the given quaternion algebra. If one of these macros is undefined, the compiler will check if `LATTICE_NUMBER` is defined and use the lattices from `LATTICE_LIST`. To be run on LLL instead or in addition to BKZ, the `BKZ_TEST` and `LLL_TEST` can be used as described above.

### Benchmark

Runs benchmarks on BKZ and LLL, with a fixed number of iterations on each lattice from `input.h`. 

Compile and run it with `make bench`

Requires `inputs.h` (with corresponding `inputs.c`or datafile), `bkz_constants.h` (with `BENCHMARK_ITERATIONS`), `lll_constants.c` and `limbnum.h`to be generated correctly and consistently.

If `DATA_FILE` and `LATTICE_DATA_NUMBER` are defined in `inputs.h`, it the script will be compiled so that it uses these, otherwise it will  check if `LATTICE_NUMBER` is defined and use the lattices from `LATTICE_LIST`. To be run on LLL instead or in addition to BKZ, the `BKZ_BENCH` and `LLL_BENCH` can be used as described above.


### Run

Runs BKZ once on each lattice in `input.h`. This should always run in the same time independently of the lattices in `inputs.h`, as long as the variables `QUATALG_P`, `LATTICE_NUMBER` and `LIMBNUM` are unchanged.

Compile and run it with `make run`. You can also compile it using `make` and run it using `./build/run`.

Requires `inputs.h` (with corresponding `inputs.c`or datafile), `bkz_constants.h` and `limbnum.h` to be generated correctly and consistently. 

If `DATA_FILE` and `LATTICE_DATA_NUMBER` are defined in `inputs.h`, it the script will be compiled so that it uses these. Reading from the file is however not guaranteed to run in file-content-independent time. Otherwise it will check if `LATTICE_NUMBER` is defined and use the lattices from `LATTICE_LIST`.

-----

## Additioal make targets for our measures

For the measures reported in the paper, some additional targets can be build and run.

### Test, benchmark, write to file

This is only available if `DATA_FILE` and `LATTICE_DATA_NUMBER` are defined in `inputs.h`. It runs tests on `NTESTS` lattices in the file, where `NTESTS`can be defined in `inputs.h` and otherwise defaults to 7. After this, benchmarks are run on `LATTICE_DATA_NUMBER` lattices in the datafile, and the results  for each lattice are written to a file named like the datafile prefixed with `bench_`. Our post_processing script `post_processor.py` can be run on this output and the datafile to produce a .csv file for [RLFT](https://github.com/tls-attacker/RTLF). 

Compile and run it with `make bench_stats`.

In case the LLL implementation from SQIsign's NIST submission and not our BKZ should be benchmarked, run `make bench_stats_lll` instead.

### Quality test: Test with verbous output

`make quality_test` compiles and runs the tests for Minkowski reduction and for the length of the shortest vector. It does not only report the result, but also prints both sides of the inequalities it tests. This is used in out quality tests. 

In order to get this verbous output from noral tests, it is sufficient to compile them with the additional flag `-DPRINT_FLAG=1`.

### ctgrind

`make bench_ctgrind` will compile and run BKZ-2 on the first lattice in the used datafile while using ctgind (that is, valgrind with the input lattice declared as uninitialized memory). Be aware that this requires to have Valgrind installed, which is not possible on ARM Macs.

### bench_nct_lll

`make bench_nct_lll` will compile and run benchmarks for the LLL implementation in `nct_intbig_lll`.

-----

## Code structure and origin by subfolders

### Files directly in `src`

Files directly in the `src`folder are either files corresponding to make targets (`bench_ctgrind.c`, `bench.c`, `bench_stats.c`, `bench_stats_lll.c`, `run.c`, `small_test.c`, `test.c`), common headers used by these (`bench.h`, `read_data.h`) or input and header files which need to be adapted frequently (`inputs.h`, `inputs.c`,`limbnum.h`, `bkz_constants.h`). The files `ctgrind.h` and `ctgrind.c` are from ctgrind and only used for the target `bench_ctgrind`. 

The make targets were all written specifically for this projects, as well as most of the variable headers and `read_data.h`. The ctgrind files are from ctgrind, and `bench.h` is an adaptation from the `bench.h` file in the NIST submission of SQIsign.

### ct_intbig

This folder contains our constant-time integer and rational arithmetic in the file `ct_intbig.c` and publishes its functions in `ct_intbig.h`. The signatureof most functions in this files matches exactly the signature of a function in the *intbig* module of the SQIsign NIST submission. This is intended. We added very few additional functions, mostly for debugging purposes. A few functions on 4x4 matrices are also comntained in these files, even though their equivalent in the NIST submission is in the quaternion module. 

Tests for the integer functions are in `test_ibz.c`, and tests for the rationals in `test_ibq.c`. `test_ct_intbig.h` declares the functions in the test library. They are only run by the target `small_test`.

The only external dependencies of this folder are GMP and `limbnum.h`.

### ct_helpers

This folder declares (in `ct_helpers.h`) and defines (in `ct_helpers.c`) functions which are not exactly integer functions, but are still in constant time and most of them are required by ou BKZ implementation (but not by LLL). Most of these functions already existed in the NIST submission code, in which case we took them over and adapted if they were not constant-time. A few functions were added to fit the requirements of BKZ. This folder only depends on `ct_intbig`.

### bkz

The files `bkz.c`and `bkz.h` contain the implementation of our constant-time BKZ-2 implementation. The code is written specifically for BKZ-2, using the tools from `ct_helpers`and `ct_intbig`. It closely follows the sage implementation, but uses the reduced norm throughout and implements the dicǘision by the gcd of the gram matrix, as long as `NO_DIVISION` is undefined at compile time.

The folder also contains some static tests for BKZ-2 and its subfunction in the files `test_bkz.h`and `test_bkz.c`. They are only run by the target `small_test`.

### lll

The `lll.c`file contains the exact code from the file `src/quaternion/ref/generic/lll.c` of the Round 1 NIST submission of SQIsign. The only difference is the include in the beginning. The LLL function it defines is declared in `lll.h`. 

A small static test for sanity-checking is available in `test_lll.c`and exported in `test_lll.h`. It is only run by the target `small_test`.

The `lll_constants.h` and `lll_constants.c` files declare and define the constants 0,1,2 and 3 which LLL requires even though our `ct_intbig` library does not implement them.

### nct_intbig

This module, with its files `nct_intbig.c` and `nct_intbig.h`,  is an exact copy of the *intbig* module of the SQIsign NIST implementation, except that all function and type names are prefixed by `nct` which stands for non-constant time. This variable-size and variable-time integer type based on GMP is only used in tests. It allows us to fix the size of our constant-time integers only in function of their use in BKZ or LLL and not on their use in the tests. This is important since some tests compute Hermite Normal Forms and therefore require larnge integers. Furthermore the use of this integer type speeds up the tests. The integers from SQIsign were used because there already exists an implementation of various quaternion algebra operations on top of it, which we can use for our tests. 

### nct_intbig_lll

Contains the original LLL implementation from the SQIsign round 1 NIST submission, modified only so that it uses our names for the nct_intbig functions (which are the sames than those used in the SQIsign NIST submission).

### test_helpers

The files `nct_helpers.h` and `nct_helpers.c` mostly contains functions from the quaternion module of the Round 1 NIST implementation of SQIsign, which are build on top of the `nct_intbig` module. 

The files `test_helpers.c` and `test_helper.h` use the tools from `nct_helpers` and `nct_intbig`. More that half of the file consists in more complex tools for operations on algebra elements and lattices which are taken over and adapted from the NIST implementation of SQIsign. The remaining functions (which are not marked as impoted in the headerfile) are either basic building blocks for these functiosn, or adaptations for specific tests. 

The files `test_code.c` and `test_code.h` provide the tests functions used in the targets `test` and `quality_test`, group them together with BKZ or LLL and make them easy to use. This file relies onm both `ct_intbig` (`lll`, `bkz`. and `ct_helpers`) as well as `nct_intbig` (`nct_helpers` and `test_helpers`), and provides a translation function from constant-time to non-constant time integers. 

The tests in file `test_test_helpers.c` and `test_test_helpers.h` only cover part of the functions in `test_helpers.h`, and most of them are taken from the unit tests of SQIsign*s quaternion module.

The file `printers.h` provides print functions for some `nct`-structures for debugging purposes. 


------

## Data folders

There are three folders with data and scripts from our experiments. Each of them contains a README giving more detail on its content.

### RTLF-data

Contains the data, benchmarks, input files, raw outputs and plots from our experiments with RTLF.

### quality-tests

Contains input data and scripts for our quality test.

### benchmark-inputs

Provides input data for our benchmarks comparing LLL with BKZ.


