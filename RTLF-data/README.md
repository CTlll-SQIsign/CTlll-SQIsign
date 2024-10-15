# RTLF tests and data

This folder contains the in- and outputs of the RTLF tests, separated in three main folders:

## datafiles

The (compressed) `datafiles` folder contains data files with 172.800 lattices generated via `create_inputs.sh` run in 12 threads in 12 different folders, each containing a copy of this repo. This creates 12 data files with 1200 lattices each per thread.

### Creating datafiles

`create_inputs.sh` must be run from the root folder, and the datafiles will be created there too.

## measure-results

We ran both SQIsignLLL (the LLL implementation from the round 1 NIST submission of SQIsign, but using our constant-time integers) and our BKZ2 on the data from the `datafiles` folder. For each we divided the data in 3 batches and run 4 threads in parallel on a machine with 6 (physical) cores, to avoid noise due to hyperthreading. The results are reported in the (compressed) folder `measure-results`.

### Running measures

To execute similar tests using BKZ2, coping `set_headers.py` and `runner.sh` to the root directory and run `runner.sh` in the root directory. It takes 3 arguments, which are the folder containing the datafiles in the same structure as the folder in `datafiles.zip`, the identifiant of the thread (determining from which subfolder of datafiles data will be read) and a folder where to write the results. An example call would be `bash runner.sh datafiles 0 results`, assuming folder datafiles is the dezipped version of datafiles and results an existing but empty folder.

It is not possible to run several threads using the same `src` and `build` folders (since they would overwrite each others header files), therefore the all python and c files in the root folder, the  Makfile and the complete `src` and `build` folder must be copied to a thread-specific folder if one wnats to run several threads at same time.

For running measues for SQIsignLLL instead of BKZ2, use `runner_lll.sh` in the root directory instead of `runner.sh`, and otherwise do everything as in the paragraph above. `runner.sh` and `runner_lll.sh` shouldn't be run simultaneously using the same `src` folder.  They can use the same datafiles, but should not use the same results folder, as they will overwrite each others' results.

## analysis

To analyze the data, we use [RTLF](https://github.com/tls-attacker/RTLF). Therefore, we need to create .csv files where our runtimes are separated into two sets depending on some property of the input data (ou lattices). We ran several such tests.

In the analysis, the 1st and the 12th file are ignored, since the 1st served for warming up the processor, and the last one was potentially less impacted by parallelism as other threads might alredy stopped while it was finished. For normtest and orthogonalitytest, the remaining files are merged first by thread, then all together using the `merger.py` script (with slight adaptations for the merge of all threads).

Once a .csv file suitable for RTLF is obtained, the plotting functions in the `plot_runtimes.py` script can be useful for a quick first look on their distribution.

### Normtest

Once benchmarks were also made, the `normtest/csv_maker.py` script can be used to separate the test data sets into two groups x and y depending on the size of the first coefficient of the basis matrix, which is related to the norm for left O_0 ideals in hermite normal form. RTLF can be run on the resulting csv files. For ease of the analysis, the .csv files obtained from the measurements in `measure_results` is provided for both our BKZ and the original LLL (using out integers).


### Orthogonality test

Once benchmarks were also made, the `orthogonalitytest/orthogonality_csv_maker.py` script can be used to separate the test data sets into two groups x and y depending on the orthogonality defect of the input lattices. RTLF can be run on the resulting csv files. For ease of the analysis, the csv files obtained from the measurements in `measure_results` is provided for both our BKZ and the original LLL (using out integers).

For the script to work, one should unpack:
* the `dumpred.zip` archive into `orthogonalitytest`;
* the content of `./RTLF-data/datafiles.zip` into the same directory;
* the content of `./RTLF-data/measure_results.zip` into the same directory;

The archive `dumpred.zip` contains pre-reduced lattices used to generate the csv files. The archive `datafiles.zip` contains not reduced lattice bases. The measures of both LLL and BKZ algorithm are compressed into `measure_results.zip`.

To recompute the reduced lattices and prepare the CSV files, call the script in the following manner: `sage orthogonality_csv_maker.sage --reduce`. If one wants to obtain the CSV files for LLL algorithm instead of BKZ, the `--lll` flag shall be specified: `sage orthogonality_csv_maker.sage --lll`.

All other `.py` scripts are for internal use and shall not be modified to avoid errors.

### Shuffle test

For the shuffletest, it is necessary to run the measure twice: Once on the input data, and once on the input data where the relevant lattices (so those in files 2 to 11) are randomly shuffled. We therefore ran this measure only on 144000 lattices (so 120000 actually analyzed ones) and on 10 threads in parallel (therefore incurring some hyperthreading, so that the measures are less clean and appear "striped" when plotted).

For the shuffling, the file `shuffletest/shuffledata.py` can be used. The CSV file which classified the 2nd measured runtimes by the runtimes of the 1st measure on the same lattice can be obtained using the `shuffletest/make_csv.py` script.

We provide the obtained `shuffletest/shuffletest.csv`, along with a (compressed) folder `shuffletest/shuffletest_data` containing the measurements and the exact shuffling of the data we used. The unshuffled data corresponds to the first 10 subfolders of `datafiles`.
