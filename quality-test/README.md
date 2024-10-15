# Inputs for the quality test section

The quality tests are fully deterministic, besides the lattice generation in `quality-tests.sh`. 

Here we provide the input files we obtained by running this script. All outputs we mention in the paper can be reproduced by:

- Commenting out the first line using sage in `quality-tests.sh`
- Copying the content of this folder into the root directory.
- Manually setting the variables in `lll_constants.c`, `limbnum.h` and `inputs.h` according to the security level.
- Instead of manual setting these fles, the provided headers can be used, as long as the scratch space in `limbnum.h` is adjusted to the machine.
- Adjusting the level in `quality-tests.sh` to what you want to run
- Coping `quality-tests.sh` and `bkz_constants.py` to the root directory and running `quality-tests.sh` there

Besides we provide the results of the quality tests in folders `lvl1results`, `lvl3results`, `lvl5results` for security levels 1,3,5. The .csv files in these folders aggregate the data and can be used for statistics and plots. Files `lvl1_ln30_nl20`, `lvl3_ln30_nl28`, `lvl5_ln30_nl37` each contain 30 lattices on which the quality tests were executed. 

For help with analyzing the results, the script `process_quality_tests.py` is provided. This script generates the aforementioned .csv files.
