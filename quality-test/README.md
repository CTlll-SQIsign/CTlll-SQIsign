# Inputs for the quality test section

The quality tests are fully deterministic, besides the lattice generation in `quality-tests.sh`. 

Here we provide the input files we obtained by running this script. All outputs we mention in the paper can be reproduced by:

- Commenting out the first line using sage in `quality-tests.sh`
- Copying the content of this folder into the root directory.
- Manually setting the variables in `lll_constants.c`, `limbnum.h` and `inputs.h` according to the security level.
- Adjusting the level in `quality-tests.sh` to what you want to run
- Coping `quality-tests.sh` and `bkz_constants.py` to the root directory and running `quality-tests.sh` there

For help with analyzing the results, the script `process_quality_tests.py` is provided.