# Supporting sage code

This directory contains the supporting proof-of-concept implementation of constant time lattice related algorithms such as BKZ-2 and close vectors enumeration. 
References to the lemmas (theorems, etc). as per the paper "Constant time lattice reduction in dimension 4 with application to SQIsign".

 - `clll.sage` contains the following methods:
   - `bitsize_bnd(B,p)` given a basis `B` and a natural number `p>=1` defining the norm computes the bitsize of integers required by the BKZ-2 algorithm (as per Lemma 2).
   - `tournum( B, p )` given a basis `B` and a natural number `p>=1` defining the norm estimates the number of BKZ-2 tours required to obtain the guarantee from the Lemma 6.
   - `test_subroutines()` performs tests of the BKZ-2 algorithm on randomly chosen lattices. 
 - `estims.sage` conducts tests of the BKZ-2 algorithm on randomly chosen lattices. It also tests if the bitsize of integers was chosen appropriately.
 - `fxp_cvp_sampler.sage` a proof-of-concept implementation of constant time close vector enumeration algorithm.

In order to run BKZ-2 tests launch the file `estims.sage` in your sage session, use `sage estims.sage`.

In order to run enumeration tests launch the file `fxp_cvp_sampler.sage` in your sage session: 
```sage fxp_cvp_sampler.sage```.
                                                                                                                                                                                                                          
