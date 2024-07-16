#!/usr/bin/env bash

NLATTICES=1200 #number of lattices per file
NFILES=12 #number of files


# input to create_ideals.sage
P=33506587778976371636384290201141387
NI=12890707586852862
NU=187181532229268607942
LN=$NLATTICES
NL=11
BT=10
LT=6


echo "Run on ${NFILES} files with ${NLATTICES} ideals in each of them with p=${P}  nl=${NL} bt=${BT} lt=${LT}"

DATAFILES=""
BENCHMARKFILES=""
for i in `seq $NFILES`
do
    FILENAME="p335_nl${NL}_${i}"
    echo "Running benchmarks on file ${FILENAME}"
    make bench_stats
    DATAFILES="${DATAFILES} ${FILENAME}"
    BENCHMARKFILES="${BENCHMARKFILES} bench_${FILENAME}"
done