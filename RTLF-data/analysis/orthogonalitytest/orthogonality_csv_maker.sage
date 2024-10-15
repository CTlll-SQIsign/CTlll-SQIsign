"""
The following files required:
serializer.py
compute_scratch.py
data_serializer.py
parse_output.py

Program assumes that the {BKZ, LLL} tests are extracted to RTLF-data directory.
Also, the dumpred.pkl should be extracted in this directory from dumpred.zip.
"""

from parse_output import *
import os, re
from fpylll import *
import time
import numpy as np
import pickle
from copy import deepcopy
from random import randrange
from math import log, sqrt
import argparse

# from ctlll import GSO_obj, lll_const, bilin
load("../../../sage/clll.sage")
# load("clll.sage")

FPLLL.set_precision(512)
LIMBSIZE_BYTES = 8
p = 33506587778976371636384290201141387
lagrange_max_steps = 15
lll_max_tours = 50

def search_files(directory, pattern):
    # Compile the regex pattern
    regex = re.compile(pattern)

    # List all files in the directory
    files = os.listdir(directory)

    # Filter files matching the pattern
    matching_files = [file for file in files if regex.search(file)]

    return matching_files

def convert_int_to_bytes(i,limbnum):
    b = i.to_bytes(limbnum*LIMBSIZE_BYTES,"little")
    return(b)

def convert_lattice_to_bytes(d,M,limbnum):
    b = convert_int_to_bytes(d,limbnum)
    for i in range(4):
        for j in range(4):
            b = b + convert_int_to_bytes(M[i][j],limbnum)
    return(b)

def write_lattice_to_file(L,limbnum,filename):
    d,M = L
    b = convert_lattice_to_bytes(d,M,limbnum)
    f = open(filename,"ab")
    f.write(b)
    f.close()
    return

def convert_bytes_to_int(bts,limbnum):
    i = int.from_bytes(bts, "little")
    return i

def convert_bytes_to_lattice( bts,limbnum,sub=True ):
    #if sub=True, no d is given
    if sub:
        new_bts = bts
        d = None
    else:
        new_bts = bts[limbnum*LIMBSIZE_BYTES:]
        d = convert_bytes_to_int( bts[:limbnum], limbnum )
    M = [ [0 for i in range(4)] for j in range(4) ]
    for i in range(4):
        for j in range(4):
            curind = 4*i + j
            M[j][i] = convert_bytes_to_int( new_bts[curind*limbnum*LIMBSIZE_BYTES:(curind+1)*limbnum*LIMBSIZE_BYTES],limbnum )
    # print(f"{(curind+1)*limbnum} | {len(bts)/limbnum}")
    return d, M

def parse_datafile(dt, nl):
    """
    The datafile is d, lat0, lat1, ... etc
    """
    bytes_per_int = 8*nl
    bytes_per_lattice = bytes_per_int*17
    # try:
    f = open(dt,"rb")
    lats_bytes = f.read()
    f.close()
    lattice_number = len(lats_bytes)//bytes_per_lattice
    L = []
    for i in range(lattice_number):
        index = i*bytes_per_lattice+bytes_per_int
        bts = lats_bytes[index:index+bytes_per_lattice]
        _, M = convert_bytes_to_lattice( bts,limbnum=nl )
        L.append( M )
        # coeff1_bytes = lats_bytes[index: index+bytes_per_int]
        # coeff1 = int.from_bytes(coeff1_bytes,"little")
        # l.append(coeff1)
    return(L)

def are_equiv(u,v):
    u, v = list(u), list(v)
    for uu in u:
        if uu in v:
            v.pop(v.index(uu))
        elif -uu in v:
            v.pop(v.index(-uu))
        else: return False
    return True

def redlat(M):
    B = matrix(M)
    C, _ = lll_const(B, p=p, logger=None)
    G = GSO_obj( C, p=p )
    G.compute_GSO(0,C.nrows())
    isred = G.is_lll_reduced(0.75,0.55)
    assert isred, f"NOT RED!"
    return C

def orth_def(B):
    Bp = matrix(B)
    n = Bp.nrows()
    D = log( abs( det( Bp ) ) ) + log(p)
    skewness = 0
    for b in Bp:
        skewness += log( bilin(b,b,p=p) )/2
    return exp( (D-skewness)/n )

def parse_benchmarkfile(bn):
    # from post_processor.py
    try:
        f = open(bn,"r")
        bn_lines = f.readlines()
        f.close()
        bn_clean_lines = [b.strip() for b in bn_lines]
        bn_measures = [int(n,10) for n in bn_clean_lines]
        return(bn_measures)
    except Exception as err:
        raise err
        print("Failure at opening or reading file {}\n".format(bn))
        exit(1)

if __name__=="__main__":
    """
    Load all 144000 lattices. Ms is a 3d array with lattices Ms[folder][file][lattice]
    """

    parser = argparse.ArgumentParser("Param_Parser")
    parser.add_argument("--reduce", help="If toggled, reduces lattices from datafiles instead of loading the reduced ones.", action='store_true')
    parser.add_argument("--lll", help="If toggled, prepares CSV for LLL instead of BKZ.", action='store_true')
    args = parser.parse_args()

    Ms = []
    name_offset = len("p335_nl11_")
    totlat = 0
    if args.reduce:
        for batchnum in range(12): #1 done
            directory = f"../../datafiles/{batchnum}/" #path to lattices
            pattern = "p*"
            filenames = search_files(directory, pattern)

            print(filenames)
            Ms.append( [] )

            for filename in filenames:
                if (not "p335_nl11_" in filename) or ("all" in filename):
                    continue
                print(f"filenum: {filename[name_offset:]} | {filename}")
                filenum = int(filename[name_offset:])
                latnum = 0
                Ms[batchnum].append( [] )
                then = time.perf_counter()
                for tmp in parse_datafile(directory+filename, nl=11):
                    M_ = tmp
                    M_ = redlat( M_ )
                    skewness = orth_def(M_)
                    Ms[batchnum][-1].append([M_,skewness])
                    totlat += 1
                    if totlat % 50 == 0:
                        print(f"{totlat} out of {144000} done")
                print(f"file {filename} processed in {time.perf_counter()-then}" )
        # - - - save progress - - -
        with open( f"dumpred.pkl", "wb" ) as file:
            pickle.dump( Ms, file )
    # - - - load progress - - -
    if not args.reduce:
        with open( f"dumpred.pkl", "rb" ) as file:
            Ms = pickle.load(file)
    # - - - Now we process benchmarks - - -
    M = deepcopy( Ms )
    for batchnum in range(12):
        directory = f"../../bkz_measures/results_dell/{batchnum}/" #processes BKZ2
        if args.lll:
            directory = f"../../lll_measures/results_lll/{batchnum}/" #processes LLL
            print(f"Proc LLL")
        else:
            directory = f"../../bkz_measures/results_dell/{batchnum}/" #processes BKZ2
            print(f"Proc BKZ")
        pattern = "(bench_p335_nl11_\\d)"
        filenames = search_files(directory, pattern)

        print(filenames)
        then = time.perf_counter()
        cntr = 0
        for filename in filenames:
            if "all" in filename:
                continue
            print(f"file:{directory+filename}, {batchnum,cntr}")
            rtlf_data = parse_benchmarkfile(directory+filename)
            print(f"len rtlf: {len(rtlf_data)} vs {len(M[batchnum][cntr])}")
            for i in range( len(M[batchnum][cntr]) ):
                _, odl = Ms[batchnum][cntr][i]
                M[batchnum][cntr][i] = [odl, rtlf_data[i]]
            cntr += 1

            print(f"file {filename} processed in {time.perf_counter()-then}" )

    # - - - make the data flat - - -
    l = []
    for i in range( len(M) ):
        for j in range(len(M[0])):
            if j==0 or j == 11:
                continue
            for k in range(len(M[0][0])):
                l.append( M[i][j][k] )
    l = sorted( l, key=lambda x: x[0] ) #move it or comment it if no sorting needed
    # - - - save raw data - - -
    with open(f"data.csv","w") as file:
        file.write( "orth_def, cpu\n" )
        for ll in l:
            file.write( f"{ll[0]}, {ll[1]}\n" )
    # - - - save RTLF input - - -
    med = np.median( [ll[0] for ll in l] )
    cntr = 0
    with open(f"xy_orth_def_data.csv","w") as file:
        file.write( "V1,V2\n" )
        for ll in l:
            if ll[0] <= med:
                file.write( f"X,{ll[1]}\n" )
            else:
                file.write( f"Y,{ll[1]}\n" )
            cntr+=1
