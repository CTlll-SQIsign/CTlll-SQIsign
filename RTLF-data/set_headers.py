#This file sets the headers as required for the runner.sh
import sys
from data_serializer import set_header_files

if __name__=="__main__":
    filename,limbnum,lattice_number,alg_p, bkz_tours,lagrange_tours , outfile= sys.argv[1], int(sys.argv[2],10), int(sys.argv[3],10), int(sys.argv[4],10), int(sys.argv[5],10), int(sys.argv[6],10),sys.argv[7]
    set_header_files(filename,limbnum, lattice_number, alg_p, bkz_tours, lagrange_tours)
    f=open("src/inputs.h","r")
    l = f.readlines()
    f.close()
    l = [e for e in l if e[:-1].strip()!="#endif"]
    l.append("#define BENCH_STATS_OUTPUTFILE \"" + outfile + "\"\n")
    l.append("\n")
    l.append("#endif\n")
    l.append("\n")
    print(l)
    f=open("src/inputs.h","w")
    for e in l:
        f.write(e)
    f.close()
