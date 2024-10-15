from serializer import create_file_headers, adapt_limbnum, set_bkz_constants, write_file, convert_number

LIMBSIZE_BITS = 64
LIMBSIZE_BYTES = 8

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

def set_header_files(filename,limbnum, lattice_number, alg_p, bkz_tours, lagrange_tours,benchmark_iterations=1):
    adapt_limbnum(limbnum)
    set_bkz_constants(bkz_tours,lagrange_tours,benchmark_iterations)
    content = "\n#define LATTICE_DATA_NUMBER {}\nconst ibz_t QUATALG_P = {};\n#define DATA_FILE \"{}\"\n".format(lattice_number,convert_number(alg_p,limbnum),filename)
    cc = create_file_headers(content,"#include \"ct_intbig/ct_intbig.h\"","INPUTS_H")
    write_file("src/inputs.h",cc)
    return

