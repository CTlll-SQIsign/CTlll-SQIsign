
from binascii import hexlify
from compute_scratch import get_scratch

#File writing

def create_file_headers(content, includes, ctsname):
    start = "#ifndef " + ctsname + "\n#define " + ctsname +"\n\n"
    start = start + includes + "\n"
    end = "\n#endif\n"
    all = start + content + end
    return(all)

def write_file(filename, complete_content):
    f = open(filename,"w")
    f.write(complete_content)
    f.close()


def create_c_file(content):
    cc = "#include \"ct_intbig/ct_intbig.h\"\n#include \"inputs.h\"\n\n" + content
    write_file("src/inputs.c",cc)

def standart_file_write(content):
    content_h, content_c = content
    cc = create_file_headers(content_h, "#include \"ct_intbig/ct_intbig.h\"\n","TEST_INPUTS_H")
    write_file("src/inputs.h",cc)
    create_c_file(content_c)

def adapt_limbnum(limbnum):
    content = "#define LIMBNUM " + str(limbnum) + "\n"
    mul, div = get_scratch(limbnum)
    content = content + "#define MUL_SCRATCH " + str(mul) + "\n"
    content = content + "#define DIV_SCRATCH " + str(div) + "\n"
    cc = create_file_headers(content,"","LIMBNUM_H")
    write_file("src/limbnum.h",cc)
    content_lll = "#include \"lll_constants.h\"\n\n"
    content_lll =  content_lll +"const ibz_t ibz_const_zero = "+convert_number(0,limbnum)+";\n"
    content_lll =  content_lll +"const ibz_t ibz_const_one = "+convert_number(1,limbnum)+";\n"
    content_lll =  content_lll +"const ibz_t ibz_const_two = "+convert_number(2,limbnum)+";\n"
    content_lll =  content_lll +"const ibz_t ibz_const_three = "+convert_number(3,limbnum)+";\n"
    write_file("src/lll/lll_constants.c",content_lll)


def set_bkz_constants(bkz_tours, lagrange_tours, benchmark_iterations=1):
    content = "#define BKZ_TOURS " + str(bkz_tours) + "\n"
    content = content + "#define LAGRANGE_TOURS " + str(lagrange_tours) + "\n"
    content = content + "#define BENCHMARK_ITERATIONS " + str(benchmark_iterations) + "\n"
    cc = create_file_headers(content,"","BKZ_CONSTANTS_H")
    write_file("src/bkz_constants.h",cc)

#Conversion

def convert_number(number,limbnum):
    limbsize = 8
    limbval = 256**limbsize-1
    l = []
    for i in range(0,limbnum):
        n = (number >> 64*i)&(limbval)
        nbytes = n.to_bytes(limbsize,"big")
        l.append("0x"+hexlify(nbytes).decode())
    all = "{"
    for s in l:
        all = all + s + ","
    all = all[:-1]+"}"
    return(all)

def convert_integer_matrix(matrix,limbnum):
    l = [[convert_number(x,limbnum) for x in line] for line in matrix]
    m = ["{"+",".join(line)+"}" for line in l]
    s = "{"+",".join(m)+"}"
    return(s)

def combine_lattices(lattice_list,p,limbnum, normlist = []):
    n = len(lattice_list)
    string_list = []
    for lattice in lattice_list:
        d,m = lattice
        ms = convert_integer_matrix(m,limbnum)
        ds = convert_number(d,limbnum)
        string_list.append("{"+ds+","+ms+"}")
    s = "{"+",".join(string_list)+"}"
    s = "const quat_lattice_t LATTICE_LIST["+ str(n)+"] = "+ s + ";\n"
    ps = "static const ibz_t QUATALG_P = " + convert_number(p, limbnum) + ";\n"
    norms = ""
    if(len(normlist)==len(lattice_list)):
        norms = "const ibz_t NORM_LIST[" + str(n) + "] = {"+",".join([convert_number(x,limbnum) for x in normlist])+"};\n"
    number = "#define LATTICE_NUMBER "+str(n) + "\n"
    sh = number + ps + "extern const quat_lattice_t LATTICE_LIST["+ str(n)+"];\n" 
    if(len(normlist)==len(lattice_list)):
        sh = sh + "extern const ibz_t NORM_LIST[" + str(n) + "];\n"
    sc = s + norms
    return(sh, sc)

def test():
    Llist = [(1,[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]),(2,[[2,0,0,0],[0,2,0,0],[0,1,1,0],[1,0,0,1]])]
    standart_file_write(combine_lattices(Llist,103,10,[1,1]))
    adapt_limbnum(10)

if __name__=="__main__":
    test()