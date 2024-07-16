from serializer import standart_file_write, combine_lattices, adapt_limbnum, set_bkz_constants
from data_serializer import set_header_files, write_lattice_to_file
from sys import argv

#Usage of this file 
#Use the main branch of the NIST submission
#Fix the memory leaks by finalizing correctly parts of the signature
#add in the file "src/quaternion/ref/generic/lideal.c", at the beginning of function "quat_lideal_reduce_basis"
#the line: "ibz_printf("prime: %Zd\n",&(alg->p)); quat_left_ideal_print(lideal);""
#compile a release build of the sqisign nist submission
#in its build/test folder, run "./sqisign_test_scheme_lvl1 > output_filename"
#then run the present file with arguments "output_filename 60"
#this sets the input lattices in a header file to those generated during the signature, and the integer size sufficiently high
#compile and run the present project
#to use the data file instead of a c file for storing the input lattices, add -data to the arguments of this file
#to use another name than "data_file" for the data file, add the argument datafile=<chosen-name>

prime = int("0xea6a4dda9518e5c5d50ccdfbd97e4c49efe85e0e09039c7ffffffffffffffff",16)
limbnum = 40
bkz_tours = 30
lagrange_tours = 10
benchmark_iterations = 1

def read_file(filename):
    f = open(filename,"r")
    l = f.readlines()
    f.close()
    return(l)

def parse_output(lines):
    matrices = []
    norms = []
    denominators = []
    i = -1
    while i < len(lines)-1:
        i = i+1
        clean = lines[i].strip()
        words = clean.split(" ")
        if(words[0]=="basis:"):
            m = []
            m.append([int(x,10) for x in words[1:]])
            m.append([int(x,10) for x in lines[i+1].strip().split(" ")])
            m.append([int(x,10) for x in lines[i+2].strip().split(" ")])
            m.append([int(x,10) for x in lines[i+3].strip().split(" ")])
            matrices.append(m)
            i = i+3
        elif(words[0]=="norm"):
            norms.append(int(words[2],10))
        elif(words[0]=="denominator:"):
            denominators.append(int(words[1],10))
        elif(words[0]=="prime:"):
            global prime
            prime = int(words[1],10)
        else:
            continue
    assert(len(matrices)==len(norms))
    assert(len(matrices)==len(denominators))
    lattices = [(denominators[i],matrices[i]) for i in range(0,len(matrices))]
    return(lattices,norms)

def read_input():
    lens = limbnum
    b_tours = bkz_tours
    l_tours = lagrange_tours
    b_iter = benchmark_iterations
    assert(len(argv)>1)
    filename = argv[1]
    if(len(argv)>2):
        lens = int(argv[2],10)
    if(len(argv)>4):
        b_tours = int(argv[3],10)
        l_tours = int(argv[4],10)
    if(len(argv)>5):
        b_iter = int(argv[5],10)
    print("Set constants to :\n Integer size " + str(64*lens) + " bits\n")
    print(str(b_tours) + " BKZ tours\n")
    print(str(l_tours) + " Lagrange tours at each call to Lagrange\n")
    print(str(b_iter) + " Interations on each lattice or benchmarking\n")
    print("Use input file " + filename + "\n")
    return(filename,lens,b_tours,l_tours,b_iter)

def parse_and_write(filename,lens,b_tours,l_tours,b_iter):
    l = read_file(filename)
    lattices, norms = parse_output(l)
    content = combine_lattices(lattices,prime,lens,norms)
    standart_file_write(content)
    adapt_limbnum(lens)
    set_bkz_constants(b_tours,l_tours,b_iter)
    return

def standard_read_write():
    l,lens,b_tours,l_tours,b_iter = read_input()
    parse_and_write(l,lens,b_tours,l_tours,b_iter)
    return

def parse_and_write_successively(filename,l,lens,b_tours,l_tours,b_iter):
    lattices, _ = parse_output(l)
    set_header_files(filename,lens,len(lattices),prime,b_tours,l_tours,b_iter)
    for lat in lattices:
        write_lattice_to_file(lat,lens,filename)
    return

def bench_read_write(filename):
    f = open(filename,"w")
    f.close()
    infile, lens,b_tours,l_tours,b_iter = read_input()
    l = read_file(infile)
    parse_and_write_successively(filename,l,lens,b_tours,l_tours,b_iter)

if __name__=="__main__":
    if("-data" in argv):
        argv = [a for a in argv if a!= "-data"]
        data_file = "data_file"
        for a in argv:
            if "=" in a and a.split("=")[0]=="datafile":
                argv = [v for v in argv if a!= v]
                data_file = a.split("=")[1]
        bench_read_write(data_file)
    else: 
        standard_read_write()