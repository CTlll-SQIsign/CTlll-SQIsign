from quaternion_tools import representInteger
from serializer import standart_file_write, combine_lattices, adapt_limbnum, set_bkz_constants
from data_serializer import write_lattice_to_file, set_header_files
from sys import argv

intro = "Expected inputs: p ni nu ln nl bt lt bi\n"
p_text = "p prime defining the quaternion algebra (should be 3 modulo 4)\n"
ni_text = "ni lower bound on the norms of the generated ideals\n"
nu_text = "nu upper bound on the norms of the generated ideals, should be above 2sqrt(p)\n"
ln_text = "ln number of ideals to be generated\n"
nl_text = "nl number of 64-bit words used to store each number, must be large enough (for example 4x than nu needs)\n"
bt_text = "bt Number of BKZ outer tours for each call to BKZ\n"
lt_text = "lt Number of lagrange tours for each call to Lagrange's 2-dimensional lattice algorithm\n"
bi_text = "bi Number of iterations on each lattice for benchmarking (optional, defaults to 1)\n"
data_text = "\nTo generate a data file, insert -data at any position between or after the arguments\n"
df_text = "To change the data file name from the default data_file, add as last argument \"datafile=<new_name>\"\n"
helptext = intro + p_text + ni_text + nu_text + ln_text + nl_text + bt_text + lt_text + bi_text + data_text

format_str = "Use the following values (use -h for help):\np  {}\nni {}\nnu {}\nln {}\nnl {}\nbt {}\nlt {}\nbi {}\n"
format_str_data = "Use the following values (use -h for help):\np  {}\nni {}\nnu {}\nln {}\nnl {}\nbt {}\nlt {}\nbi {}\ndata file name: {}\n"

def setup(p):
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O0 = B.quaternion_order([1,i,(i+j)/2,(1+k)/2])
    return(B,O0)

def prepare_hnf(old):
    lines = 4
    cols = 4
    new = [[0 for i in range(cols)] for j in range(lines)]
    for i in range(cols):
        for j in range(lines):
	        new[j][cols-i-1] = old[i][lines-j-1]
    return(matrix(ZZ,new))

def recover_hnf(old):
    lines = 4
    cols = 4
    new = [[0 for i in range(lines)] for j in range(cols)]
    for i in range(lines):
        for j in range(cols):
	        new[j][i] = old[lines-1-i][cols-j-1]
    return(matrix(ZZ,new))

def hnf(mat):
    #modify to compute a HNF over columns
    mp = prepare_hnf(mat)
    m = mp.echelon_form(algorithm="pari")
    mf = recover_hnf(m)
    return mf

#norm should be odd
def get_ideal(B,O0,norm, llp):
    p = B.ramified_primes()[0]
    cofactor = max(2,(p//norm+1).nbits()*llp)
    x = representInteger(norm*2**(cofactor),B, 2*p.nbits())
    if(x!=False):
        I = O0*x+O0*norm
        return(I)
    return(False)

def get_ideals(p,norm_inf,norm_sub,number):
    ideal_list = []
    B,O0 = setup(p)
    it = 0
    llp = p.bit_length().bit_length()
    while it <number:
        it = it+1
        norm = 0
        while norm%2==0:
            norm = randrange(norm_inf,norm_sub)
        res = get_ideal(B,O0,norm,llp)
        if(res!=False):
            ideal_list.append(res)
            print("generated ideal "+str(it))
        else:
            it = it-1
            print("failed ideal "+str(it))
    return(ideal_list)

def transform_ideal_list(ideal_list):
    norm_list = []
    matrix_list = []
    for I in ideal_list:
        norm_list.append(int(I.norm()))
        M = I.basis_matrix().transpose()
        d = M.denominator()
        Mi = d*M
        Mi = Mi.change_ring(ZZ)
        Mi = hnf(Mi)
        Mints = [[int(mi) for mi in line] for line in Mi]
        matrix_list.append((int(d),Mints))
    #matrix_list.sort(key = lambda x : x[1][0]) #sorting w.r.t. the first matrix entry for stat. test
    return(matrix_list,norm_list)

def test():
    p = 103
    Is = get_ideals(p,10,200,2)
    Ms,Ns = transform_ideal_list(Is)
    s = combine_lattices(Ms,int(p),2,Ns)
    print(s)

def create_list_as_header(p,norm_inf,norm_max,lattice_number,limbnum,bkz_tours,lagrange_tours,bench_iterations):
    Is = get_ideals(p,norm_inf,norm_max,lattice_number)
    Ms,Ns = transform_ideal_list(Is)
    for l in range(len(Ms)):
        d,M = Ms[l]
        Mq = matrix(QQ,M)
        Mqt = ((1/d)*Mq).transpose()
        B = QuaternionAlgebra(p)
        i,j,k = B.gens()
        x = [Mqt[l,0]+Mqt[l,1]*i+Mqt[l,2]*j+Mqt[l,3]*k for l in range(0,4)]
        J = B.ideal(x)
        assert(Is[l]==J)
    s = combine_lattices(Ms,int(p),int(limbnum),Ns)
    standart_file_write(s)
    adapt_limbnum(limbnum)
    set_bkz_constants(bkz_tours,lagrange_tours,bench_iterations)

def create_data_file(p,norm_inf,norm_max,lattice_number,limbnum,bkz_tours,lagrange_tours,bench_iterations, filename):
    B,O0 = setup(p)
    it = 0
    llp = p.bit_length().bit_length()
    f = open(filename,"w")
    f.close()
    set_header_files(filename,limbnum,lattice_number,p,bkz_tours,lagrange_tours,bench_iterations)
    while it <lattice_number:
        it = it+1
        norm = 0
        while norm%2==0:
            norm = randrange(norm_inf,norm_max)
        res = get_ideal(B,O0,norm,llp)
        if(res!=False):
            ms, _ = transform_ideal_list([res])
            m = ms[0]
            write_lattice_to_file(m,limbnum,filename)
            print("generated ideal "+str(it))
        else:
            it = it-1
            print("failed ideal "+str(it))
    print("Done")
    return

def parse_inputs_then_run():
    if(len(argv)>1):
        if(argv[1]=="-h"):
            print(helptext)
            return
    data = False
    v = argv
    if("-data" in argv):
        v = [a for a in v if a!= "-data"]
        data = True
    if(len(v)<=7):
        print("Too few arguments, use -h for help\n")
        return
    bi = 1
    p = int(v[1],10)
    norm_inf = int(v[2],10)
    norm_max = int(v[3],10)
    lattice_number = int(v[4],10)
    limbnum = int(v[5],10)
    bkz_tours = int(v[6],10)
    lagrange_tours = int(v[7],10)
    if(data):
        data_file = "data_file"
        if(len(v)>9):
            bi = int(v[8],10)
            assert(v[9].split("=")[0]=="datafile")
            data_file = v[9].split("=")[1]
        if(len(v)==9):
            if("=" in v[8]):
                assert(v[8].split("=")[0]=="datafile")
                data_file = v[8].split("=")[1]
            else: 
                bi = int(v[8],10)
        create_data_file(p,norm_inf,norm_max,lattice_number,limbnum,bkz_tours,lagrange_tours,bi,data_file)
        print(format_str_data.format(p,norm_inf,norm_max,lattice_number,limbnum,bkz_tours,lagrange_tours,bi,data_file))
    else:
        if(len(v)>8):
            bi = int(v[8],10)
        create_list_as_header(p,norm_inf,norm_max,lattice_number,limbnum,bkz_tours,lagrange_tours,bi)
        print(format_str.format(p,norm_inf,norm_max,lattice_number,limbnum,bkz_tours,lagrange_tours,bi))

def test():
    create_list_as_header(103,200,400,1,10)


if __name__ == "__main__":
    parse_inputs_then_run()