from sys import argv

intro_str = "Call this file with:\npython post_processor.py <nl> -datafiles <dt1> <dt2> ... -benchmarkfiles <bn1> <bn2> ...\nWhere:\n\n"

nl_str = "- nl is the number of 64-bit limbs representing each number in the datafiles.\n  It is identical to the value of LIMBNUM in limbnum.h when the benchmarks where run.\n"
dt_str = "- dt1, dt2, ... are names of files in the same directory as this script. \n  They are all data files containing lattices generated for the same value of LIMBNUM given in nl.\n  At least one such data file is required.\n"
bn_str = "- bn1, bn2, ... are names of files in the same directory as this script.\n  They must contain outputs of bench_stats on dt1, dt2, ... in this order, and there must be as many such fils as data files.\n\n"
helpstr = intro_str + nl_str + dt_str + bn_str

def parse_arguments():
    if(len(argv)==1) and (argv[0]=="-h"):
        print(helpstr)
    else: 
        nl = 0
        dts = []
        bns = []
        try:
            nl = int(argv[1],10)
            assert(argv[2]=="-datafiles")
            i = 3
            while(argv[i]!="-benchmarkfiles"):
                i = i+1
            dts = argv[3:i]
            i = i-2
            bns = argv[i+3:]
            assert(len(bns)==len(dts))
            assert(len(bns)>0)
            return(nl,dts,bns)
        except:
            print("Wrong input. Use -h to get more information.\n")
            exit(1)


def parse_datafile(dt, nl):
    bytes_per_int = 8*nl
    bytes_per_lattice = bytes_per_int*17
    try:
        f = open(dt,"rb")
        lats_bytes = f.read()
        f.close()
        lattice_number = len(lats_bytes)//bytes_per_lattice
        l = []
        for i in range(lattice_number):
            index = i*bytes_per_lattice+bytes_per_int
            coeff1_bytes = lats_bytes[index: index+bytes_per_int]
            coeff1 = int.from_bytes(coeff1_bytes,"little")
            l.append(coeff1)
        return(l)
    except: 
        print("Failure at opening or reading file {}\n".format(dt))
        exit(1)

def parse_benchmarkfile(bn):
    try:
        f = open(bn,"r")
        bn_lines = f.readlines()
        f.close()
        bn_clean_lines = [b.strip() for b in bn_lines]
        bn_measures = [int(n,10) for n in bn_clean_lines]
        return(bn_measures)
    except: 
        print("Failure at opening or reading file {}\n".format(bn))
        exit(1)

def accumulate_data(nl,dts,bns):
    n = len(dts)
    acc = [[] for i in range(n)]
    for i in range(n):
        dt = parse_datafile(dts[i],nl)
        bn = parse_benchmarkfile(bns[i])
        assert(len(dt)==len(bn))
        acc[i] = [(dt[j],bn[j]) for j in range(len(bn))]
    res = []
    for i in range(n):
        res = res + acc[i]
    return(res)

def tuple_first(a):
    a1,_ = a
    return(a1)

def order_data(acc):
    res = sorted(acc,key=tuple_first)
    return(res)

def label_data(acc):
    ordered = order_data(acc)
    middle_index = (len(acc)//2)-1
    median = ordered[middle_index][0]
    print(ordered)
    print(median)
    small = [b for a,b in ordered if a <=median]
    large = [b for a,b in ordered if a > median]
    return(small,large)

def print_labeled_to_csv(small,large,filename):
    try:
        f = open(filename,"w")
        f.write("V1,V2\n")
        for m in small:
            f.write("X,{}\n".format(m))
        for m in large:
            f.write("Y,{}\n".format(m))
        f.close()
    except: 
        print("Failed writing to {}".format(filename))
        exit(1)
    return

def run():
    nl,dts,bns = parse_arguments()
    acc = accumulate_data(nl,dts,bns)
    small,large = label_data(acc)
    print_labeled_to_csv(small,large,"measures.csv")
    return

    
if __name__=="__main__":
    run()