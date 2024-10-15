from shuffledata import parse_all_data
from copy import deepcopy

def parse_benchmark(filename):
    try:
        f = open(filename,"rb")
        l = f.readlines()
        lints = [int(i.strip(),10) for i in l]
        f.close()
        assert(len(lints) == 1200)
        return(lints)
    except: 
        print("Failure at opening or reading file {}\n".format(filename))
        exit(1)

def parse_all_benchmarks(foldername,fileprefix):
    l = [[parse_benchmark(foldername+"/" + str(i)+"/"  + fileprefix + str(j)) for j in range(1,13)] for i in range (0,10)]
    return(l)

def align_datasets_and_bench(data,bench,idx1,idx2):
    assert(len(data)==len(bench))
    assert(len(data[0])==len(bench[0]))
    assert(len(data[0])>idx2)
    assert(len(data[0])>idx1)
    res = []
    for i in range(10):
        for j in range(idx1,idx2+1):
            res = res + list(zip(data[i][j],bench[i][j]))
    assert(len(res)==1200*(idx2-idx1+1)*10)
    return(res)

def tuple_first(t):
    a,b = t
    return(a)

def sort_by_first(l):
    return(sorted(l,key=tuple_first))


def tuple_second(t):
    a,b = t
    return(b)

def sort_on_second(l):
    return(sorted(l,key=tuple_second))

def align_benchs(zip1,zip2):
    z1 = sort_by_first(zip1)
    z2 = sort_by_first(zip2)
    r1 = [t[1] for t in z1]
    r2 = [t[1] for t in z2]
    r = zip(r1,r2)
    return(r)

def label_sets(matched):
    s = sort_by_first(matched)
    idx = len(s)//2
    l1 = [("X",t[1]) for t in s[:idx]]
    l2 = [("Y",t[1]) for t in s[idx:]]
    res = l1+l2
    assert(len(res)==1200*10*10)
    return(res)   

def writeCSV(l ,filename):
    h = "V1,V2\n"
    for a,b in l:
        h = h + a + "," + str(b) + "\n"
    f = open(filename,"w")
    f.write(h)
    f.close()
    return

def reorder_as_measure(l,b2):
    b = []
    for i in range (0,10):
        for j in range(2,12):
            b = b+b2[i][j]
    z = [(i,b[i]) for i in range(len(b))]
    zo = sort_on_second(z)
    lo = sort_on_second(l)
    assert(len(lo)==len(zo))
    for i in range(len(lo)):
        assert(lo[i][1]==zo[i][1])
    group = [(zo[i][0],lo[i]) for i in range (len(zo))]
    ordered = sort_by_first(group)
    uncount = [o[1] for o in ordered]
    return(uncount)

def all(data1,bench1,data2,bench2,benchprefix,dataprefix, idx1,idx2,nl, outfile):
    b1 = parse_all_benchmarks(bench1,benchprefix)
    b2 = parse_all_benchmarks(bench2,benchprefix)
    d1 = parse_all_data(data1,dataprefix,nl)
    d2 = parse_all_data(data2,dataprefix,nl)
    zip1 = align_datasets_and_bench(d1,b1,idx1,idx2)
    zip2 = align_datasets_and_bench(d2,b2,idx1,idx2)
    a = align_benchs(zip1,zip2)
    l = label_sets(a)
    o = reorder_as_measure(l,b2)
    writeCSV(o,outfile)
    return

if __name__=="__main__":
    all("../datafiles","shuffletest_data/results_10_threads_unshuffled","shuffletest_data/shuffled_data","shuffletest_data/results_10_threads_shuffled","bench_p335_nl11_","p335_nl11_",2,11,11,"shuffletest.csv")
