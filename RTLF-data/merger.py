
#indices = [i for i in range(2, 12)]# + [i for i in range(21, 37)]

indices = [i for i in range(0, 10)] #indices of files in each folder

for ctr in [i for i in range(0, 10)]: #folder counter

    folder = "results_lll_10_cores/"+str(ctr)+"/"
    #folder = "results_bkz_10_cores_2/all/"  #uncomment for the last merge


    #merge data files
    for i in indices:
        filename_lattices = folder+"p335_nl11_"+str(i)+"_all"
        #with open(folder+"p335_nl11_"+str(ctr)+"_all", 'ab') as out, open(filename_lattices, 'rb') as inp:
        with open(folder+"p335_nl11_all_all", 'ab') as out, open(filename_lattices, 'rb') as inp:
            out.write(inp.read())


    #merge benchmarks
    benchmarks = ""
    for i in indices:
        filename_benchmarks = folder+"bench_p335_nl11_"+str(i)+"_all"
        with open(filename_benchmarks, 'r') as fp:
            benchmarks += fp.read()

    #with open(folder+"bench_p335_nl11_"+str(ctr)+"_all", 'w') as fp:
    with open(folder+"bench_p335_nl11_all_all", 'w') as fp:
        fp.write(benchmarks)




