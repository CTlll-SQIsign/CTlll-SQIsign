

def read_lines(filename):
    f=open(filename,"r")
    l = f.readlines()
    f.close()
    return(l)

def get_counts(l):
    bns = []
    for line in l:
        try: 
            li = line.strip()
            halves = li.split(" ")
            if(len(halves)==2):
                b = int(halves[1].strip(),10)
                bns.append(b)
        except: 
            continue
    return(bns)

def average(l):
    s =0
    for c in l:
        s = s+c
    s = s/len(l)
    return(s)

def line_average(filename):
    l = read_lines(filename)
    bns = get_counts(l)
    avg = average(bns)
    print(avg)
    return(avg)

import sys

if __name__=="__main__":
    line_average(sys.argv[1])