import random

#this script assumes that there are 10 folders called 0-9 and 12 data files per folder noted 1 to 12 with a prefix, of which 1st and last have a special an separated role
#All files should contain the same number of lattices

FILES_USED = 10
FILES_ALL = 12
THREADS = 10

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
            index = i*bytes_per_lattice
            lat_bytes = lats_bytes[index: index+bytes_per_lattice]
            l.append(lat_bytes)
        return(l)
    except: 
        print("Failure at opening or reading file {}\n".format(dt))
        exit(1)

def parse_all_data(foldername,fileprefix,nl):
    l = [[parse_datafile(foldername+"/" + str(i)+"/"  + fileprefix + str(j), nl) for j in range(1,FILES_ALL+1)] for i in range (0,THREADS)]
    return(l)

#mixes all data between 2 and USED_FILES+1, while all data in 1 and after USED_FILES+1 is mixed separately
def mix_data(l):
    middle = []
    borders = []
    lattices_per_file = len(l[0][0])
    for i in range(0,THREADS):
        for j in range(1,FILES_USED+1):
            middle = middle + l[i][j]
        borders = borders + l[i][0]+l[i][11]
    random.shuffle(middle,random.SystemRandom().random)
    random.shuffle(borders,random.SystemRandom().random)
    res = [[]for i in range(0,THREADS)]
    for i in range(0,THREADS):
        res[i].append(borders[(FILES_ALL-FILES_USED)*i*lattices_per_file:((FILES_ALL-FILES_USED)*i+1)*lattices_per_file])
        for j in range(0,FILES_USED):
            res[i].append(middle[(FILES_USED*i+j)*lattices_per_file:(FILES_USED*i+1+j)*lattices_per_file])
        for _ in range(FILES_ALL-FILES_USED-1):
            res[i].append(borders[((FILES_ALL-FILES_USED)*i+1)*lattices_per_file:((FILES_ALL-FILES_USED)*i+(FILES_ALL-FILES_USED))*lattices_per_file])
    return(res)

def write_all_files(foldername, fileprefix, l):
    for i in range(0,THREADS):
        prefix = foldername + "/" + str(i)+"/" + fileprefix
        for j in range(0,FILES_ALL):
            f = open(prefix+str(1+j),"wb")
            text = b""
            for b in l[i][j]:
                text = text + b
            f.write(text)
            f.close()
    return


def all(prefix, datafolder, output, nl):
    l = parse_all_data(datafolder, prefix,nl)
    r = mix_data(l)
    write_all_files(output, prefix,r)

def test(prefix, datafolder, output, nl):
    li = parse_all_data(datafolder,prefix,nl)
    ri = parse_all_data(output,prefix,nl)
    l = []
    r = []
    for i in range(0,THREADS):
        for j in range(0,FILES_ALL):
            l = l + li[i][j]
            r = r + ri[i][j]
    assert(len(l)==len(r))
    assert(len(r)==THREADS*FILES_ALL*1200)
    l.sort()
    r.sort()
    for i in range(len(l)):
        assert(l[i]==r[i])
    #only internal
    li = parse_all_data(datafolder,prefix,nl)
    ri = parse_all_data(output,prefix,nl)
    l = []
    r = []
    for i in range(0,THREADS):
        for j in range(1,FILES_USED+1):
            l = l + li[i][j]
            r = r + ri[i][j]
    assert(len(l)==len(r))
    assert(len(r)==THREADS*FILES_USED*1200)
    l.sort()
    r.sort()
    for i in range(len(l)):
        assert(l[i]==r[i])
    return

if __name__ == "__main__":
    all("p335_nl11_", "../../datafiles","shuffled_data",11)
    test("p335_nl11_", "../../datafiles","shuffled_data",11)