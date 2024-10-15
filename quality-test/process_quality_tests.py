import matplotlib.pyplot as plt


bt = [i for i in range(1,21)]
lt = [i for i in range(1,13)]

btlt = [0]*(len(bt))
ctr = 0
for i in range(len(bt)):
    btlt[i] = [0]*(len(lt))
    for j in range(len(lt)):
        (btlt[i])[j] = (i+1, j+1)

ln = 30 #number of lattices

def make_int(s):
    res = s.split("/")
    if (len(res) == 1): return int(res[0])
    else:
        nom = int(res[0])
        denom = int(res[1])
    return nom//denom


def process_file(prefix, bt, lt):
    data_mink = []
    data_lem6 = []
    filename=prefix+'_bt'+str(bt)+'_lt'+str(lt)
    file = open(filename,'r')
    filedata = file.readlines()

    for i in range(1, len(filedata)-2, 3):
        if filedata[i].strip()=="Minkowski":
            rhs_s = (filedata[i+2].split("=")[1]).strip()
            rhs_v = make_int((rhs_s))
            lhs_s = (filedata[i+1].split("=")[1]).strip()
            lhs_v = make_int((lhs_s))
            data_mink.append([rhs_v, lhs_v])
        elif filedata[i].strip()=="Lemma6":
            rhs_s = (filedata[i+1].split("=")[1]).strip()
            rhs_v = make_int((rhs_s))
            lhs_s = (filedata[i+2].split("=")[1]).strip()
            lhs_v = make_int((lhs_s))
            data_lem6.append([rhs_v, lhs_v])
        else:
            print("Error in reading ", filename)
            return [],[]
    return data_mink, data_lem6

def compute_average_approx_factor(ar):
    """
    :param ar: array of pairs [rhs, lhs]
    :return: average of rhs/lhs
    """
    approx = 0
    ctr = 0
    for el in ar:
        if el[1]<=el[0]:
            approx+=el[0]/el[1]
        ctr = ctr+1

    return approx/(ctr)


def num_failures(ar):
    n = 0
    for el in ar:
        if el[1]>el[0]: n = n+1
    return n


def process_all_files(prefix):

    all_data = []
    stat_data = []
    for list_ in btlt:
        for (bt_, lt_) in list_:
            ar1, ar2 = process_file(prefix, bt_, lt_)
            assert(len(ar1)==ln)
            assert(len(ar2)==ln)

            #ar1 = array of Minkowski 2nd thm. check. As it is the product of 4 norms, we take the 4th root
            #av1 = compute_average_approx_factor(ar1)
            fail1 = num_failures(ar1)

            #ar1 = array of Lemma 6 raised to 8th power, hence the 8th root
            #av2 =  compute_average_approx_factor(ar2)
            fail2 = num_failures(ar2)

            all_data.append([bt_, lt_, ar1, ar2])
            #stat_data.append([bt_, lt_, av1**(1./4), fail1,  av2**(1./8), fail2])
            stat_data.append([bt_, lt_, fail1, fail2])
    return all_data, stat_data

def make_csv_from_stats(stat_data):
    """
    :param stat_data: array of tuples [bt, lt,  num_failures(mink),   num_failures(lem6)]
    :return: csv files for each bt of the form lt, num_failures(mink)
    """
    ctr_out = 0
    i = 0
    prev_len = len(btlt[0])
    for list_ in btlt:
        filename='lvl1results/test/lvl1_p239_nl20_b'+str(list_[0][0])+'.csv'
        with open(filename, 'w') as fp:
            fp.write("lt,fail\n")
            ctr_in = 0
            for el in list_:
                fp.write(str(el[1])+','+str(stat_data[i][2])+'\n')
                i+=1
                ctr_in+=1
        prev_len = len(btlt[ctr_out])
        ctr_out = ctr_out +1


def make_plots(stat_data):
    ctr_out = 0
    for bt_ in bt:
        xpoints = []
        ypoints = []
        for ctr_in in range(len(lt)):
            i = ctr_out*len(lt)+ctr_in
            xpoints.append(stat_data[i][1])
            ypoints.append(stat_data[i][3])
        plt.plot(xpoints, ypoints, label='bkz '+str(bt_))
        ctr_out = ctr_out +1
    plt.legend()
    plt.show()


if __name__=="__main__":

    all, stats = process_all_files("lvl1results/lvl1_p239_nl20")
    make_csv_from_stats(stats)
