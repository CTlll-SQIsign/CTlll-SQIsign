import matplotlib.pyplot as plt


bt = [i for i in range(3,16)]
lt = [i for i in range(5,13)]

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
    for bt_ in bt:
        for lt_ in lt:
            ar1, ar2 = process_file(prefix, bt_, lt_)
            #print(bt_, lt_, len(ar1), len(ar2))
            assert(len(ar1)==ln)
            assert(len(ar2)==ln)

            #ar1 = array of Minkowski 2nd thm. check. As it is the product of 4 norms, we take the 4th root
            av1 = compute_average_approx_factor(ar1)
            fail1 = num_failures(ar1)

            #ar1 = array of Lemma 6 raised to 8th power, hence the 8th root
            av2 =  compute_average_approx_factor(ar2)
            fail2 = num_failures(ar2)

            all_data.append([bt_, lt_, ar1, ar2])
            stat_data.append([bt_, lt_, av1**(1./4), fail1,  av2**(1./8), fail2])
    return all_data, stat_data

def make_csv_from_stats(stat_data):
    """
    :param stat_data: array of tuples [bt, lt, rhs/lhs(mink), num_failures(mink),  rhs/lhs(lem6), num_failures(lem6)]
    :return: csv files for each bt of the form lt, rhs/lhs(mink), num_failures(mink)
    """
    ctr_out = 0
    for bt_ in bt:
        filename='minkstatslvl3'+str(bt_)+'.csv'
        with open(filename, 'w') as fp:
            fp.write("lt,gamma,fail\n")
            for ctr_in in range(len(lt)):
                i = ctr_out*len(lt)+ctr_in
                fp.write(str(stat_data[i][1])+','+str(stat_data[i][2])+','+str(stat_data[i][3])+'\n')
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
    #all, stats = process_all_files("lvl5/lvl5_p239_nl37")
    #all, stats = process_all_files("lvl1_p239_nl20")
    all, stats = process_all_files("lvl3_p596_nl28")
    make_plots((stats))
    make_csv_from_stats(stats)













