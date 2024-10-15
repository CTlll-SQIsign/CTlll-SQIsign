import matplotlib.pyplot as plt

def read_csv(filename):
    f = open(filename,"r")
    lines = f.readlines()[1:]
    f.close()
    return(lines)

def plot_uncolored(infile, outfile):
    lines = read_csv(infile)
    x = [int(l.split(",")[1].strip(),10)/1000000000 for l in lines]
    y = [i for i in range (len(x))]
    _, ax = plt.subplots(figsize=(10, 7))  
    ax.set_ylabel("Gigacycles",fontsize=20)         
    ax.scatter(y, x,marker="x")  
    plt.show()
    plt.savefig(outfile)


def plot_uncolored(infile, outfile):
    lines = read_csv(infile)
    x = [int(l.split(",")[1].strip(),10)/1000000000 for l in lines]
    y = [i for i in range (len(x))]
    _, ax = plt.subplots(figsize=(10, 7))  
    ax.set_ylabel("Gigacycles",fontsize=20)         
    ax.scatter(y, x,marker="x")
    plt.savefig(outfile)

def plot_uncolored_double(infile1, infile2, outfile,fixed_scale=False):
    lines1 = read_csv(infile1)
    lines2 = read_csv(infile2)
    y1 = [int(l.split(",")[1].strip(),10)/1000000000 for l in lines1]
    C1 = [ "orange" for _ in lines1]
    y2 = [int(l.split(",")[1].strip(),10)/1000000000 for l in lines2]
    C2 = ["green"  for _ in lines2]
    assert(len(y1)==len(y2))
    x = [i for i in range (len(y1))]
    _, ax = plt.subplots(figsize=(20, 8))
    if(fixed_scale):
        ax.set_autoscaley_on(False)
        ax.set_yticks([0,1,2,3,4,5,6,7])
        ax.set_xticks([0,40000,80000,120000, 160000])
        ax.tick_params(labelsize=16)
        ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
        ax.set_xlabel("Increasing $N_1$",fontsize=20)  
    ax.set_ylabel("Runtime in Gigacycles",fontsize=20)         
    ax.scatter(x, y1,c=C1,marker="x")        
    ax.scatter(x, y2,c=C2,marker="x")
    plt.savefig(outfile)

def plot_colored(infile, outfile, fixed_scale=False):
    lines = read_csv(infile)
    y = [int(l.split(",")[1].strip(),10)/1000000000 for l in lines]
    C = ["green" if l.split(",")[0].strip()=="X" else "orange" for l in lines]
    x = [i for i in range (len(y))]
    _, ax = plt.subplots(figsize=(20, 8))
    if(fixed_scale):
        ax.set_autoscaley_on(False)
        ax.set_yticks([0,1,2,3,4,5,6,7])
        ax.set_xticks([0,40000,80000,120000, 160000])
        ax.tick_params(labelsize=16)
        ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
        ax.set_xlabel("Increasing $N_1$",fontsize=20) 
    ax.set_ylabel("Runtime in Gigacycles",fontsize=20)         
    ax.scatter(x, y,c=C,marker="x")
    plt.savefig(outfile)

def plot_colored_double(infile1, infile2, outfile, fixed_scale=False):
    lines1 = read_csv(infile1)
    lines2 = read_csv(infile2)
    y1 = [int(l.split(",")[1].strip(),10)/1000000000 for l in lines1]
    C1 = ["green" if l.split(",")[0].strip()=="X" else "orange" for l in lines1]
    y2 = [int(l.split(",")[1].strip(),10)/1000000000 for l in lines2]
    C2 = ["blue" if l.split(",")[0].strip()=="X" else "red" for l in lines2]
    assert(len(y1)==len(y2))
    x = [i for i in range (len(y1))]
    _, ax = plt.subplots(figsize=(20, 8))
    if(fixed_scale):
        ax.set_autoscaley_on(False)
        ax.set_yticks([0,1,2,3,4,5,6,7])
        ax.set_xticks([0,40000,80000,120000, 160000])
        ax.tick_params(labelsize=16)
        ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
        ax.set_xlabel("Increasing $N_1$",fontsize=20)
    ax.set_ylabel("Runtime in Gigacycles",fontsize=20)         
    ax.scatter(x, y1,c=C1,marker="x")        
    ax.scatter(x, y2,c=C2,marker="x")
    plt.savefig(outfile)

if __name__=="__main__":
    plot_colored("analysis/normtest/normtest_bkz.csv","plotBKZ.png",True)
    plot_colored("analysis/normtest/normtest_lll.csv","plotLLL.png",True)
    plot_colored_double("analysis/normtest/normtest_bkz.csv","analysis/normtest/normtest_lll.csv","plotBoth.png",True)
    plot_uncolored_double("analysis/normtest/normtest_bkz.csv","analysis/normtest/normtest_lll.csv","plotBothUncolored.png",True)
