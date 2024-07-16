from serializer import set_bkz_constants
from sys import argv

intro = "Expected inputs: bkz_tours(int) lagrange_tours(int)\n"
helptext = intro

def parse_inputs_and_create_header():
    if(len(argv)>1):
        if(argv[1]=="-h"):
            print(helptext)
            return
    data = False
    v = argv
    if(not len(v)==3):
        print("Wrong number of inputs, use -h for help\n")
        return
    bkz_tours = int(v[1],10)
    lagrange_tours = int(v[2],10)
    set_bkz_constants(bkz_tours, lagrange_tours)


if __name__ == "__main__":
    parse_inputs_and_create_header()


