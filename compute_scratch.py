from subprocess import check_output

def get_scratch(limbnum):
    compile = "gcc get_scratch_size.c -lgmp -o build/get_scratch".split(" ")
    assert(b''==check_output(compile))
    result = check_output("./build/get_scratch "+str(limbnum),shell=True)
    text = result.decode().strip()
    lines = text.split("\n")
    mul = int(lines[1].strip().split(" ")[-1].strip())
    div = int(lines[2].strip().split(" ")[-1].strip())
    return(mul,div)
