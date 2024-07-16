//Compile with gcc get_scratch_size.c -lgmp
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv []){
    int limbnum = atoi(argv[1]);
    int mul = mpn_sec_mul_itch(limbnum,limbnum);
    int div = mpn_sec_div_qr_itch(limbnum,limbnum);
    printf("LIMBNUM %d\n",limbnum);
    printf("MUL_SCRATCH %d\nDIV_SCRATCH %d\n",mul,div);
    return(0);
}