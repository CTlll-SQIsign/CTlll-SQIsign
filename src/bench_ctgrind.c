//#include "bench.h"
#include "inputs.h"
#include "bkz/bkz.h"
#include "bkz_constants.h"
#include <stdio.h>
#include <inttypes.h>
#include "ctgrind.h"

#if defined(DATA_FILE) && defined(LATTICE_DATA_NUMBER)
#include "read_data.h"

int main(void) {
    ibz_mat_4x4_t red;
    quat_lattice_t lat;
    quat_alg_t alg;
    quat_alg_init_set(&alg,&QUATALG_P);
    quat_lattice_init(&lat);
    ibz_mat_4x4_init(&red);
    read_lattice_from_file(&lat,0);
    //quat_lattice_finalize(&lat);
    ct_poison(&lat, sizeof(quat_lattice_t));
    quat_lattice_bkz(&red, &lat,BKZ_TOURS,LAGRANGE_TOURS,&alg);
    ct_unpoison(&lat,sizeof(quat_lattice_t));
    //quat_lattice_init(&lat);
    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_alg_finalize(&alg);
}

#else
int main(int argc, char* argv[]){
    printf("input.h not correctly set\n");
    return(0);
}
#endif