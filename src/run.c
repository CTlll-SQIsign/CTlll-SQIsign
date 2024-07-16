#include "bkz/bkz.h"
#include "inputs.h"
#include "bkz_constants.h"

#if defined(LATTICE_NUMBER)

//FOR CT testing (should run in same time for fixed LATTICE_NUMBER, QUATALG_P, bkz_iterations and lagrange_iterations)
void quat_ct_bkz_on_list(ibz_mat_4x4_t *red, int lattice_number,const quat_lattice_t lat[lattice_number],const quat_alg_t *alg, int bkz_iterations, int lagrange_iterations){
    for(int i = 0; i<lattice_number;i++){
        quat_lattice_bkz(red, &(lat[i]),bkz_iterations,lagrange_iterations,alg);
    }
}


//run test on input 1, ct runnner on input 0 or no input
int main(int argc, char* argv[]){
    int lagrange_iterations = LAGRANGE_TOURS;
    int bkz_iterations = BKZ_TOURS;
    quat_alg_t alg;
    ibz_mat_4x4_t red;
    ibz_mat_4x4_init(&red);
    quat_alg_init_set(&alg,&QUATALG_P);
    quat_ct_bkz_on_list(&red,LATTICE_NUMBER,LATTICE_LIST,&alg,bkz_iterations,lagrange_iterations);
    ibz_mat_4x4_finalize(&red);
    quat_alg_finalize(&alg);
}


#elif (defined(DATA_FILE)&&defined(LATTICE_DATA_NUMBER))
#include "read_data.h"

uint64_t quat_run_bkz_from_file(ibz_mat_4x4_t *red, quat_lattice_t *lat, const quat_alg_t *alg){
    for(int i = 0; i<LATTICE_DATA_NUMBER;i++){
        read_lattice_from_file(lat,i);
        quat_lattice_bkz(red, lat,BKZ_TOURS,LAGRANGE_TOURS,alg);
    }
}

int main(int argc, char* argv[]){
    quat_alg_t alg;
    quat_lattice_t lat;
    ibz_mat_4x4_t red;
    ibz_mat_4x4_init(&red);
    quat_lattice_init(&lat);
    quat_alg_init_set(&alg,&QUATALG_P);
    quat_run_bkz_from_file(&red,&lat,&alg);
    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_alg_finalize(&alg);
    return(0);
}

#else
#include <stdio.h>
int main(int argc, char* argv[]){
    printf("input.h not correctly set for run.c\n");
    return(1);
}
#endif