#include "../ct_intbig/ct_intbig.h"
#include "lll.h"
#include <stdio.h>


//int quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const ibz_t *q, int precision);
int test_lattice_lll(){
    int res = 0;
    quat_lattice_t lat, test;
    ibz_mat_4x4_t red;
    ibz_t num, denom, q;
    ibq_t coeff;
    ibz_init(&num);
    ibz_init(&denom);
    ibz_init(&q);
    ibq_init(&coeff);
    ibz_mat_4x4_init(&red);
    ibz_mat_4x4_init(&lat.basis);
    ibz_init(&lat.denom);

    // set lattice
    ibz_set(&lat.denom, 60);
    ibz_set(&lat.basis[0][0], 3);
    ibz_set(&lat.basis[1][0], 7);
    ibz_set(&lat.basis[0][1], 1);
    ibz_set(&lat.basis[3][1], -6);
    ibz_set(&lat.basis[1][2], 12);
    ibz_set(&lat.basis[2][2], 5);
    ibz_set(&lat.basis[0][3], -19);
    ibz_set(&lat.basis[3][3], 3);
    
    ibz_set(&q,103);
    res = res || quat_lattice_lll(&red,&lat,&q,10);
    
    if (res != 0){
        printf("    Test lattice_lll failed\n");
    }
    ibz_finalize(&num);
    ibz_finalize(&denom);
    ibz_finalize(&q);
    ibq_finalize(&coeff);
    ibz_mat_4x4_finalize(&red);
    ibz_mat_4x4_finalize(&lat.basis);
    ibz_finalize(&lat.denom);
    return(res);
}


int lll_tests(){
    printf("Run lll test:\n");
    int res = test_lattice_lll();
    if(res){
        printf("LLL test failed\n");
    } else {
        printf("All passed\n");
    }
    return(res);
}
