
#include <stdio.h>
#include <assert.h>
#include "test_code.h"

void ct_to_nct_integer_translation(nct_ibz_t *nct, const ibz_t *ct){
    char str_ct[PRINTSIZE];
    ibz_t a;
    ibz_init(&a);
    ibz_abs(&a,ct);
    ibz_to_str(str_ct,&a);
    nct_ibz_set_from_str(nct,str_ct,10);
    if(ibz_cmp(&a,ct)!=0) nct_ibz_neg(nct,nct);
    ibz_finalize(&a);
}

int nct_quat_lattice_check_all(nct_quat_lattice_t *red, const nct_quat_lattice_t *inpt, const nct_ibz_t *norm, const nct_quat_alg_t *alg, int print_flag){
    int res, res_min, res_len,res_eq;
    nct_ibz_t n;
    nct_ibz_init(&n);
    nct_ibz_copy(&(red->denom),&((*inpt).denom));
    nct_ibz_mul(&n,norm,&(red->denom));
    nct_ibz_mul(&n,&n,&(red->denom));
    res_min =  nct_quat_test_helper_is_reduced(&(red->basis),&n,alg,print_flag);
    if (res_min&&!print_flag) printf("Failed test for Minkowski bound\n");
    res_len = nct_quat_test_length(red,alg,print_flag);
    if (res_len&&!print_flag) printf("Failed vector length test\n");
    nct_quat_lattice_hnf(red);
    res_eq =  !nct_quat_lattice_equal(red,inpt);
    if (res_eq&&!print_flag) printf("Failed equality test\n");
    nct_ibz_finalize(&n);
    res = res_len | res_eq | res_min;
    return(res);
}

int translate_and_test_reduction(const quat_lattice_t *red, const quat_lattice_t *inpt, const quat_alg_t *alg,int print_flag){
    nct_quat_alg_t n_alg;
    nct_quat_lattice_t n_lat;
    nct_quat_lattice_t n_red;
    nct_ibz_t n_norm;
    nct_ibz_init(&n_norm);
    nct_quat_lattice_init(&n_lat);
    nct_quat_lattice_init(&n_red);
    ct_to_nct_integer_translation(&n_norm,&(alg->p));
    nct_quat_alg_init_set(&n_alg,&n_norm);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ct_to_nct_integer_translation(&(n_red.basis[i][j]),&(red->basis[i][j]));
            ct_to_nct_integer_translation(&(n_lat.basis[i][j]),&(inpt->basis[i][j]));
        }
    }
    ct_to_nct_integer_translation(&(n_red.denom),&(red->denom));
    ct_to_nct_integer_translation(&(n_lat.denom),&(inpt->denom));
    nct_quat_lattice_ideal_norm(&n_norm,&n_lat,&n_alg);
    int res = nct_quat_lattice_check_all(&n_red,&n_lat,&n_norm,&n_alg,print_flag);
    nct_ibz_finalize(&n_norm);
    nct_quat_lattice_finalize(&n_lat);
    nct_quat_lattice_finalize(&n_red);
    nct_quat_alg_finalize(&n_alg);
    return(res);
}

int quat_test_lll_on_lattice(const quat_lattice_t *lat,const quat_alg_t *alg){
    quat_lattice_t red;
    int res;
    quat_lattice_init(&red);
    ibz_copy(&(red.denom),&(lat->denom));
    quat_lattice_lll(&(red.basis),lat,&(alg->p),1);
    res = translate_and_test_reduction(&red,lat,alg,0);
    quat_lattice_finalize(&red);
    return(res);
}

int quat_test_bkz_on_lattice(const quat_lattice_t *lat,const quat_alg_t *alg, int bkz_iterations, int lagrange_iterations, int print_flag){
    int res;
    quat_lattice_t red;
    quat_lattice_init(&red);
    quat_lattice_bkz(&(red.basis),lat,bkz_iterations,lagrange_iterations,alg);
    res = translate_and_test_reduction(&red,lat,alg,print_flag);
    quat_lattice_finalize(&red);
    return(res);
}
