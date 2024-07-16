
#include "ct_helpers.h"
#include <assert.h>
#include <stdio.h>

void quat_alg_init_set(quat_alg_t *alg, const ibz_t *p){
    ibz_init(&(*alg).p);
    ibz_mat_4x4_init(&(*alg).gram);
    ibz_copy(&(*alg).p, p);
    ibz_set(&(*alg).gram[0][0], 1);
    ibz_set(&(*alg).gram[1][1], 1);
    ibz_copy(&(*alg).gram[2][2], p);
    ibz_copy(&(*alg).gram[3][3], p);
}

void quat_alg_init_set_ui(quat_alg_t *alg, unsigned int p) {
    ibz_t bp;
    ibz_init(&bp);
    ibz_set(&bp, p);
    quat_alg_init_set(alg, &bp);
    ibz_finalize(&bp);
}

void quat_alg_finalize(quat_alg_t *alg){
    ibz_finalize(&(*alg).p);
    ibz_mat_4x4_finalize(&(*alg).gram);
}

void quat_alg_elem_init(quat_alg_elem_t *elem){
    ibz_vec_4_init(&(*elem).coord);
    ibz_init(&(*elem).denom);
    ibz_set(&(*elem).denom, 1);
}
void quat_alg_elem_finalize(quat_alg_elem_t *elem){
    ibz_vec_4_finalize(&(*elem).coord);
    ibz_finalize(&(*elem).denom);
}

void quat_lattice_init(quat_lattice_t *lat){
    ibz_mat_4x4_init(&(*lat).basis);
    ibz_init(&(*lat).denom);
    ibz_set(&(*lat).denom, 1);
}
void quat_lattice_finalize(quat_lattice_t *lat){
    ibz_finalize(&(*lat).denom);
    ibz_mat_4x4_finalize(&(*lat).basis);
}

void ibq_round(ibz_t *rounded, const ibq_t *q){
    ibz_t num, den;
    ibz_init(&num);
    ibz_init(&den);
    ibz_set(&den,2);
    ibq_num(&num,q);
    ibz_mul(&num,&den,&num);
    ibq_denom(&den,q);
    ibz_add(&num,&den,&num);
    ibz_div_floor(rounded,&den,&num,&den);
    ibz_set(&den,2);
    ibz_div_floor(rounded,&den,rounded,&den);
    ibz_finalize(&num);
    ibz_finalize(&den);
}

void ibz_vec_4_bilinear(ibz_t *prod, const ibz_vec_4_t *a, const ibz_vec_4_t *b, const quat_alg_t *alg){
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_mul(prod,&((*a)[0]),&((*b)[0]));
    ibz_mul(&tmp,&((*a)[1]),&((*b)[1]));
    ibz_add(prod,prod,&tmp);
    ibz_mul(&tmp,&((*a)[2]),&((*b)[2]));
    ibz_mul(&tmp,&tmp,&(alg->p));
    ibz_add(prod,prod,&tmp);
    ibz_mul(&tmp,&((*a)[3]),&((*b)[3]));
    ibz_mul(&tmp,&tmp,&(alg->p));
    ibz_add(prod,prod,&tmp);
    ibz_finalize(&tmp);
}

// set integer vector coefficients to given integers
void ibz_vec_4_set(ibz_vec_4_t *vec, int a, int b, int c, int d){
    ibz_set(&((*vec)[0]),a);
    ibz_set(&((*vec)[1]),b);
    ibz_set(&((*vec)[2]),c);
    ibz_set(&((*vec)[3]),d);
}

// set integer matrix coefficients to given integers, in order 00,01,10,11
void quat_test_helper_ibz_mat_2x2_set(ibz_mat_2x2_t *mat, int a, int b, int c, int d){
    ibz_set(&((*mat)[0][0]),a);
    ibz_set(&((*mat)[0][1]),b);
    ibz_set(&((*mat)[1][0]),c);
    ibz_set(&((*mat)[1][1]),d);
}

// set rational to fraction of two small integers
void quat_test_helper_ibq_set_i(ibq_t *q, int n, int d){
    ibz_t num,denom;
    ibz_init(&num);
    ibz_init(&denom);
    ibz_set(&num,n);
    ibz_set(&denom,d);
    ibq_set(q,&num,&denom);
    ibz_finalize(&num);
    ibz_finalize(&denom);
}

// test equality of integer vectors, eturn 0 if equal, 1 if not
int quat_test_helper_ibz_vec_4_equal(const ibz_vec_4_t *vec, const ibz_vec_4_t *cmp){
    int res = 0;
    res = res | (ibz_cmp(&((*vec)[0]),&((*cmp)[0]))!=0);
    res = res | (ibz_cmp(&((*vec)[1]),&((*cmp)[1]))!=0);
    res = res | (ibz_cmp(&((*vec)[2]),&((*cmp)[2]))!=0);
    res = res | (ibz_cmp(&((*vec)[3]),&((*cmp)[3]))!=0);
    return(res);
}

// compares integer to small int via ibz_cmp
int quat_test_helper_ibz_equal_i(const ibz_t *x, int cmp){
    ibz_t c;
    ibz_init(&c);
    ibz_set(&c,cmp);
    int res = ibz_cmp(x,&c);
    ibz_finalize(&c);
    return(res);
}

// compares integer to fraction of small ints via ibq_cmp
int quat_test_helper_ibq_equal_i(const ibq_t *x, int num, int denom){
    ibq_t cmp;
    ibq_init(&cmp);
    quat_test_helper_ibq_set_i(&cmp,num,denom);
    int res = ibq_cmp(x,&cmp);
    ibq_finalize(&cmp);
    return(res);
}

//return 0 if equal, 1 otherwise
int ibz_mat_4x4_equal(const ibz_mat_4x4_t *a, const ibz_mat_4x4_t *b){
    int res = 0;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            res = res | ibz_cmp(&((*a)[i][j]),&((*b)[i][j]));
        }
    }
    return(res!=0);
}

//return 0 if equal, 1 otherwise
int ibq_mat_4x4_equal(const ibq_mat_4x4_t *a, const ibq_mat_4x4_t *b){
    int res = 0;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            res = res | ibq_cmp(&((*a)[i][j]),&((*b)[i][j]));
        }
    }
    return(res!=0);
}

//return 0 if reduced enough, 1 otherwise
int quat_test_helper_is_reduced(const ibz_mat_4x4_t *mat, const ibz_t *norm_ideal, const quat_alg_t *alg){
    ibz_mat_4x4_t trans;
    ibz_t prod, norm, ref;
    ibz_init(&prod);
    ibz_init(&norm);
    ibz_init(&ref);
    ibz_mat_4x4_init(&trans);
    ibz_mat_4x4_transpose(&trans,mat);
    ibz_set(&prod,1);
    for (int i = 0; i < 4; i++){
        ibz_vec_4_bilinear(&(norm),&(trans[i]),&(trans[i]),alg);
        ibz_mul(&prod,&prod,&norm);
    }
    //never certain of the factors four in this story
    //prod = 16*product of norms
    //ibz_set(&norm,16);
    //ibz_mul(&prod,&prod,&norm);
    //ref = 4*p^2*NI^4
    ibz_mul(&ref,&(alg->p),&(alg->p));
    ibz_set(&norm,4);
    ibz_mul(&ref,&ref,&norm);
    ibz_mul(&norm,norm_ideal,norm_ideal);
    ibz_mul(&norm,&norm,&norm);
    ibz_mul(&ref,&norm,&ref);
    int res = (0<ibz_cmp(&prod,&ref));
    ibz_finalize(&prod);
    ibz_finalize(&norm);
    ibz_finalize(&ref);
    ibz_mat_4x4_finalize(&trans);
    return(res);
}

void ibz_vec_4_set_ibz(ibz_vec_4_t *vec, int64_t coord0, int64_t coord1, int64_t coord2, int64_t coord3){
    ibz_set(&((*vec)[0]),coord0);
    ibz_set(&((*vec)[1]),coord1);
    ibz_set(&((*vec)[2]),coord2);
    ibz_set(&((*vec)[3]),coord3);
}

