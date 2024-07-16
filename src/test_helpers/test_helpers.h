
#ifndef TESTHELPERS_H
#define TESTHELPERS_H

//should remove this depedency and add the new ints
#include "nct_helpers.h"

//extended gcd, not ct.
void nct_ibz_xgcd_non_ct(nct_ibz_t *gcd, nct_ibz_t *u, nct_ibz_t *v, const nct_ibz_t *n, const nct_ibz_t *d);

//return 0 if correct, 1 otherwise. Does not check intervals of coefficients
int nct_quat_test_xgcd_verify(const nct_ibz_t *gcd, const nct_ibz_t *u, const nct_ibz_t *v, const nct_ibz_t *d, const nct_ibz_t *n);


int nct_ibz_vec_4_scalar_div(nct_ibz_vec_4_t *quot, const nct_ibz_t *scalar, const nct_ibz_vec_4_t *vec);
int nct_ibz_mat_4x4_scalar_div(nct_ibz_mat_4x4_t *quot, const nct_ibz_t *scalar, const nct_ibz_mat_4x4_t *mat);

//copied from SQIsign

//linear combination
void nct_ibz_vec_4_linear_combination(nct_ibz_vec_4_t *lc, const nct_ibz_t *coeff_a, const nct_ibz_vec_4_t  *vec_a, const nct_ibz_t *coeff_b, const nct_ibz_vec_4_t *vec_b);

//Small helper for integers
void nct_ibz_mod_not_zero(nct_ibz_t *res, const nct_ibz_t *x, const nct_ibz_t *mod);

//centered and rather positive then negative
void nct_ibz_centered_mod(nct_ibz_t *remainder, const nct_ibz_t *a, const nct_ibz_t *mod);

// if c, res = x, else res = y
void nct_ibz_conditional_assign(nct_ibz_t *res, const nct_ibz_t *x, const nct_ibz_t *y, int c);

//understand how to put the non-zero on the 1st one
void nct_ibz_xgcd_with_u_not_0(nct_ibz_t *d, nct_ibz_t *u, nct_ibz_t *v, const nct_ibz_t *x, const nct_ibz_t *y);

void nct_ibz_rounded_div(nct_ibz_t *q, const nct_ibz_t *a, const nct_ibz_t *b);

// HNF functions
int nct_ibz_mat_4x4_is_hnf(const nct_ibz_mat_4x4_t *mat);

//centered mod
void nct_ibz_vec_4_linear_combination_mod(nct_ibz_vec_4_t *lc, const nct_ibz_t *coeff_a, const nct_ibz_vec_4_t  *vec_a, const nct_ibz_t *coeff_b, const nct_ibz_vec_4_t *vec_b, const nct_ibz_t *mod);

void nct_ibz_vec_4_copy_mod(nct_ibz_vec_4_t *res, const nct_ibz_vec_4_t *vec, const nct_ibz_t *mod); 

// no need to center this, and not 0
void nct_ibz_vec_4_scalar_mul_mod(nct_ibz_vec_4_t *prod, const nct_ibz_t *scalar, const nct_ibz_vec_4_t *vec, const nct_ibz_t *mod);

//Algorithm used is the one at number 2.4.8 in Henri Cohen's "A Course in Computational Algebraic Number Theory" (Springer Verlag, in series "Graduate texts in Mathematics") from 1993
// assumes nct_ibz_xgcd outputs u,v which are small in absolute value (as described in the book)
void nct_ibz_mat_4xn_hnf_mod_core(nct_ibz_mat_4x4_t *hnf, int generator_number, const nct_ibz_vec_4_t (*generators)[generator_number], const nct_ibz_t *mod);

//take over whole algebra.c file and slightly adapt, removed useless functions
void nct_quat_alg_equal_denom(nct_quat_alg_elem_t *res_a, nct_quat_alg_elem_t *res_b, const nct_quat_alg_elem_t *a, const nct_quat_alg_elem_t *b);
void nct_quat_alg_sub(nct_quat_alg_elem_t *res, const nct_quat_alg_elem_t *a, const nct_quat_alg_elem_t *b);
void nct_quat_alg_normalize(nct_quat_alg_elem_t *x);
int nct_quat_alg_elem_is_zero(const nct_quat_alg_elem_t *x);
int nct_quat_alg_elem_equal(const nct_quat_alg_elem_t *a, const nct_quat_alg_elem_t *b);
void nct_quat_alg_elem_copy_ibz(nct_quat_alg_elem_t *elem, const nct_ibz_t *denom, const nct_ibz_t *coord0,const nct_ibz_t *coord1,const nct_ibz_t *coord2,const nct_ibz_t *coord3);
void nct_quat_alg_elem_set(nct_quat_alg_elem_t *elem, int64_t denom, int64_t coord0, int64_t coord1, int64_t coord2, int64_t coord3);


//copied over whole lattice.c file from NIST code and slightly adapted, removed useless functions
int nct_quat_lattice_equal(const nct_quat_lattice_t *lat1, const nct_quat_lattice_t *lat2);
int nct_quat_lattice_inclusion(const nct_quat_lattice_t *sublat, const nct_quat_lattice_t *overlat);
void nct_quat_lattice_reduce_denom(nct_quat_lattice_t *reduced, const nct_quat_lattice_t *lat);
void nct_quat_lattice_add(nct_quat_lattice_t *res, const nct_quat_lattice_t *lat1, const nct_quat_lattice_t *lat2);
int nct_quat_lattice_contains_without_alg(nct_ibz_vec_4_t *coord, const nct_quat_lattice_t *lat, const nct_quat_alg_elem_t *x);
void nct_quat_lattice_index(nct_ibz_t *index, const nct_quat_lattice_t *sublat, const nct_quat_lattice_t *overlat);
void nct_quat_lattice_hnf(nct_quat_lattice_t *lat);

// Additional functions specifically for tests
void nct_quat_lattice_rfactor(nct_ibq_t *rfactor, nct_ibq_t *b0norm, const nct_ibz_mat_4x4_t *mat, const nct_quat_alg_t *alg);
int nct_quat_test_length(const nct_quat_lattice_t *red, const nct_quat_alg_t *alg, int print_flag);
int nct_quat_lattice_ideal_norm(nct_ibz_t *norm, const nct_quat_lattice_t *lat, const nct_quat_alg_t *alg);

#endif
