#ifndef CTHELPERS_H
#define CTHELPERS_H

#include "../ct_intbig/ct_intbig.h"

typedef struct quat_alg {
    ibz_t p; ///< Prime number, must be = 3 mod 4.
    ibz_mat_4x4_t gram; ///< Gram matrix of the norm form
} quat_alg_t;

typedef struct quat_alg_elem {
  ibz_t denom; ///< Denominator by which all coordinates are divided (big integer, must not be 0)
  ibz_vec_4_t coord; ///< Numerators of the 4 coordinates of the quaternion algebra element in basis (1,i,j,ij)
} quat_alg_elem_t;


void quat_alg_init_set(quat_alg_t *alg, const ibz_t *p);
void quat_alg_init_set_ui(quat_alg_t *alg, unsigned int p);
void quat_alg_finalize(quat_alg_t *alg);

void quat_alg_elem_init(quat_alg_elem_t *elem);
void quat_alg_elem_finalize(quat_alg_elem_t *elem);

void quat_lattice_init(quat_lattice_t *lat);
void quat_lattice_finalize(quat_lattice_t *lat);

void ibq_round(ibz_t *rounded, const ibq_t *q);

//bilinear form associated to the norm of alg
void ibz_vec_4_bilinear(ibz_t *prod, const ibz_vec_4_t *a, const ibz_vec_4_t *b, const quat_alg_t *alg);

// set integer vector coefficients to given integers
void ibz_vec_4_set(ibz_vec_4_t *vec, int a, int b, int c, int d);

// set integer matrix coefficients to given integers, in order 00,01,10,11
void quat_test_helper_ibz_mat_2x2_set(ibz_mat_2x2_t *mat, int a, int b, int c, int d);

// set rational to fraction of two small integers
void quat_test_helper_ibq_set_i(ibq_t *q, int n, int d);

// test equality of integer vectors, return 0 if equal, 1 if not
int quat_test_helper_ibz_vec_4_equal(const ibz_vec_4_t *vec, const ibz_vec_4_t *cmp);

// compares integer to small int via ibz_cmp
int quat_test_helper_ibz_equal_i(const ibz_t *x, int cmp);

// compares integer to fraction of small ints via ibq_cmp
int quat_test_helper_ibq_equal_i(const ibq_t *x, int num, int denom);

//return 0 if equal, 1 otherwise
int ibz_mat_4x4_equal(const ibz_mat_4x4_t *a, const ibz_mat_4x4_t *b);

//return 0 if equal, 1 otherwise
int ibq_mat_4x4_equal(const ibq_mat_4x4_t *a, const ibq_mat_4x4_t *b);

//return 0 if reduced enough, 1 otherwise
int quat_test_helper_is_reduced(const ibz_mat_4x4_t *mat, const ibz_t *norm_ideal, const quat_alg_t *alg);

//Set to given ints
void ibz_vec_4_set_ibz(ibz_vec_4_t *vec, int64_t coord0, int64_t coord1, int64_t coord2, int64_t coord3);

//returns 1 if norm is integer, 0 else
int quat_alg_norm(ibz_t *norm,const quat_alg_elem_t *elem, const quat_alg_t *alg);

int quat_lattice_ideal_norm(ibz_t *norm, const quat_lattice_t *lat, const quat_alg_t *alg);
#endif