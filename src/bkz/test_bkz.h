#ifndef BKZ_TEST_H
#define BKZ_TEST_H

#include "bkz.h"

int quat_bkz_matrix_equal(const quat_bkz_matrix_t *a,const quat_bkz_matrix_t *b);

//return 0 if equal to bkz matrix setup from is basis, 1 otherwise
int quat_test_helper_bkz_test_integrity(const quat_bkz_matrix_t *g, const quat_alg_t *alg);


//tests

//void quat_bkz_matrix_init(bkz_matrix_t *g);
//void quat_bkz_matrix_finalize(bkz_matrix_t *g);
int quat_test_bkz_matrix_finit();

//void ibq_round(ibz_t *rounded, ibq_t *q);
int quat_test_bkz_ibq_round();

//void ibz_vec_2x4_mul(ibz_vec_4_t *a, ibz_vec_4_t *b,const ibz_mat_2x2_t *U);
int quat_test_bkz_vec_2x4_mul();

//void ibz_vec_4_bilinear(ibz_t *prod, const ibz_vec_4_t *a, const ibz_vec_4_t *b, const quat_alg_t *alg);
int quat_test_bkz_bilinear();

//void quat_bkz_matrix_set(quat_bkz_matrix_t *g, const ibz_mat_4x4_t *mat, quat_alg_t *alg);
int quat_test_bkz_set_matrix();

//void quat_bkz_update_matrix(quat_bkz_matrix_t *g, int start, int end, quat_alg_t *alg);
int quat_test_bkz_update_matrix();

//void quat_bkz_size_reduce(quat_bkz_matrix_t *g, int start, int end, quat_alg_t *alg);
int quat_test_bkz_size_reduce();

//void quat_bkz_update_after_lagrange(quat_bkz_matrix_t *g, const ibz_mat_2x2_t *U, int index);
int quat_test_bkz_update_after_lagrange();

//void quat_bkz_lagrange_reduction_gram(ibz_mat_2x2_t *U, ibz_mat_2x2_t *G, int lagrange_tours);
int quat_test_bkz_lagrange_reduction_gram();

//void quat_lattice_bkz(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, int bkz_tours, int lagrange_tours, quat_alg_t *alg);
int quat_test_bkz_lattice_bkz();

//run all of the above tests
int bkz_tests();

#endif