#ifndef NCTHELPERS_H
#define NCTHELPERS_H

#include "../nct_intbig/nct_intbig.h"

/** @brief Type for matrices of 4x4 big integers
 * 
 * @typedef ibz_mat_4x4_t
*/
typedef nct_ibz_t nct_ibz_mat_4x4_t[4][4];

/** @brief Type for matrices of 2x2 big integers
 * 
 * @typedef ibz_mat_2x2_t
*/
typedef nct_ibz_t nct_ibz_mat_2x2_t[2][2];

/** @brief Type for vectors of 4 big integers
 * 
 * @typedef ibz_vec_4_t
*/
typedef nct_ibz_t nct_ibz_vec_4_t[4];


/** @brief Type for matrices of 4x4 rationals
 * 
 * @typedef ibq_mat_4x4_t
*/
typedef nct_ibq_t nct_ibq_mat_4x4_t[4][4];

/** @brief Type for matrices of 2x2 rationals
 * 
 * @typedef ibq_mat_4x4_t
*/
typedef nct_ibq_t nct_ibq_mat_2x2_t[2][2];

typedef struct nct_quat_alg {
    nct_ibz_t p; ///< Prime number, must be = 3 mod 4.
    nct_ibz_mat_4x4_t gram; ///< Gram matrix of the norm form
} nct_quat_alg_t;

typedef struct nct_quat_alg_elem {
  nct_ibz_t denom; ///< Denominator by which all coordinates are divided (big integer, must not be 0)
  nct_ibz_vec_4_t coord; ///< Numerators of the 4 coordinates of the quaternion algebra element in basis (1,i,j,ij)
} nct_quat_alg_elem_t;

/** @brief Type for lattices in dimension 4
 * 
 * @typedef quat_lattice_t
 * 
 * @struct quat_lattice
 * 
 * Represented as a rational (`frac`) times an integreal lattice (`basis`)
 * 
 * The basis is in hermite normal form, and its columns divided by its denominator are elements of the quaternion algebra, represented in basis (1,i,j,ij) where i^2 = -1, j^2 = -p.
 * 
 * All lattices must have full rank (4)
*/
typedef struct nct_quat_lattice {
  nct_ibz_t denom; ///< Denominator by which the basis is divided (big integer, must not be 0)
  nct_ibz_mat_4x4_t basis; ///< Integer basis of the lattice in hermite normal form (its columns divided by denom are algebra elements in the usual basis)
} nct_quat_lattice_t;



void nct_ibz_mat_4x4_init(nct_ibz_mat_4x4_t *x);
void nct_ibz_mat_2x2_init(nct_ibz_mat_2x2_t *x);
void nct_ibq_mat_2x2_init(nct_ibq_mat_2x2_t *x);
void nct_ibq_mat_4x4_init(nct_ibq_mat_4x4_t *x);
void nct_ibz_vec_4_init(nct_ibz_vec_4_t *x);

void nct_ibz_mat_4x4_finalize(nct_ibz_mat_4x4_t *x);
void nct_ibz_mat_2x2_finalize(nct_ibz_mat_2x2_t *x);
void nct_ibq_mat_2x2_finalize(nct_ibq_mat_2x2_t *x);
void nct_ibq_mat_4x4_finalize(nct_ibq_mat_4x4_t *x);
void nct_ibz_vec_4_finalize(nct_ibz_vec_4_t *x);


void nct_quat_alg_init_set(nct_quat_alg_t *alg, const nct_ibz_t *p);
void nct_quat_alg_init_set_ui(nct_quat_alg_t *alg, unsigned int p);
void nct_quat_alg_finalize(nct_quat_alg_t *alg);

void nct_quat_alg_elem_init(nct_quat_alg_elem_t *elem);
void nct_quat_alg_elem_finalize(nct_quat_alg_elem_t *elem);

void nct_quat_lattice_init(nct_quat_lattice_t *lat);
void nct_quat_lattice_finalize(nct_quat_lattice_t *lat);


void nct_ibz_mat_4x4_transpose(nct_ibz_mat_4x4_t *transposed, const nct_ibz_mat_4x4_t *mat);
void nct_ibz_mat_4x4_copy(nct_ibz_mat_4x4_t *copy, const nct_ibz_mat_4x4_t *copied);

//bilinear form associated to the norm of alg
void nct_ibz_vec_4_bilinear(nct_ibz_t *prod, const nct_ibz_vec_4_t *a, const nct_ibz_vec_4_t *b, const nct_quat_alg_t *alg);

// set integer vector coefficients to given integers
void nct_ibz_vec_4_set(nct_ibz_vec_4_t *vec, int a, int b, int c, int d);

// set rational to fraction of two small integers
void nct_quat_test_helper_nct_ibq_set_i(nct_ibq_t *q, int n, int d);

// test equality of integer vectors, return 0 if equal, 1 if not
int nct_quat_test_helper_nct_ibz_vec_4_equal(const nct_ibz_vec_4_t *vec, const nct_ibz_vec_4_t *cmp);

// compares integer to small int via nct_ibz_cmp
int nct_quat_test_helper_nct_ibz_equal_i(const nct_ibz_t *x, int cmp);

//return 0 if equal, 1 otherwise
int nct_ibz_mat_4x4_equal(const nct_ibz_mat_4x4_t *a, const nct_ibz_mat_4x4_t *b);

//return 0 if reduced enough, 1 otherwise
int nct_quat_test_helper_is_reduced(const nct_ibz_mat_4x4_t *mat, const nct_ibz_t *norm_ideal, const nct_quat_alg_t *alg, int print_flag);

//As in NIST code
void nct_ibz_content(nct_ibz_t *content, const nct_ibz_vec_4_t *v);
void nct_ibz_vec_4_set_ibz(nct_ibz_vec_4_t *vec, int64_t coord0, int64_t coord1, int64_t coord2, int64_t coord3);
void nct_ibz_vec_4_copy(nct_ibz_vec_4_t *new, const nct_ibz_vec_4_t  *vec);
void nct_ibz_vec_4_copy_ibz(nct_ibz_vec_4_t *res, const nct_ibz_t *coord0,const nct_ibz_t *coord1,const nct_ibz_t *coord2,const nct_ibz_t *coord3);
void nct_ibz_vec_4_negate(nct_ibz_vec_4_t *neg, const nct_ibz_vec_4_t  *vec);
void nct_ibz_vec_4_add(nct_ibz_vec_4_t *res, const nct_ibz_vec_4_t *a, const nct_ibz_vec_4_t *b);
void nct_ibz_vec_4_sub(nct_ibz_vec_4_t *res, const nct_ibz_vec_4_t *a, const nct_ibz_vec_4_t *b);
int nct_ibz_vec_4_is_zero(const nct_ibz_vec_4_t *x);
void nct_ibz_vec_4_scalar_mul(nct_ibz_vec_4_t *prod, const nct_ibz_t *scalar, const nct_ibz_vec_4_t *vec);
void nct_ibz_mat_4x4_identity(nct_ibz_mat_4x4_t *id);
int nct_ibz_mat_4x4_is_identity(const nct_ibz_mat_4x4_t *mat);
void nct_ibz_mat_4x4_zero(nct_ibz_mat_4x4_t *zero);
void nct_ibz_mat_4x4_scalar_mul(nct_ibz_mat_4x4_t *prod, const nct_ibz_t *scalar, const nct_ibz_mat_4x4_t *mat);
void nct_ibz_mat_4x4_eval(nct_ibz_vec_4_t  *res, const nct_ibz_mat_4x4_t *mat, const nct_ibz_vec_4_t *vec);
void nct_ibz_mat_4x4_gcd(nct_ibz_t *gcd, const nct_ibz_mat_4x4_t *mat);
int nct_ibz_mat_4x4_triangular_equal(const nct_ibz_mat_4x4_t *mat1, const nct_ibz_mat_4x4_t *mat2);
void nct_ibz_mat_4x4_triangular_scalar_mul(nct_ibz_mat_4x4_t *prod, const nct_ibz_t *scalar, const nct_ibz_mat_4x4_t *mat);
void nct_ibz_mat_4x4_triangular_eval(nct_ibz_vec_4_t  *res, const nct_ibz_mat_4x4_t *mat, const nct_ibz_vec_4_t *vec);//internal helper functions
void nct_ibz_mat_4x4_mul(nct_ibz_mat_4x4_t *res, const nct_ibz_mat_4x4_t *a, const nct_ibz_mat_4x4_t *b);
int nct_ibz_mat_4x4_inv_with_det_as_denom(nct_ibz_mat_4x4_t *inv, nct_ibz_t *det, const nct_ibz_mat_4x4_t *mat);

#endif