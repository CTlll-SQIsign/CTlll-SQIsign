#ifndef BKZ_H
#define BKZ_H

#include "../ct_helpers/ct_helpers.h"

/** @brief BKZ data structure
 *
 * Stores varying values of the basis, the gram matrixG , mu and r during a BKZ-2 run.
 */
typedef struct quat_bkz_matrix {
    ibz_mat_4x4_t basis;
    ibz_mat_4x4_t G;
    ibq_mat_4x4_t mu;
    ibq_mat_4x4_t r;
} quat_bkz_matrix_t;

/** @brief Initializes a bkz matrix structure
 *
 * @param g The matrix to be initialized (must be uninitialized when input)
 */
void quat_bkz_matrix_init(quat_bkz_matrix_t *g);

/** @brief Finalizes a bkz matrix structure
 *
 * @param g The matrix to be finalized (must be initialized when input)
 */
void quat_bkz_matrix_finalize(quat_bkz_matrix_t *g);

/** @brief (a,b)=U(a,b): Multiply 2x4 matrix by 2x2 matrix
 *
 * @param a In-and Output: 1st vector of the 2x4 matrix
 * @param b In-and Output: 2nd vector of the 2x4 matrix
 * @param U The 2x2 matrix by which the matrix (a,b) is multiplied
 */
void ibz_vec_2x4_mul(ibz_vec_4_t *a, ibz_vec_4_t *b,const ibz_mat_2x2_t *U);

/** @brief Sets all fields of a bkz matrix for input basis mat
 *
 * @param g Output: The bkz matrix to be set
 * @param mat Matrix containing a basis of a (rank 4) lattice
 * @param alg The quaternion algebra whose norm is used for reducing
 */
void quat_bkz_matrix_set(quat_bkz_matrix_t *g, const ibz_mat_4x4_t *mat, const quat_alg_t *alg);

/** @brief Partial BKZ recumputation after basis an G are modified
 * 
 * Updates all fields which are only related to the basis elements before index end of a bkz matrix after its basis and gram matrix G was modified between indices start and end.
 *
 * @param g Output: The bkz matrix to be updated
 * @param start Index of the first possibly modified entry
 * @param end Index of the first entry outside of the subspace for which we want to update the matrix. Must be larger than start.
 */
void quat_bkz_update_matrix(quat_bkz_matrix_t *g, int start, int end);

/** @brief Size-reduces a bkz matrix and maintains all fields
 *
 * @param g Output: The bkz matrix to be size-reduced
 * @param start Index of the first entry to get size-reduced
 * @param end Index of the entry after the last entry to get size-reduced
 */
void quat_bkz_size_reduce(quat_bkz_matrix_t *g, int start, int end);

/** @brief Updates bkz matrix after lagrange
 *
 * @param g Output: The bkz matrix to be adapted
 * @param U transformation matrix output by lagrange
 * @param index Index of the first of the two neighboring vectors in the basis given by the matrix g on which lagrange was called to obtain U
 */
void quat_bkz_update_after_lagrange(quat_bkz_matrix_t *g, const ibz_mat_2x2_t *U, int index);

/** @brief Computes the transformation matrix for lagrange-reducing a matrix
 *
 * @param U Output: Ontained transformation matrix
 * @param G Gram matrix of he two vectors input to lagrange
 * @param lagrange_tours Number of loop iterations. Should be at least the logarithm of the largest imput norm to giuarantee that a shortest basis is obtained
 */
void quat_bkz_lagrange_reduction_gram(ibz_mat_2x2_t *U, ibq_mat_2x2_t *G, int lagrange_tours);

/** @brief Constant-time reduction of a rank 4 lattice in dimension 4
 * @param red Output: Reduced basis, on the same denominator than lattice
 * @param lattice Lattice input
 * @param bkz_tours Outer tours
 * @param lagrange_tours Number of tours in each call to the lagrange subroutine (called 3*bkz_tours times)
 * @param alg Quaternion algebra providing the norm for the reduction
 */
void quat_lattice_bkz(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, int bkz_tours, int lagrange_tours, const quat_alg_t *alg);

#endif