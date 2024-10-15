#ifndef LLL_H
#define LLL_H

#include "../nct_intbig/nct_intbig.h"



/** @brief Type for matrices of 4x4 big integers
 * 
 * @typedef ibz_mat_4x4_t
*/
typedef nct_ibz_t nct_ibz_mat_4x4_t[4][4];

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


static inline void nct_ibz_mat_4x4_init(nct_ibz_mat_4x4_t *mat){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_init(&(*mat)[i][j]);
        }
    }
}
static inline void nct_ibz_mat_4x4_finalize(nct_ibz_mat_4x4_t *mat){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            nct_ibz_finalize(&(*mat)[i][j]);
        }
    }
}

static inline void nct_ibz_mat_4x4_copy(nct_ibz_mat_4x4_t *copy, const nct_ibz_mat_4x4_t *copied){
    for(int i = 0; i < 4; i ++){
        for(int j = 0; j < 4; j ++){
            nct_ibz_copy(&((*copy)[i][j]),&((*copied)[i][j]));
        }
    }
}

static inline void nct_ibz_mat_4x4_transpose(nct_ibz_mat_4x4_t *transposed, const nct_ibz_mat_4x4_t *mat){
    nct_ibz_mat_4x4_t work;
    nct_ibz_mat_4x4_init(&work);
    for(int i = 0; i < 4; i ++){
        for(int j = 0; j < 4; j ++){
            nct_ibz_copy(&(work[i][j]),&((*mat)[j][i]));
        }
    }
    nct_ibz_mat_4x4_copy(transposed,&work);
    nct_ibz_mat_4x4_finalize(&work);
}



/// @brief LLL reduction on 4-dimensional lattice
/// Implements Algorithm 2.6.3 from Henri Cohen's "A Course in Computational Algebraic Number Theory"
/// @param red 
/// @param lattice 
/// @return 
int nct_quat_lattice_lll(nct_ibz_mat_4x4_t *red, const nct_quat_lattice_t *lattice, const nct_ibz_t *q, int precision);

#endif