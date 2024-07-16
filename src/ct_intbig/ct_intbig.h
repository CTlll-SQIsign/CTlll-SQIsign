#ifndef INTBIG_H
#define INTBIG_H

#include <gmp.h>
#include <stdint.h>

#include "../limbnum.h"
#define STRSIZE 20*LIMBNUM
#define PRINTSIZE 220*LIMBNUM
#define GCD_CONST_0 64*LIMBNUM //d
#define GCD_CONST_1 ((49*GCD_CONST_0+57)/17) //m
#define GCD_CONST_2 GCD_CONST_1+ GCD_CONST_0// d + m

/** @defgroup ibz_all Signed big integers
 * @{
*/

/** @defgroup ibz_t Precise number types
 * @{
*/

/** @brief Type for signed long integers
 * 
 * @typedef ibz_t
 * 
 * For integers of arbitrary size, used by intbig module, using gmp
*/
typedef mp_limb_t ibz_t[LIMBNUM];

/** @brief Type for fractions of integers
 * 
 * @typedef ibq_t
 * 
 * For fractions of integers of arbitrary size, used by intbig module, using gmp
*/
typedef ibz_t ibq_t[2];


/** @brief Type for matrices of 4x4 big integers
 * 
 * @typedef ibz_mat_4x4_t
*/
typedef ibz_t ibz_mat_4x4_t[4][4];

/** @brief Type for matrices of 2x2 big integers
 * 
 * @typedef ibz_mat_2x2_t
*/
typedef ibz_t ibz_mat_2x2_t[2][2];

/** @brief Type for vectors of 4 big integers
 * 
 * @typedef ibz_vec_4_t
*/
typedef ibz_t ibz_vec_4_t[4];


/** @brief Type for matrices of 4x4 rationals
 * 
 * @typedef ibq_mat_4x4_t
*/
typedef ibq_t ibq_mat_4x4_t[4][4];

/** @brief Type for matrices of 2x2 rationals
 * 
 * @typedef ibq_mat_4x4_t
*/
typedef ibq_t ibq_mat_2x2_t[2][2];
/** @}
*/


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
typedef struct quat_lattice {
  ibz_t denom; ///< Denominator by which the basis is divided (big integer, must not be 0)
  ibz_mat_4x4_t basis; ///< Integer basis of the lattice in hermite normal form (its columns divided by denom are algebra elements in the usual basis)
} quat_lattice_t;


/********************************************************************/

/* constructors/destructors (possibly no-ops) */

/** @defgroup ibz_c Constructors and Destructors
 * @{
*/

void ibz_init(ibz_t *x);
void ibq_init(ibq_t *x);
void ibz_mat_4x4_init(ibz_mat_4x4_t *x);
void ibz_mat_2x2_init(ibz_mat_2x2_t *x);
void ibq_mat_2x2_init(ibq_mat_2x2_t *x);
void ibq_mat_4x4_init(ibq_mat_4x4_t *x);
void ibz_vec_4_init(ibz_vec_4_t *x);

void ibz_finalize(ibz_t *x);
void ibq_finalize(ibq_t *x);
void ibz_mat_4x4_finalize(ibz_mat_4x4_t *x);
void ibz_mat_2x2_finalize(ibz_mat_2x2_t *x);
void ibq_mat_2x2_finalize(ibq_mat_2x2_t *x);
void ibq_mat_4x4_finalize(ibq_mat_4x4_t *x);
void ibz_vec_4_finalize(ibz_vec_4_t *x);

/** @}
*/

/** @defgroup ibz_za Basic integer arithmetic
 * @{
*/

/** @brief sum=a+b
*/
void ibz_add(ibz_t *sum, const ibz_t *a, const ibz_t *b);

/** @brief diff=a-b
*/
void ibz_sub(ibz_t *diff, const ibz_t *a, const ibz_t *b);

/** @brief prod=a*b
*/
void ibz_mul(ibz_t *prod, const ibz_t *a, const ibz_t *b);

/** @brief if cond {a,b = b,a}
*/
void ibz_conditional_swap(ibz_t *a, ibz_t *b, int cond);

/** @brief Euclidean division of a by b
 *
 * Computes quotient, remainder so that remainder+quotient*b = a where 0<=|remainder|<|b|
 * The quotient is rounded towards minus infinity.
 */
void ibz_div_floor(ibz_t *q, ibz_t *r, const ibz_t *n, const ibz_t *d);

/** @brief r = a mod b
 * 
 * Assumes valid inputs
 * The sign of the divisor is ignored, the result is always non-negative
*/
void ibz_mod(ibz_t *r, const ibz_t *a, const ibz_t *b);


/** @brief Compare a and b
 * 
 * @returns a positive value if a > b, zero if a = b, and a negative value if a < b
*/
int ibz_cmp(const ibz_t *a, const ibz_t *b);

/** @brief Test if x is 0
 * 
 * @returns 1 if x=0, 0 otherwise
*/
int ibz_is_zero(const ibz_t *x);

/** @brief set i to value x
*/
void ibz_set(ibz_t *i, int64_t x);

int ibz_set_from_str(ibz_t *i, const char *str, int base);

/** @brief Copy value into target
*/
void ibz_copy(ibz_t *target, const ibz_t *value);

/** @brief Exchange values of a and b
 */
void ibz_swap(ibz_t *a, ibz_t *b);

/** @brief sum=a+b
*/
void ibq_add(ibq_t *sum, const ibq_t *a, const ibq_t *b);

/** @brief diff=a-b
*/
void ibq_sub(ibq_t *diff, const ibq_t *a, const ibq_t *b);

/** @brief neg=-x
 */
void ibq_neg(ibq_t *neg, const ibq_t *x);

/** @brief abs=|x|
 */
void ibq_abs(ibq_t *abs, const ibq_t *x);

/** @brief prod=a*b
*/
void ibq_mul(ibq_t* prod, const ibq_t *a, const ibq_t *b);

/** @brief quot = a/b
 * @param quot Output a/b
 * @param a
 * @param b must not be 0
*/
void ibq_div(ibq_t *quot, const ibq_t *a, const ibq_t *b);

/** @brief Compare a and b
 * 
 * @returns a positive value if a > b, zero if a = b, and a negative value if a < b
*/
int ibq_cmp(const ibq_t *a, const ibq_t *b);

/** @brief Test if x is 0
 * 
 * @returns 1 if x=0, 0 otherwise
*/
int ibq_is_zero(const ibq_t *x);

/** @brief Set q to a/b if b not 0
 *
 * @returns 1 if b not 0 and q is set, 0 otherwise
*/
int ibq_set(ibq_t *q, const ibz_t *a, const ibz_t *b);

/** @brief Copy value into target
*/
void ibq_copy(ibq_t *target, const ibq_t *value);

/** @brief Exchange values of a and b
 */
void ibq_swap(ibq_t *a, ibq_t *b);

/** @brief Denominator of x
*/
void ibq_denom(ibz_t *d, const ibq_t *x);

/** @brief Numerator of x
*/
void ibq_num(ibz_t *n, const ibq_t *x);

/**
 * @brief Converts a fraction q to an integer y, if q is an integer.
 * 
 * @returns 1 if z is an integer, 0 if not
*/
int ibq_to_ibz(ibz_t *z, const ibq_t *q);

/** @}
*/

/** @defgroup ibz_n Number theory functions
 * @{
*/

/**
 * @brief Greatest common divisor
 *
 * @param gcd Output: Set to the gcd of a and b
 * @param a
 * @param b
 */
void ibz_gcd(ibz_t *gcd, const ibz_t *a, const ibz_t *b);
/** @}
*/

void ibz_mat_4x4_transpose(ibz_mat_4x4_t *transposed, const ibz_mat_4x4_t *mat);

void ibz_mat_4x4_copy(ibz_mat_4x4_t *copy, const ibz_mat_4x4_t *copied);
// end of ibz_all
/** @}
*/

/** @defgroup ibz_help Internal helpers
 * @{
*/
void ibz_to_str(char *result, const ibz_t *z);
void ibz_printf_1(char *format, const ibz_t *z);
void ibq_printf_1(char *format, const ibq_t *q);
void ibz_print_limbs(const ibz_t *z);
void ibz_neg(ibz_t *n, const ibz_t *x);
void ibz_abs(ibz_t *a, const ibz_t *x);
void ibq_red(ibq_t *red);
void ibz_bitwise_and(ibz_t *res, const ibz_t *a, const ibz_t *b);
void ibz_internal_truncate(ibz_t *f, int t);
void ibz_internal_divsteps(ibz_t *a, ibz_t *b);
/** @}
*/
#endif

