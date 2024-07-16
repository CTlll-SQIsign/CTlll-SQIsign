/** @file
 * 
 * @authors Luca De Feo, Sina Schaeffler
 * 
 * @brief Declarations for big integers in the reference implementation
*/

#ifndef NCT_INTBIG_H
#define NCT_INTBIG_H

#include <gmp.h>
#include <stdint.h>
//#include <tutil.h>

/** @defgroup nct_ibz_all Signed big integers
 * @{
*/

/** @defgroup nct_ibz_t Precise number types
 * @{
*/

/** @brief Type for signed long integers
 * 
 * @typedef nct_ibz_t
 * 
 * For integers of arbitrary size, used by intbig module, using gmp
*/
typedef mpz_t nct_ibz_t ;

/** @brief Type for fractions of integers
 * 
 * @typedef nct_ibq_t
 * 
 * For fractions of integers of arbitrary size, used by intbig module, using gmp
*/
typedef mpq_t nct_ibq_t ;


/** @brief Type for vector of 2 big integers
 * 
 * @typedef nct_ibz_vec_2_t
*/
typedef nct_ibz_t nct_ibz_vec_2_t[2];

/** @}
*/


/********************************************************************/

/** @defgroup nct_ibz_c Constants
 * @{
*/

/**
 * Constant zero
*/
extern const nct_ibz_t nct_ibz_const_zero;

/**
 * Constant one
*/
extern const nct_ibz_t nct_ibz_const_one;

/**
 * Constant two
*/
extern const nct_ibz_t nct_ibz_const_two;

/**
 * Constant three
*/
extern const nct_ibz_t nct_ibz_const_three;

/** @}
*/

/* constructors/destructors (possibly no-ops) */

/** @defgroup nct_ibz_c Constructors and Destructors
 * @{
*/

void nct_ibz_init(nct_ibz_t *x);
void nct_ibq_init(nct_ibq_t *x);
void nct_ibz_vec_2_init(nct_ibz_vec_2_t *vec);


void nct_ibz_finalize(nct_ibz_t *x);
void nct_ibq_finalize(nct_ibq_t *x);
void nct_ibz_vec_2_finalize(nct_ibz_vec_2_t *vec);

/** @brief overwrites memory before freeing it
*/
void nct_ibz_secure_finalize(nct_ibz_t *x);

/** @brief overwrites memory before freeing it
*/
void nct_ibq_secure_finalize(nct_ibq_t *x);

/** @brief overwrites memory before freeing it
*/
void nct_ibz_vec_2_secure_finalize(nct_ibz_vec_2_t *vec);

/** @}
*/

/** @defgroup nct_ibz_za Basic integer arithmetic
 * @{
*/

/** @brief sum=a+b
*/
void nct_ibz_add(nct_ibz_t *sum, const nct_ibz_t *a, const nct_ibz_t *b);

/** @brief diff=a-b
*/
void nct_ibz_sub(nct_ibz_t *diff, const nct_ibz_t *a, const nct_ibz_t *b);

/** @brief prod=a*b
*/
void nct_ibz_mul(nct_ibz_t *prod, const nct_ibz_t *a, const nct_ibz_t *b);

/** @brief prod=a*2^exp
*/
void nct_ibz_mul_2exp(nct_ibz_t *prod, const nct_ibz_t *a, uint64_t exp);

/** @brief neg=-a
*/
void nct_ibz_neg(nct_ibz_t *neg, const nct_ibz_t* a);

/** @brief abs=|a|
*/
void nct_ibz_abs(nct_ibz_t *abs, const nct_ibz_t* a);

/** @brief Euclidean division of a by b
 * 
 * Computes quotient, remainder so that remainder+quotient*b = a where 0<=|remainder|<|b|
 * The quotient is rounded towards zero.
*/
void nct_ibz_div(nct_ibz_t *quotient, nct_ibz_t *remainder, const nct_ibz_t *a, const nct_ibz_t *b);

/** @brief Euclidean division of a by b
 *
 * Computes quotient, remainder so that remainder+quotient*b = a where 0<=|remainder|<|b|
 * The quotient is rounded towards minus infinity.
 */
void nct_ibz_div_floor(nct_ibz_t *q, nct_ibz_t *r, const nct_ibz_t *n, const nct_ibz_t *d);

/** @brief Euclidean division of a by 2^exp
 * 
 * Computes a right shift of abs(a) by exp bits, then sets sign(quotient) to sign(a).
 * 
 * Division and rounding is as in nct_ibz_div.
*/
void nct_ibz_div_2exp(nct_ibz_t *quotient, const nct_ibz_t *a, uint64_t exp);

/** @brief r = a mod b
 * 
 * Assumes valid inputs
 * The sign of the divisor is ignored, the result is always non-negative
*/
void nct_ibz_mod(nct_ibz_t *r, const nct_ibz_t *a, const nct_ibz_t *b);

/** @brief Test if a = 0 mod b
*/
int nct_ibz_divides(const nct_ibz_t *a, const nct_ibz_t *b);

/** @brief pow=x^e
 * 
 * Assumes valid inputs, The case 0^0 yields 1.
*/
void nct_ibz_pow(nct_ibz_t *pow, const nct_ibz_t *x, int64_t e);

/** @brief pow=(x^e) mod m
 * 
 * Assumes valid inputs
*/
void nct_ibz_pow_mod(nct_ibz_t *pow, const nct_ibz_t *x, const nct_ibz_t *e, const nct_ibz_t *m);

/** @brief Compare a and b
 * 
 * @returns a positive value if a > b, zero if a = b, and a negative value if a < b
*/
int nct_ibz_cmp(const nct_ibz_t *a, const nct_ibz_t *b);

/** @brief Test if x is 0
 * 
 * @returns 1 if x=0, 0 otherwise
*/
int nct_ibz_is_zero(const nct_ibz_t *x);

/** @brief Test if x is 1
 * 
 * @returns 1 if x=1, 0 otherwise
*/
int nct_ibz_is_one(const nct_ibz_t *x);

/** @brief Test if x is even
 * 
 * @returns 1 if x is even, 0 otherwise
*/
int nct_ibz_is_even(const nct_ibz_t *x);

/** @brief set i to value x
*/
void nct_ibz_set(nct_ibz_t *i, int64_t x);

/** @brief set i to integer ontained in string when read as number in base
 * 
 * Base should be 10 or 16, and the number should be written without ponctuation or whitespaces
 * 
 * Case for base 16 does not matter
 * 
 * @returns  1 if the string could be converted to an integer, 0 otherwise
 */
int nct_ibz_set_from_str(nct_ibz_t *i, const char *str, int base);

/** @brief Copy value into target
*/
void nct_ibz_copy(nct_ibz_t *target, const nct_ibz_t *value);

/** @brief Exchange values of a and b
 */
void nct_ibz_swap(nct_ibz_t *a, nct_ibz_t *b);

/** @brief Copy dig array to target, given digits and the length of the dig array
 * 
 *  @param target Target nct_ibz_t element
 *  @param dig array of digits
 *  @param dig_len length of the digits array
*/
//void nct_ibz_copy_digits(nct_ibz_t *target, const digit_t *dig, int dig_len);
//#define nct_ibz_copy_digit_array(I, T) do { nct_ibz_copy_digits((I), (T), sizeof(T)/sizeof(*(T))); } while (0)

/** @brief Copy an nct_ibz_t to target digit_t array.
 *  Restrictions: ibz >= 0 and target must hold sufficient elements to hold ibz 
 * 
 *  @param target Target digit_t array
 *  @param ibz ibz source nct_ibz_t element
*/
//void nct_ibz_to_digits(digit_t *target, const nct_ibz_t *ibz);
//#define nct_ibz_to_digit_array(T, I) do { memset((T), 0, sizeof(T)); nct_ibz_to_digits((T), (I)); } while (0)

/** @brief get int64_t equal to the lowest bits of i
*/
int64_t nct_ibz_get(const nct_ibz_t *i);

//void nct_ibz_printf(const char* format, ...);
#define nct_ibz_printf gmp_printf

/** @brief generate random value in [a, b]
 *  assumed that a >= 0 and b >= 0 and a < b
 * @returns 1 on success, 0 on failiure
*/
//int nct_ibz_rand_interval(nct_ibz_t *rand, const nct_ibz_t *a, const nct_ibz_t *b);

/** @brief generate random value in [a, b]
 *  assumed that a >= 0, b >= 0 and a < b
 * @returns 1 on success, 0 on failiure
*/
//int nct_ibz_rand_interval_i(nct_ibz_t *rand, int64_t a, int64_t b);

/** @brief generate random value in [-m, m]
 *  assumed that m > 0 and bitlength of m < 64 bit
 * @returns 1 on success, 0 on failiure
*/
//int nct_ibz_rand_interval_minm_m(nct_ibz_t *rand, int64_t m);


/** @brief Bitsize of a.
 * 
 *  @returns Bitsize of a.
 * 
*/
int nct_ibz_bitsize(const nct_ibz_t *a);


/* etc....*/

/** @}
*/

/** @defgroup nct_ibz_qa Basic fraction arithmetic
 * @{
*/

/** @brief sum=a+b
*/
void nct_ibq_add(nct_ibq_t *sum, const nct_ibq_t *a, const nct_ibq_t *b);

/** @brief diff=a-b
*/
void nct_ibq_sub(nct_ibq_t *diff, const nct_ibq_t *a, const nct_ibq_t *b);

/** @brief neg=-x
 */
void nct_ibq_neg(nct_ibq_t *neg, const nct_ibq_t *x);

/** @brief abs=|x|
 */
void nct_ibq_abs(nct_ibq_t *abs, const nct_ibq_t *x);

/** @brief prod=a*b
*/
void nct_ibq_mul(nct_ibq_t* prod, const nct_ibq_t *a, const nct_ibq_t *b);

/** @brief inv=1/x
 * 
 * @returns 0 if x is 0, 1 if inverse exists and was computed
*/
int nct_ibq_inv(nct_ibq_t *inv, const nct_ibq_t *x);

/** @brief quot = a/b
 * @param quot Output a/b
 * @param a
 * @param b must not be 0
*/
void nct_ibq_div(nct_ibq_t *quot, const nct_ibq_t *a, const nct_ibq_t *b);

/** @brief Compare a and b
 * 
 * @returns a positive value if a > b, zero if a = b, and a negative value if a < b
*/
int nct_ibq_cmp(const nct_ibq_t *a, const nct_ibq_t *b);

/** @brief Test if x is 0
 * 
 * @returns 1 if x=0, 0 otherwise
*/
int nct_ibq_is_zero(const nct_ibq_t *x);

/** @brief Test if x is 1
 *
 * @returns 1 if x=1, 0 otherwise
*/
int nct_ibq_is_one(const nct_ibq_t *x);

/** @brief Set q to a/b if b not 0
 *
 * @returns 1 if b not 0 and q is set, 0 otherwise
*/
int nct_ibq_set(nct_ibq_t *q, const nct_ibz_t *a, const nct_ibz_t *b);

/** @brief Set x to 0
 * 
 * Assumes x is initialized
*/
void nct_ibq_set_zero(nct_ibq_t *x);

/** @brief Set x to 1
 * 
 * Assumes x is initialized
*/
void nct_ibq_set_one(nct_ibq_t *x);

/** @brief Copy value into target
*/
void nct_ibq_copy(nct_ibq_t *target, const nct_ibq_t *value);

/** @brief Exchange values of a and b
 */
void nct_ibq_swap(nct_ibq_t *a, nct_ibq_t *b);

/** @brief Denominator of x
*/
void nct_ibq_denom(nct_ibz_t *d, const nct_ibq_t *x);

/** @brief Numerator of x
*/
void nct_ibq_num(nct_ibz_t *n, const nct_ibq_t *x);

/** @brief Checks if q is an integer
 *  
 * @returns 1 if yes, 0 if not
*/
int nct_ibq_is_ibz(const nct_ibq_t *q);

/**
 * @brief Converts a fraction q to an integer y, if q is an integer.
 * 
 * @returns 1 if z is an integer, 0 if not
*/
int nct_ibq_to_ibz(nct_ibz_t *z, const nct_ibq_t *q);

/** @}
*/

/** @defgroup nct_ibz_n Number theory functions
 * @{
*/


/**
 * @brief Probabilistic primality test
 *
 * @param n The number to test
 * @param reps Number of Miller-Rabin repetitions. The more, the slower and the less likely are false positives
 * @return 1 if probably prime, 0 if certainly not prime, 2 if certainly prime
 * 
 * Using GMP's implementation:
 * 
 * From GMP's documentation: "This function performs some trial divisions, a Baillie-PSW probable prime test, then reps-24 Miller-Rabin probabilistic primality tests."
 */
int nct_ibz_probab_prime(const nct_ibz_t *n, int reps);

/**
 * @brief Greatest common divisor
 *
 * @param gcd Output: Set to the gcd of a and b
 * @param a
 * @param b
 */
void nct_ibz_gcd(nct_ibz_t *gcd, const nct_ibz_t *a, const nct_ibz_t *b);

/**
 * @brief GCD and Bézout coefficients u, v such that ua + bv = gcd
 *
 * @param gcd Output: Set to the gcd of a and b
 * @param u Output: integer such that ua+bv=gcd
 * @param v Output: Integer such that ua+bv=gcd
 * @param a
 * @param b
 */
void nct_ibz_xgcd(nct_ibz_t *gcd, nct_ibz_t *u, nct_ibz_t *v, const nct_ibz_t *a, const nct_ibz_t *b);

/**
 * @brief GCD, Bézout coefficients u, v such that ua + bv = gcd, and annihilators s, t such that sa + bt = 0
 *
 * @param gcd Output: Set to the gcd of a and b
 * @param s Output: integer such that sa+bt=0
 * @param t Output: Integer such that sa+bt=0
 * @param u Output: integer such that ua+bv=gcd
 * @param v Output: Integer such that ua+bv=gcd
 * @param a
 * @param b
 */
void nct_ibz_xgcd_ann(nct_ibz_t *gcd, nct_ibz_t *s, nct_ibz_t *t, nct_ibz_t *u, nct_ibz_t *v, const nct_ibz_t *a, const nct_ibz_t *b);


/**
 * @brief Modular inverse
 *
 * @param inv Output: Set to the integer in [0,mod[ such that a*inv = 1 mod (mod) if it exists
 * @param a
 * @param mod
 * @returns 1 if inverse exists and was computed, 0 otherwise
 */
int nct_ibz_invmod(nct_ibz_t *inv, const nct_ibz_t *a, const nct_ibz_t *mod);

/**
 * @brief Chinese remainders
 *
 * @param crt Output: Set so that crt = a mod (mod_a), crt = b mod (mod_b)
 * @param a
 * @param b
 * @param mod_a
 * @param mod_b
 */
void nct_ibz_crt(nct_ibz_t *crt, const nct_ibz_t *a, const nct_ibz_t *b, const nct_ibz_t *mod_a, const nct_ibz_t *mod_b);

/**
 * @brief Kronecker symbol of a mod b
 *
 * @returns Kronecker symbol of a mod b
 * @param a
 * @param b
 * 
 * Uses GMP's implementation
 */
int nct_ibz_kronecker(const nct_ibz_t *a, const nct_ibz_t *b);

/**
 * @brief Jacobi symbol of a mod odd
 *
 * @returns jacobi symbol of a mod odd
 * @param a
 * @param odd assumed odd
 * 
 * Uses GMP's implementation
 * 
 * If output is -1, a is a not a square mod odd
 */
int nct_ibz_jacobi(const nct_ibz_t *a, const nct_ibz_t *odd);

/**
 * @brief Legendre symbol of a mod p
 *
 * @returns Legendre symbol of a mod p
 * @param a
 * @param p assumed prime
 * 
 * Uses GMP's implementation
 * 
 * If output is 1, a is a square mod p, if -1, not. If 0, it is divisible by p
 */
int nct_ibz_legendre(const nct_ibz_t *a, const nct_ibz_t *p);



/**
 * @brief Integer square root of a perfect square
 *
 * @returns 1 if an integer square root of a exists and was computed, 0 otherwise
 * @param sqrt Output: Set to a integer square root of a if any exist
 * @param a number of which an integer square root is searched
 */
int nct_ibz_sqrt(nct_ibz_t *sqrt, const nct_ibz_t *a);

/**
 * @brief Floor of Integer square root
 *
 * @param sqrt Output: Set to the floor of an integer square root
 * @param a number of which a floor of an integer square root is searched
 */
void nct_ibz_sqrt_floor(nct_ibz_t *sqrt, const nct_ibz_t *a);

/**
 * @brief Square root modulo a prime
 *
 * @returns 1 if square root of a mod p exists and was computed, 0 otherwise
 * @param sqrt Output: Set to a square root of a mod p if any exist
 * @param a number of which a square root mod p is searched
 * @param p assumed prime
 */
int nct_ibz_sqrt_mod_p(nct_ibz_t *sqrt, const nct_ibz_t *a, const nct_ibz_t *p);

/**
 * @brief Square root modulo a the double of a given prime
 *
 * @returns 1 if square root of a mod (2p) exists and was computed, 0 otherwise
 * @param sqrt Output: Set to a square root of a mod (2p) if any exist
 * @param a number of which a square root mod (2p) is searched
 * @param p assumed prime
 */
int nct_ibz_sqrt_mod_2p(nct_ibz_t *sqrt, const nct_ibz_t *a, const nct_ibz_t *p);

/** @}
*/

// end of nct_ibz_all
/** @}
*/
#endif
