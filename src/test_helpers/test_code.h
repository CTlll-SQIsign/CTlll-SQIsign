#include "test_helpers.h"
#include "../bkz/bkz.h"
#include "../lll/lll.h"

void ct_to_nct_integer_translation(nct_ibz_t *nct, const ibz_t *ct);
int nct_quat_lattice_check_all(nct_quat_lattice_t *red, const nct_quat_lattice_t *inpt, const nct_ibz_t *norm, const nct_quat_alg_t *alg, int print_flag);
int translate_and_test_reduction(const quat_lattice_t *red, const quat_lattice_t *inpt, const quat_alg_t *alg, int print_flag);
int quat_test_bkz_on_lattice(const quat_lattice_t *lat ,const quat_alg_t *alg, int bkz_iterations, int lagrange_iterations, int print_flag);
int quat_test_lll_on_lattice(const quat_lattice_t *lat,const quat_alg_t *alg);
