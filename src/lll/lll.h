#ifndef LLL_H
#define LLL_H

#include "../ct_intbig/ct_intbig.h"
#include "lll_constants.h"

/// @brief LLL reduction on 4-dimensional lattice
/// Implements Algorithm 2.6.3 from Henri Cohen's "A Course in Computational Algebraic Number Theory"
/// @param red 
/// @param lattice 
/// @return 
int quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const ibz_t *q, int precision);

#endif