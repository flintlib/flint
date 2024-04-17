/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "mpn_mod.h"

/* todo: retune this when arithmetic is more optimized */

/* cutoffs for classical -> waksman and waksman -> multi_mod */
/* equal cutoffs means waksman is not used */
static const short mat_mul_cutoff_tab[][2] = {
  {88, 88},   /* bits = 80 */
  {230, 230},   /* bits = 96 */
  {220, 220},   /* bits = 112 */
  {160, 160},   /* bits = 128 */
  {10, 39},   /* bits = 144 */
  {9, 57},   /* bits = 160 */
  {17, 57},   /* bits = 176 */
  {74, 74},   /* bits = 192 */
  {8, 64},   /* bits = 208 */
  {8, 64},   /* bits = 224 */
  {8, 92},   /* bits = 240 */
  {17, 63},   /* bits = 256 */
  {8, 69},   /* bits = 272 */
  {8, 65},   /* bits = 288 */
  {8, 91},   /* bits = 304 */
  {14, 67},   /* bits = 320 */
  {8, 75},   /* bits = 336 */
  {7, 77},   /* bits = 352 */
  {7, 93},   /* bits = 368 */
  {10, 79},   /* bits = 384 */
  {6, 77},   /* bits = 400 */
  {6, 91},   /* bits = 416 */
  {6, 93},   /* bits = 432 */
  {10, 73},   /* bits = 448 */
  {6, 75},   /* bits = 464 */
  {6, 87},   /* bits = 480 */
  {6, 87},   /* bits = 496 */
  {9, 71},   /* bits = 512 */
  {6, 70},   /* bits = 528 */
  {5, 77},   /* bits = 544 */
  {5, 84},   /* bits = 560 */
  {6, 67},   /* bits = 576 */
  {5, 75},   /* bits = 592 */
  {5, 74},   /* bits = 608 */
  {4, 78},   /* bits = 624 */
  {6, 63},   /* bits = 640 */
  {4, 68},   /* bits = 656 */
  {4, 68},   /* bits = 672 */
  {4, 75},   /* bits = 688 */
  {6, 59},   /* bits = 704 */
  {4, 66},   /* bits = 720 */
  {4, 66},   /* bits = 736 */
  {4, 71},   /* bits = 752 */
  {5, 59},   /* bits = 768 */
  {4, 59},   /* bits = 784 */
  {4, 66},   /* bits = 800 */
  {4, 63},   /* bits = 816 */
  {4, 58},   /* bits = 832 */
  {4, 57},   /* bits = 848 */
  {4, 61},   /* bits = 864 */
  {4, 61},   /* bits = 880 */
  {4, 55},   /* bits = 896 */
  {4, 57},   /* bits = 912 */
  {4, 64},   /* bits = 928 */
  {4, 65},   /* bits = 944 */
  {4, 57},   /* bits = 960 */
  {4, 61},   /* bits = 976 */
  {4, 62},   /* bits = 992 */
  {4, 65},   /* bits = 1008 */
  {4, 53},   /* bits = 1024 */
};

int
mpn_mod_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong ar = A->r;

    if (ar <= 3)
        return gr_mat_mul_classical(C, A, B, ctx);

    slong ac = A->c;
    slong bc = B->c;

    slong tab_i, cutoff_waksman, cutoff_multi_mod;
    slong bits = MPN_MOD_CTX_MODULUS_BITS(ctx);

    tab_i = (bits - FLINT_BITS - 1) / 16;

    cutoff_waksman = mat_mul_cutoff_tab[tab_i][0];
    cutoff_multi_mod = mat_mul_cutoff_tab[tab_i][1];

    if (ar < cutoff_waksman || ac < cutoff_waksman || bc < cutoff_waksman)
        return gr_mat_mul_classical(C, A, B, ctx);

    if (ar < cutoff_multi_mod || ac < cutoff_multi_mod || bc < cutoff_multi_mod)
        return mpn_mod_mat_mul_waksman(C, A, B, ctx);

    /* special case: near the two-limb boundary, strassen beats multi_mod on a single thread */
#ifndef FLINT_USES_BLAS
    if (bits >= 113 && bits <= 128 && flint_get_num_available_threads() == 1)
        return gr_mat_mul_strassen(C, A, B, ctx);
#endif

    return mpn_mod_mat_mul_multi_mod(C, A, B, ctx);
}
