/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"

void nmod_mat_mul_nmod_vec(
    mp_limb_t * c,
    const nmod_mat_t A,
    const mp_limb_t * b, slong blen)
{
    nmod_t mod = A->mod;
    slong i, j;
    slong len = FLINT_MIN(A->c, blen);
    int nlimbs = _nmod_vec_dot_bound_limbs(len, mod);

    for (i = A->r - 1; i >= 0; i--)
    {
        const mp_limb_t * Ai = A->rows[i];
        NMOD_VEC_DOT(c[i], j, len, Ai[j], b[j], mod, nlimbs);
    }
}

