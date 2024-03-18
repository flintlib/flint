/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_mat.h"

/* TODO try delaying the reductions */
void nmod_mat_nmod_vec_mul(
    mp_limb_t * c,
    const mp_limb_t * a, slong alen,
    const nmod_mat_t B)
{
    slong i;
    slong len = FLINT_MIN(B->r, alen);
    slong ncols = B->c;

    /* scalar_addmul wants non-empty */
    if (ncols < 1)
        return;

    if (len > 0)
        _nmod_vec_scalar_mul_nmod(c, B->rows[0], ncols, a[0], B->mod);
    else
        _nmod_vec_zero(c, ncols);

    for (i = 1; i < len; i++)
        _nmod_vec_scalar_addmul_nmod(c, B->rows[i], ncols, a[i], B->mod);
}

void nmod_mat_nmod_vec_mul_ptr(
    mp_limb_t * const * c,
    const mp_limb_t * const * a, slong alen,
    const nmod_mat_t B)
{
    slong i;
    slong len = FLINT_MIN(B->r, alen);
    slong ncols = B->c;
    mp_limb_t * aa, * cc;
    TMP_INIT;

    TMP_START;

    aa = TMP_ARRAY_ALLOC(len, mp_limb_t);
    cc = TMP_ARRAY_ALLOC(ncols, mp_limb_t);

    for (i = 0; i < len; i++)
        aa[i] = a[i][0];

    nmod_mat_nmod_vec_mul(cc, aa, len, B);

    for (i = 0; i < ncols; i++)
        c[i][0] = cc[i];

    TMP_END;
}
