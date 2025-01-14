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
    ulong * c,
    const ulong * a, slong alen,
    const nmod_mat_t B)
{
    slong i;
    slong len = FLINT_MIN(B->r, alen);
    slong ncols = B->c;

    /* scalar_addmul wants non-empty */
    if (ncols < 1)
        return;

    if (len > 0)
        _nmod_vec_scalar_mul_nmod(c, nmod_mat_entry_ptr(B, 0, 0), ncols, a[0], B->mod);
    else
        _nmod_vec_zero(c, ncols);

    for (i = 1; i < len; i++)
        _nmod_vec_scalar_addmul_nmod(c, nmod_mat_entry_ptr(B, i, 0), ncols, a[i], B->mod);
}

void nmod_mat_nmod_vec_mul_ptr(
    ulong * const * c,
    const ulong * const * a, slong alen,
    const nmod_mat_t B)
{
    slong i;
    slong len = FLINT_MIN(B->r, alen);
    slong ncols = B->c;
    ulong * aa, * cc;
    TMP_INIT;

    TMP_START;

    aa = TMP_ARRAY_ALLOC(len, ulong);
    cc = TMP_ARRAY_ALLOC(ncols, ulong);

    for (i = 0; i < len; i++)
        aa[i] = a[i][0];

    nmod_mat_nmod_vec_mul(cc, aa, len, B);

    for (i = 0; i < ncols; i++)
        c[i][0] = cc[i];

    TMP_END;
}
