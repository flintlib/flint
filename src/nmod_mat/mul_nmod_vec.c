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

void nmod_mat_mul_nmod_vec(
    ulong * c,
    const nmod_mat_t A,
    const ulong * b, slong blen)
{
    const slong len = FLINT_MIN(A->c, blen);
    const dot_params_t params = _nmod_vec_dot_params(len, A->mod);

    for (slong i = 0; i < A->r; i++)
        c[i] = _nmod_vec_dot(A->rows[i], b, len, A->mod, params);
}

void nmod_mat_mul_nmod_vec_ptr(
    ulong * const * c,
    const nmod_mat_t A,
    const ulong * const * b, slong blen)
{
    slong i;
    slong len = FLINT_MIN(A->c, blen);
    slong nrows = A->r;
    ulong * bb, * cc;
    TMP_INIT;

    TMP_START;

    bb = TMP_ARRAY_ALLOC(len, ulong);
    cc = TMP_ARRAY_ALLOC(nrows, ulong);

    for (i = 0; i < len; i++)
        bb[i] = b[i][0];

    nmod_mat_mul_nmod_vec(cc, A, bb, len);

    for (i = 0; i < nrows; i++)
        c[i][0] = cc[i];

    TMP_END;
}
