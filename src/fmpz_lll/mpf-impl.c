/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "mpf-impl.h"
#include "fmpz.h"

mpf * _mpf_vec_init(slong len, flint_bitcnt_t prec)
{
    slong i;

    mpf * vec = flint_malloc(len * sizeof(mpf));

    for (i = 0; i < len; i++)
        mpf_init2(vec + i, prec);

    return vec;
}

void _mpf_vec_clear(mpf * vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        mpf_clear(vec + i);
    flint_free(vec);
}

void _mpf_vec_set_fmpz_vec(mpf * appv, const fmpz * vec, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        fmpz_get_mpf(appv + i, vec + i);
}

void _mpf_vec_norm(mpf_t res, const mpf * vec, slong len)
{
    slong i;
    mpf_t tmp;
    mpf_init(tmp);

    flint_mpf_set_ui(res, 0);
    for (i = 0; i < len; i++)
    {
        mpf_mul(tmp, vec + i, vec + i);
        mpf_add(res, res, tmp);
    }

    mpf_clear(tmp);
}

void _mpf_vec_norm2(mpf_t res, const mpf * vec, slong len, flint_bitcnt_t prec)
{
    slong i;
    mpf_t tmp;
    mpf_init2(tmp, prec);

    flint_mpf_set_ui(res, 0);
    for (i = 0; i < len; i++)
    {
        mpf_mul(tmp, vec + i, vec + i);
        mpf_add(res, res, tmp);
    }

    mpf_clear(tmp);
}

int
_mpf_vec_dot2(mpf_t res, const mpf * vec1, const mpf * vec2, slong len2, flint_bitcnt_t prec)
{
    slong i;
    int r = 0;
    mpf_t tmp, tmp2;
    mpf_init2(tmp, prec);
    mpf_init2(tmp2, prec);

    flint_mpf_set_ui(res, 0);
    for (i = 0; i < len2; i++)
    {
        mpf_mul(tmp, vec1 + i, vec2 + i);
        mpf_add(res, res, tmp);
    }

    _mpf_vec_norm(tmp, vec1, len2);
    _mpf_vec_norm(tmp2, vec2, len2);
    mpf_mul(tmp, tmp, tmp2);
    mpf_div_2exp(tmp, tmp, prec);
    mpf_mul(tmp2, res, res);

    if (mpf_cmp(tmp2, tmp) > 0)
        r = 1;

    mpf_clear(tmp);
    mpf_clear(tmp2);

    return r;
}

void
mpf_mat_init(mpf_mat_t mat, slong rows, slong cols, flint_bitcnt_t prec)
{
    if (rows != 0 && cols != 0)
    {
        slong i;
        mat->entries = flint_malloc(flint_mul_sizes(rows, cols) * sizeof(mpf));
        mat->rows = flint_malloc(rows * sizeof(mpf *));

        for (i = 0; i < rows * cols; i++)
            mpf_init2(mat->entries + i, prec);
        for (i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    else
    {
       mat->entries = NULL;
       mat->rows = NULL;
    }

    mat->r = rows;
    mat->c = cols;
    mat->prec = prec;
}

void mpf_mat_clear(mpf_mat_t mat)
{
    if (mat->entries)
    {
        slong i;
        for (i = 0; i < mat->r * mat->c; i++)
            mpf_clear(mat->entries + i);
        flint_free(mat->entries);
        flint_free(mat->rows);
    }
}
