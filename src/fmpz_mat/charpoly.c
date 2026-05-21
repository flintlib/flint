/*
    Copyright (C) 2012, 2013 Sebastian Pancratz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr.h"
#include "gr_mat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_poly.h"

void _fmpz_mat_charpoly_berkowitz(fmpz *cp, const fmpz_mat_t mat)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz(ctx);
    GR_MUST_SUCCEED(_gr_mat_charpoly_berkowitz(cp, (const gr_mat_struct *) mat, ctx));
}

void fmpz_mat_charpoly_berkowitz(fmpz_poly_t cp, const fmpz_mat_t mat)
{
     fmpz_poly_fit_length(cp, mat->r + 1);
    _fmpz_poly_set_length(cp, mat->r + 1);
    _fmpz_mat_charpoly_berkowitz(cp->coeffs, mat);
}

#define MAT(ii, jj) fmpz_mat_entry(x, ii, jj)

static void _fmpz_mat_charpoly_small_2x2(fmpz *rop, const fmpz_mat_t x)
{
    fmpz_add   (rop + 1, MAT(0, 0), MAT(1, 1));
    fmpz_neg   (rop + 1, rop + 1);
    fmpz_mul   (rop + 0, MAT(0, 0), MAT(1, 1));
    fmpz_submul(rop + 0, MAT(0, 1), MAT(1, 0));
}

static void _fmpz_mat_charpoly_small_3x3(fmpz *rop, const fmpz_mat_t x)
{
    fmpz a[2];
    fmpz_init(a + 0);
    fmpz_init(a + 1);

    fmpz_mul(   a + 0,   MAT(1, 0), MAT(2, 1));
    fmpz_submul(a + 0,   MAT(1, 1), MAT(2, 0));
    fmpz_mul(   rop + 0, a + 0,    MAT(0, 2));
    fmpz_neg(   rop + 0, rop + 0);
    fmpz_mul(   rop + 1, MAT(2, 0), MAT(0, 2));
    fmpz_neg(   rop + 1, rop + 1);

    fmpz_mul(   a + 0,   MAT(1, 2), MAT(2, 0));
    fmpz_submul(a + 0,   MAT(1, 0), MAT(2, 2));
    fmpz_submul(rop + 0, a + 0,    MAT(0, 1));
    fmpz_submul(rop + 1, MAT(1, 0), MAT(0, 1));

    fmpz_mul(   a + 0,   MAT(1, 1), MAT(2, 2));
    fmpz_add(   a + 1,   MAT(1, 1), MAT(2, 2));
    fmpz_neg(   a + 1,   a + 1);
    fmpz_submul(a + 0,   MAT(1, 2), MAT(2, 1));

    fmpz_submul(rop + 0, a + 0,    MAT(0, 0));
    fmpz_submul(rop + 1, a + 1,    MAT(0, 0));
    fmpz_add(   rop + 1, rop + 1,  a + 0);
    fmpz_sub(   rop + 2, a + 1,    MAT(0, 0));

    fmpz_clear(a + 0);
    fmpz_clear(a + 1);
}

#undef MAT

static int
want_berkowitz(slong n, const fmpz_mat_t op)
{
    if (n <= 8)
        return 1;

    slong bits = fmpz_mat_max_bits(op);
    bits = FLINT_ABS(bits);

    return n < (slong) FLINT_BIT_COUNT(bits);
}

void _fmpz_mat_charpoly(fmpz * rop, const fmpz_mat_t op)
{
    const slong n = op->r;

    if (n <= 3)
    {
        fmpz_one(rop + n);

        if (n == 1)
            fmpz_neg(rop, op->entries);
        else if (n == 2)
            _fmpz_mat_charpoly_small_2x2(rop, op);
        else if (n == 3)
            _fmpz_mat_charpoly_small_3x3(rop, op);
        return;
    }

    if (fmpz_mat_is_zero(op))
    {
        _fmpz_vec_zero(rop, n);
        fmpz_one(rop + n);
        return;
    }

    if (want_berkowitz(n, op))
        _fmpz_mat_charpoly_berkowitz(rop, op);
    else
        _fmpz_mat_charpoly_modular(rop, op);
}

void fmpz_mat_charpoly(fmpz_poly_t cp, const fmpz_mat_t mat)
{
    if (mat->r != mat->c)
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_charpoly).  Non-square matrix.\n");

     fmpz_poly_fit_length(cp, mat->r + 1);
    _fmpz_poly_set_length(cp, mat->r + 1);
    _fmpz_mat_charpoly(cp->coeffs, mat);
}

