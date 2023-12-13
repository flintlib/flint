/*
    Copyright (C) 2011, 2012 Fredrik Johansson
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr_poly.h"

void
_fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(fmpz * vs, const fmpz * poly,
    slong plen, fmpz_poly_struct * const * tree, slong len, const fmpz_mod_ctx_t ctx)
{
    slong height, i, j, pow, left;
    slong tree_height;
    fmpz_t temp, inv;
    fmpz * t, * u, * pb, * pc, * swap;
    fmpz_poly_struct * pa;

    fmpz_init(temp);
    fmpz_init(inv);

    /* avoid worrying about some degenerate cases */
    if (len < 2 || plen < 2)
    {
        if (len == 1)
        {
            fmpz_mod_neg(temp, tree[0]->coeffs, ctx);
            _fmpz_mod_poly_evaluate_fmpz(vs, poly, plen, temp, ctx);
        } else if (len != 0 && plen == 0)
            _fmpz_vec_zero(vs, len);
        else if (len != 0 && plen == 1)
            for (i = 0; i < len; i++)
                fmpz_set(vs + i, poly);

        fmpz_clear(temp);
        return;
    }

    t = _fmpz_vec_init(2*len);
    u = _fmpz_vec_init(2*len);

    left = len;

    /* Initial reduction. We allow the polynomial to be larger
       or smaller than the number of points. */
    height = FLINT_BIT_COUNT(plen - 1) - 1;
    tree_height = FLINT_CLOG2(len);
    while (height >= tree_height)
        height--;
    pow = WORD(1) << height;

    for (i = j = 0; i < len; i += pow, j++)
    {
        pa = tree[height] + j;
        fmpz_mod_inv(inv, pa->coeffs + pa->length - 1, ctx);
        _fmpz_mod_poly_rem(t + i, poly, plen, pa->coeffs, pa->length, inv, ctx);
    }

    for (i = height - 1; i >= 0; i--)
    {
        pow = WORD(1) << i;
        left = len;
        pa = tree[i];
        pb = t;
        pc = u;

        left = len;
        while (left >= 2 * pow)
        {
            fmpz_mod_inv(inv, pa->coeffs + pa->length - 1, ctx);
            _fmpz_mod_poly_rem(pc, pb, 2 * pow, pa->coeffs, pa->length, inv, ctx);

            pa++;
            fmpz_mod_inv(inv, pa->coeffs + pa->length - 1, ctx);
            _fmpz_mod_poly_rem(pc + pow, pb, 2 * pow, pa->coeffs, pa->length, inv, ctx);

            pa++;
            pb += 2 * pow;
            pc += 2 * pow;
            left -= 2 * pow;
        }

        if (left > pow)
        {
            fmpz_mod_inv(inv, pa->coeffs + pa->length - 1, ctx);
            _fmpz_mod_poly_rem(pc, pb, left, pa->coeffs, pa->length, inv, ctx);

            pa ++;
            fmpz_mod_inv(inv, pa->coeffs + pa->length - 1, ctx);
            _fmpz_mod_poly_rem(pc + pow, pb, left, pa->coeffs, pa->length, inv, ctx);
        }
        else if (left > 0)
            _fmpz_vec_set(pc, pb, left);

        swap = t;
        t = u;
        u = swap;
    }

    fmpz_clear(temp);
    fmpz_clear(inv);

    _fmpz_vec_set(vs, t, len);

    _fmpz_vec_clear(t, 2*len);
    _fmpz_vec_clear(u, 2*len);
}

void _fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys, const fmpz * poly, slong plen,
    const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_evaluate_vec_fast(ys, poly, plen, xs, n, gr_ctx));
}

void fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys, const fmpz_mod_poly_t poly,
                            const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_poly_evaluate_fmpz_vec_fast(ys, poly->coeffs,
                               poly->length, xs, n, ctx);
}

void
_fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys, const fmpz * coeffs, slong len,
    const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
{
    slong i;
    for (i = 0; i < n; i++)
        _fmpz_mod_poly_evaluate_fmpz(ys + i, coeffs, len, xs + i, ctx);
}

void fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys, const fmpz_mod_poly_t poly,
                           const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_poly_evaluate_fmpz_vec_iter(ys, poly->coeffs,
                               poly->length, xs, n, ctx);
}

void
_fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys, const fmpz * coeffs,
                        slong len, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
{
    if (len < FMPZ_MOD_POLY_EVALUATE_FMPZ_VEC)
        _fmpz_mod_poly_evaluate_fmpz_vec_iter(ys, coeffs, len, xs, n, ctx);
    else
        _fmpz_mod_poly_evaluate_fmpz_vec_fast(ys, coeffs, len, xs, n, ctx);
}

void fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys, const fmpz_mod_poly_t poly,
                            const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_poly_evaluate_fmpz_vec(ys, poly->coeffs, poly->length, xs, n, ctx);
}
