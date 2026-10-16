/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"

void
fmpz_mpoly_q_factored_neg(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_q_factored_set(res, x, ctx);
    fmpz_mpoly_neg(fmpz_mpoly_q_factored_numref(res), fmpz_mpoly_q_factored_numref(res), ctx);
}

/*
    (xn / xd) * (yn / yd) where xn is coprime to xd and yn is coprime to yd.

    The irreducible factors of the denominators are merged in sorted order.
    A factor p occurring in both xd and yd cannot divide xn or yn,
    so its exponents are simply added. A factor occurring only in xd can
    only cancel against yn, and a factor occurring only in yd can only
    cancel against xn; these cancellations are found by trial division.
*/
void
fmpz_mpoly_q_factored_mul(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)
{
    const fmpz_mpoly_factor_struct * xd = fmpz_mpoly_q_factored_denref(x);
    const fmpz_mpoly_factor_struct * yd = fmpz_mpoly_q_factored_denref(y);
    fmpz_mpoly_t xn, yn;
    fmpz_mpoly_factor_t D;
    fmpz_t e, g1, g2;
    slong i, j;
    int c;

    if (fmpz_mpoly_q_factored_is_zero(x, ctx) || fmpz_mpoly_q_factored_is_zero(y, ctx))
    {
        fmpz_mpoly_q_factored_zero(res, ctx);
        return;
    }

    /* Fast path: multiplying two polynomials */
    if (xd->num == 0 && yd->num == 0 && fmpz_is_one(xd->constant) && fmpz_is_one(yd->constant))
    {
        fmpz_mpoly_mul(fmpz_mpoly_q_factored_numref(res), fmpz_mpoly_q_factored_numref(x),
                                                        fmpz_mpoly_q_factored_numref(y), ctx);
        fmpz_mpoly_factor_one(fmpz_mpoly_q_factored_denref(res), ctx);
        return;
    }

    fmpz_mpoly_init(xn, ctx);
    fmpz_mpoly_init(yn, ctx);
    fmpz_mpoly_factor_init(D, ctx);
    fmpz_mpoly_factor_fit_length(D, xd->num + yd->num, ctx);
    fmpz_init(e);
    fmpz_init(g1);
    fmpz_init(g2);

    fmpz_mpoly_set(xn, fmpz_mpoly_q_factored_numref(x), ctx);
    fmpz_mpoly_set(yn, fmpz_mpoly_q_factored_numref(y), ctx);

    i = j = 0;
    while (i < xd->num || j < yd->num)
    {
        if (j >= yd->num)
            c = -1;
        else if (i >= xd->num)
            c = 1;
        else
            c = fmpz_mpoly_cmp(xd->poly + i, yd->poly + j, ctx);

        if (c < 0)
        {
            fmpz_set(e, xd->exp + i);
            _fmpz_mpoly_q_factored_divide_out(yn, e, xd->poly + i, ctx);
            if (!fmpz_is_zero(e))
                _fmpz_mpoly_q_factored_append(D, xd->poly + i, e, ctx);
            i++;
        }
        else if (c > 0)
        {
            fmpz_set(e, yd->exp + j);
            _fmpz_mpoly_q_factored_divide_out(xn, e, yd->poly + j, ctx);
            if (!fmpz_is_zero(e))
                _fmpz_mpoly_q_factored_append(D, yd->poly + j, e, ctx);
            j++;
        }
        else
        {
            fmpz_add(e, xd->exp + i, yd->exp + j);
            _fmpz_mpoly_q_factored_append(D, xd->poly + i, e, ctx);
            i++;
            j++;
        }
    }

    /* Cross-cancel integer content against the constant denominators. */
    _fmpz_vec_content2(g1, xn->coeffs, xn->length, yd->constant);
    _fmpz_vec_content2(g2, yn->coeffs, yn->length, xd->constant);

    if (!fmpz_is_one(g1))
    {
        fmpz_mpoly_scalar_divexact_fmpz(xn, xn, g1, ctx);
        fmpz_divexact(e, yd->constant, g1);
    }
    else
    {
        fmpz_set(e, yd->constant);
    }

    if (!fmpz_is_one(g2))
    {
        fmpz_mpoly_scalar_divexact_fmpz(yn, yn, g2, ctx);
        fmpz_divexact(D->constant, xd->constant, g2);
    }
    else
    {
        fmpz_set(D->constant, xd->constant);
    }

    fmpz_mul(D->constant, D->constant, e);

    fmpz_mpoly_mul(fmpz_mpoly_q_factored_numref(res), xn, yn, ctx);
    fmpz_mpoly_factor_swap(fmpz_mpoly_q_factored_denref(res), D, ctx);

    fmpz_mpoly_clear(xn, ctx);
    fmpz_mpoly_clear(yn, ctx);
    fmpz_mpoly_factor_clear(D, ctx);
    fmpz_clear(e);
    fmpz_clear(g1);
    fmpz_clear(g2);
}
