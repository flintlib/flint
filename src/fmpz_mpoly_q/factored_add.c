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

/* m = m * p^e, where m_is_one indicates that m is currently 1 */
static void
_mul_pow(fmpz_mpoly_t m, int * m_is_one, const fmpz_mpoly_t p, const fmpz_t e,
    fmpz_mpoly_t t, const fmpz_mpoly_ctx_t ctx)
{
    if (!fmpz_mpoly_pow_fmpz(t, p, e, ctx))
        flint_throw(FLINT_ERROR, "fmpz_mpoly_q_factored_add: exponent too large\n");

    if (*m_is_one)
    {
        fmpz_mpoly_swap(m, t, ctx);
        *m_is_one = 0;
    }
    else
    {
        fmpz_mpoly_mul(m, m, t, ctx);
    }
}

/*
    xn / xd + yn / yd with xn coprime to xd and yn coprime to yd.

    We bring both terms to the common denominator L = lcm(xd, yd), formed
    by taking the maximum exponent of each irreducible factor. If p has
    different exponents in xd and yd, then p does not divide the new
    numerator xn (L/xd) + yn (L/yd) (exactly one of the two terms is
    divisible by p). Thus only factors with equal exponents in xd and yd
    need to be checked for cancellation. The same argument shows that
    any common integer content must divide gcd(cx, cy) where cx, cy are
    the constant parts of xd, yd.
*/
static void
_fmpz_mpoly_q_factored_add_sub(fmpz_mpoly_q_factored_t res,
    const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y,
    int subtract, const fmpz_mpoly_ctx_t ctx)
{
    const fmpz_mpoly_factor_struct * xd = fmpz_mpoly_q_factored_denref(x);
    const fmpz_mpoly_factor_struct * yd = fmpz_mpoly_q_factored_denref(y);
    const fmpz_mpoly_struct * xn = fmpz_mpoly_q_factored_numref(x);
    const fmpz_mpoly_struct * yn = fmpz_mpoly_q_factored_numref(y);
    fmpz_mpoly_t mx, my, t, u, N;
    int mx_is_one, my_is_one;
    fmpz_mpoly_factor_t D;
    slong * cand;
    slong i, j, k, ncand;
    fmpz_t e, g, h;
    int c;

    if (fmpz_mpoly_q_factored_is_zero(x, ctx))
    {
        if (subtract)
            fmpz_mpoly_q_factored_neg(res, y, ctx);
        else
            fmpz_mpoly_q_factored_set(res, y, ctx);
        return;
    }

    if (fmpz_mpoly_q_factored_is_zero(y, ctx))
    {
        fmpz_mpoly_q_factored_set(res, x, ctx);
        return;
    }

    /* Fast path: adding two polynomials */
    if (xd->num == 0 && yd->num == 0 && fmpz_is_one(xd->constant) && fmpz_is_one(yd->constant))
    {
        if (subtract)
            fmpz_mpoly_sub(fmpz_mpoly_q_factored_numref(res), xn, yn, ctx);
        else
            fmpz_mpoly_add(fmpz_mpoly_q_factored_numref(res), xn, yn, ctx);
        fmpz_mpoly_factor_one(fmpz_mpoly_q_factored_denref(res), ctx);
        return;
    }

    fmpz_mpoly_init(mx, ctx);
    fmpz_mpoly_init(my, ctx);
    fmpz_mpoly_init(t, ctx);
    fmpz_mpoly_init(u, ctx);
    fmpz_mpoly_init(N, ctx);
    fmpz_mpoly_factor_init(D, ctx);
    fmpz_mpoly_factor_fit_length(D, xd->num + yd->num, ctx);
    fmpz_init(e);
    fmpz_init(g);
    fmpz_init(h);
    cand = flint_malloc(sizeof(slong) * (FLINT_MIN(xd->num, yd->num) + 1));
    mx_is_one = my_is_one = 1;
    ncand = 0;

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
            _fmpz_mpoly_q_factored_append(D, xd->poly + i, xd->exp + i, ctx);
            _mul_pow(my, &my_is_one, xd->poly + i, xd->exp + i, t, ctx);
            i++;
        }
        else if (c > 0)
        {
            _fmpz_mpoly_q_factored_append(D, yd->poly + j, yd->exp + j, ctx);
            _mul_pow(mx, &mx_is_one, yd->poly + j, yd->exp + j, t, ctx);
            j++;
        }
        else
        {
            c = fmpz_cmp(xd->exp + i, yd->exp + j);

            if (c < 0)
            {
                fmpz_sub(e, yd->exp + j, xd->exp + i);
                _mul_pow(mx, &mx_is_one, xd->poly + i, e, t, ctx);
                _fmpz_mpoly_q_factored_append(D, yd->poly + j, yd->exp + j, ctx);
            }
            else if (c > 0)
            {
                fmpz_sub(e, xd->exp + i, yd->exp + j);
                _mul_pow(my, &my_is_one, xd->poly + i, e, t, ctx);
                _fmpz_mpoly_q_factored_append(D, xd->poly + i, xd->exp + i, ctx);
            }
            else
            {
                cand[ncand++] = D->num;
                _fmpz_mpoly_q_factored_append(D, xd->poly + i, xd->exp + i, ctx);
            }

            i++;
            j++;
        }
    }

    /* Constants: L = lcm(cx, cy) = cx * (cy / g) */
    fmpz_gcd(g, xd->constant, yd->constant);

    /* t = xn * mx * (cy / g) */
    if (mx_is_one)
        fmpz_mpoly_set(t, xn, ctx);
    else
        fmpz_mpoly_mul(t, xn, mx, ctx);
    if (!fmpz_equal(g, yd->constant))
    {
        fmpz_divexact(e, yd->constant, g);
        fmpz_mpoly_scalar_mul_fmpz(t, t, e, ctx);
    }

    /* u = yn * my * (cx / g) */
    if (my_is_one)
        fmpz_mpoly_set(u, yn, ctx);
    else
        fmpz_mpoly_mul(u, yn, my, ctx);
    if (!fmpz_equal(g, xd->constant))
    {
        fmpz_divexact(e, xd->constant, g);
        fmpz_mpoly_scalar_mul_fmpz(u, u, e, ctx);
        fmpz_mul(D->constant, e, yd->constant);
    }
    else
    {
        fmpz_set(D->constant, yd->constant);
    }

    if (subtract)
        fmpz_mpoly_sub(N, t, u, ctx);
    else
        fmpz_mpoly_add(N, t, u, ctx);

    if (fmpz_mpoly_is_zero(N, ctx))
    {
        fmpz_mpoly_q_factored_zero(res, ctx);
    }
    else
    {
        /* Cancel common polynomial factors. */
        for (k = 0; k < ncand; k++)
            _fmpz_mpoly_q_factored_divide_out(N, D->exp + cand[k], D->poly + cand[k], ctx);

        if (ncand != 0)
            _fmpz_mpoly_q_factored_sort_merge(D, ctx);

        /* Cancel common content. */
        _fmpz_vec_content2(h, N->coeffs, N->length, g);
        if (!fmpz_is_one(h))
        {
            fmpz_mpoly_scalar_divexact_fmpz(N, N, h, ctx);
            fmpz_divexact(D->constant, D->constant, h);
        }

        fmpz_mpoly_swap(fmpz_mpoly_q_factored_numref(res), N, ctx);
        fmpz_mpoly_factor_swap(fmpz_mpoly_q_factored_denref(res), D, ctx);
    }

    fmpz_mpoly_clear(mx, ctx);
    fmpz_mpoly_clear(my, ctx);
    fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_clear(u, ctx);
    fmpz_mpoly_clear(N, ctx);
    fmpz_mpoly_factor_clear(D, ctx);
    fmpz_clear(e);
    fmpz_clear(g);
    fmpz_clear(h);
    flint_free(cand);
}

void
fmpz_mpoly_q_factored_add(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_q_factored_add_sub(res, x, y, 0, ctx);
}

void
fmpz_mpoly_q_factored_sub(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_q_factored_add_sub(res, x, y, 1, ctx);
}
