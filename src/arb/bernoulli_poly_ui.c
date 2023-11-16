/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_bernoulli_poly_ui(arb_t res, ulong n, const arb_t x, slong prec)
{
    arb_t s, t, c, x2;
    ulong m, k;
    int negate;

    if (n == 0)
    {
        arb_one(res);
        return;
    }

    if (n == 1)
    {
        arb_mul_2exp_si(res, x, 1);
        arb_sub_ui(res, res, 1, prec);
        arb_mul_2exp_si(res, res, -1);
        return;
    }

    /* small integer x */
    if (arb_is_int(x) && arf_cmpabs_ui(arb_midref(x), n) < 0 && n < WORD_MAX)
    {
        if (arf_sgn(arb_midref(x)) >= 0)
        {
            m = arf_get_si(arb_midref(x), ARF_RND_DOWN);
            negate = 0;
        }
        else
        {
            m = UWORD(1) - arf_get_si(arb_midref(x), ARF_RND_DOWN);
            negate = n % 2;
        }

        arb_init(t);
        arb_zero(res);

        /* todo: use a dedicated power sum function */
        for (k = 1; k < m; k++)
        {
            arb_ui_pow_ui(t, k, n - 1, prec);
            arb_add(res, res, t, prec);
        }

        arb_mul_ui(res, res, n, prec);
        arb_bernoulli_ui(t, n, prec);
        arb_add(res, res, t, prec);
        if (negate)
            arb_neg(res, res);

        arb_clear(t);
        return;
    }

    /* assuming small n simplifies the code that follows */
    if (n >> (FLINT_BITS / 2) || !arb_is_finite(x))
    {
        arb_indeterminate(res);
        return;
    }

    arb_init(s);
    arb_init(t);
    arb_init(c);
    arb_init(x2);

    arb_mul(x2, x, x, prec);

    /* s = x^2 - x n / 2 */
    arb_mul_ui(s, x, n, prec);
    arb_mul_2exp_si(s, s, -1);
    arb_sub(s, x2, s, prec);

    /* c = n (n-1) / 2;  s = s + c / 6 */
    arb_set_ui(c, n * (n - 1));
    arb_mul_2exp_si(c, c, -1);
    arb_div_ui(t, c, 6, prec);
    arb_add(s, s, t, prec);

    for (k = 4; k <= n; k += 2)
    {
        /* c = binomial(n,k) */
        arb_mul_ui(c, c, (n + 1 - k) * (n + 2 - k), prec);
        arb_div_ui(c, c, k * (k - 1), prec);

        /* s = s x^2 + b_k c */
        arb_mul(s, s, x2, prec);
        arb_bernoulli_ui(t, k, prec);
        arb_addmul(s, t, c, prec);
    }

    if (n >= 3 && n % 2)
        arb_mul(s, s, x, prec);

    arb_swap(res, s);

    arb_clear(s);
    arb_clear(t);
    arb_clear(c);
    arb_clear(x2);
}

