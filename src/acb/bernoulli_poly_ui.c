/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_bernoulli_poly_ui(acb_t res, ulong n, const acb_t x, slong prec)
{
    acb_t s, x2;
    arb_t t, c;
    ulong k;

    if (n == 0)
    {
        acb_one(res);
        return;
    }

    if (n == 1)
    {
        acb_mul_2exp_si(res, x, 1);
        acb_sub_ui(res, res, 1, prec);
        acb_mul_2exp_si(res, res, -1);
        return;
    }

    if (acb_is_real(x))
    {
        arb_bernoulli_poly_ui(acb_realref(res), n, acb_realref(x), prec);
        arb_zero(acb_imagref(res));
        return;
    }

    /* assuming small n simplifies the code that follows */
    if (n >> (FLINT_BITS / 2) || !acb_is_finite(x))
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(s);
    acb_init(x2);
    arb_init(t);
    arb_init(c);

    acb_mul(x2, x, x, prec);

    /* s = x^2 - x n / 2 */
    acb_mul_ui(s, x, n, prec);
    acb_mul_2exp_si(s, s, -1);
    acb_sub(s, x2, s, prec);

    /* c = n (n-1) / 2;  s = s + c / 6 */
    arb_set_ui(c, n * (n - 1));
    arb_mul_2exp_si(c, c, -1);
    arb_div_ui(t, c, 6, prec);
    acb_add_arb(s, s, t, prec);

    for (k = 4; k <= n; k += 2)
    {
        /* c = binomial(n,k) */
        arb_mul_ui(c, c, (n + 1 - k) * (n + 2 - k), prec);
        arb_div_ui(c, c, k * (k - 1), prec);

        /* s = s x^2 + b_k c */
        acb_mul(s, s, x2, prec);
        arb_bernoulli_ui(t, k, prec);
        arb_mul(t, t, c, prec);
        acb_add_arb(s, s, t, prec);
    }

    if (n >= 3 && n % 2)
        acb_mul(s, s, x, prec);

    acb_swap(res, s);

    acb_clear(s);
    acb_clear(x2);
    arb_clear(t);
    arb_clear(c);
}

