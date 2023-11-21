/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

int
_acb_lambertw_check_branch(const acb_t w, const fmpz_t k, slong prec)
{
    arb_t t, u, v, a, b;
    int res = 0;

    arb_init(t);
    arb_init(u);
    arb_init(v);
    arb_init(a);
    arb_init(b);

    /* t = x sinc(y), v = -cos(y) */
    if (arb_is_exact(acb_imagref(w)))
    {
        if (arb_is_zero(acb_imagref(w)))
        {
            arb_one(t);
            arb_one(v);
        }
        else
        {
            arb_sin_cos(t, v, acb_imagref(w), prec);
            arb_div(t, t, acb_imagref(w), prec);
        }
    }
    else
    {
        arb_sinc(t, acb_imagref(w), prec);
        arb_cos(v, acb_imagref(w), prec);
    }
    arb_mul(t, t, acb_realref(w), prec);
    arb_neg(v, v);

    /* u = y / pi, with conjugate relation for k */
    arb_const_pi(u, prec);
    arb_div(u, acb_imagref(w), u, prec);
    if (fmpz_sgn(k) < 0)
        arb_neg(u, u);

    if (fmpz_is_zero(k))
    {
        /* -1 < u < 1 and t > v */
        arb_set_si(a, -1);
        arb_set_si(b, 1);

        if (arb_gt(u, a) && arb_lt(u, b) && arb_gt(t, v))
        {
            res = 1;
        }
    }
    else
    {
        arb_set_fmpz(a, k);
        arb_abs(a, a);
        arb_mul_2exp_si(a, a, 1);
        arb_add_ui(b, a, 1, prec);
        arb_sub_ui(a, a, 2, prec);

        /* if u > 2|k|-2 and u < 2|k|+1 */
        if (arb_gt(u, a) && arb_lt(u, b))
        {
            arb_add_ui(a, a, 1, prec);
            arb_sub_ui(b, b, 1, prec);

            /* if u > 2|k|-1 and u < 2|k| */
            if (arb_gt(u, a) && arb_lt(u, b))
            {
                res = 1;
            }
            else if (arb_lt(u, b) && arb_lt(t, v)) /* u < 2|k| and t < v */
            {
                res = 1;
            }
            else if (arb_gt(u, a) && arb_gt(t, v)) /* u > 2|k|-1 and t > v */
            {
                res = 1;
            }
        }
    }

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    arb_clear(a);
    arb_clear(b);

    return res;
}

int
acb_lambertw_check_branch(const acb_t w, const fmpz_t k, slong prec)
{
    if (prec > 64 && _acb_lambertw_check_branch(w, k, 32))
        return 1;

    return _acb_lambertw_check_branch(w, k, prec);
}

