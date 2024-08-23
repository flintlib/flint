/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

/* These functions are more of a rule-of-thumb kind to increase performance. */

static slong
acb_theta_ql_split(const arb_mat_t cho)
{
    slong g = arb_mat_nrows(cho);
    arb_t cmp;
    slong k;

    arb_init(cmp);

    for (k = g - 1; k >= 1; k--)
    {
        arb_mul_2exp_si(cmp, arb_mat_entry(cho, k - 1, k - 1),
            FLINT_MAX(1, 6 + k - 2 * g));
        if (arb_lt(cmp, arb_mat_entry(cho, k, k)))
        {
            break;
        }
    }

    arb_clear(cmp);
    return k;
}

slong
acb_theta_ql_nb_steps(slong * split, const acb_theta_ctx_tau_t ctx, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t x, t;
    slong s, res;

    arb_init(x);
    arb_init(t);

    s = acb_theta_ql_split(acb_theta_ctx_cho(ctx));
    *split = s;

    arb_sqr(x, arb_mat_entry(acb_theta_ctx_cho(ctx), s, s), lp);
    arb_const_log2(t, lp);
    arb_div(x, x, t, lp);
    arb_div_si(x, x, prec, lp);
    arb_log(x, x, lp);
    arb_div(x, x, t, lp);

    if (!arb_is_finite(x) || arf_cmpabs_2exp_si(arb_midref(x), FLINT_BITS - 4) > 0)
    {
        arb_clear(x);
        arb_clear(t);
        return 0;
    }

    res =  -arf_get_si(arb_midref(x), ARF_RND_NEAR);
    if (s == 0)
    {
        if (g == 1)
        {
            res -= 7;
        }
        else if (g == 2)
        {
            res -= 3;
        }
        else if (g <= 5)
        {
            res -= 1;
        }
    }
    else
    {
        res += 1;
    }
    res = FLINT_MAX(0, res);

    arb_clear(x);
    arb_clear(t);
    return res;
}
