/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_modular.h"

void
acb_modular_fundamental_domain_approx_arf(psl2z_t g,
    const arf_t xx, const arf_t yy, const arf_t one_minus_eps, slong prec)
{
    slong i;
    arf_t x, y, t;
    fmpz_t n;

    psl2z_one(g);

    /* we must be in the upper half-plane */
    if (!arf_is_finite(xx) || !arf_is_finite(yy) || arf_sgn(yy) <= 0)
        return;

    arf_init(x);
    arf_init(y);
    arf_init(t);
    fmpz_init(n);

    arf_set_round(x, xx, prec, ARF_RND_DOWN);
    arf_set_round(y, yy, prec, ARF_RND_DOWN);

    for (i = 0; i < 10 + prec / 4; i++)
    {
        /* too large shift */
        if (arf_cmpabs_2exp_si(x, prec) > 0)
        {
            psl2z_one(g);
            break;
        }

        if (fmpz_bits(&g->a) > prec || fmpz_bits(&g->b) > prec ||
            fmpz_bits(&g->c) > prec || fmpz_bits(&g->d) > prec)
        {
            psl2z_one(g);
            break;
        }

        /* shift */
        if (arf_cmpabs_2exp_si(x, -1) > 0)
        {
            arf_get_fmpz(n, x, ARF_RND_NEAR);
            arf_sub_fmpz(x, x, n, prec, ARF_RND_DOWN);
            fmpz_submul(&g->a, &g->c, n);
            fmpz_submul(&g->b, &g->d, n);
            continue;
        }

        arf_mul(t, x, x, prec, ARF_RND_DOWN);
        arf_addmul(t, y, y, prec, ARF_RND_DOWN);

        /* inversion */
        if (arf_cmp(t, one_minus_eps) < 0)
        {
            arf_div(x, x, t, prec, ARF_RND_DOWN);
            arf_neg(x, x);
            arf_div(y, y, t, prec, ARF_RND_DOWN);

            fmpz_swap(&g->a, &g->c);
            fmpz_swap(&g->b, &g->d);
            fmpz_neg(&g->a, &g->a);
            fmpz_neg(&g->b, &g->b);

            continue;
        }

        /* we're good */
        break;
    }

    if (fmpz_sgn(&g->c) < 0 || (fmpz_is_zero(&g->c) && fmpz_sgn(&g->d) < 0))
    {
        fmpz_neg(&g->a, &g->a);
        fmpz_neg(&g->b, &g->b);
        fmpz_neg(&g->c, &g->c);
        fmpz_neg(&g->d, &g->d);
    }

    arf_clear(x);
    arf_clear(y);
    arf_clear(t);
    fmpz_clear(n);
}

