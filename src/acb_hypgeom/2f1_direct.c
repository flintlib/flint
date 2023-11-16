/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_2f1_direct(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int regularized, slong prec)
{
    /* 2F1R(a,b,-n,z) = (a)_(n+1) * (b)_(n+1) * z^(n+1) / (n+1)!
                           * 2F1(a+n+1, b+n+1, n+2, z) */
    if (regularized && acb_is_int(c)
            && arf_sgn(arb_midref(acb_realref(c))) <= 0)
    {
        if ((acb_is_int(a) && arb_is_nonpositive(acb_realref(a)) &&
            arf_cmp(arb_midref(acb_realref(a)), arb_midref(acb_realref(c))) >= 0) ||
            (acb_is_int(b) && arb_is_nonpositive(acb_realref(b)) &&
            arf_cmp(arb_midref(acb_realref(b)), arb_midref(acb_realref(c))) >= 0))
        {
            acb_zero(res);
        }
        else
        {
            acb_t n, n1, t, u, v;
            acb_ptr aa;
            int p, q;

            acb_init(n);
            acb_init(n1);
            acb_init(t);
            acb_init(u);
            acb_init(v);
            aa = _acb_vec_init(4);

            acb_neg(n, c);
            acb_add_ui(n1, n, 1, prec);

            acb_add(aa, a, n1, prec);
            acb_add(aa + 1, b, n1, prec);
            acb_add_ui(aa + 2, n1, 1, prec);

            if (acb_is_one(aa))
            {
                p = q = 1;
                acb_swap(aa, aa + 1);
            }
            else if (acb_is_one(aa + 1))
            {
                p = q = 1;
            }
            else
            {
                p = q = 2;
                acb_one(aa + 3);
            }

            acb_hypgeom_pfq_direct(t, aa, p, aa + 2, q, z, -1, prec);

            /* z^(n+1) */
            acb_pow(u, z, n1, prec);
            acb_mul(t, t, u, prec);

            acb_rising(u, a, n1, prec);
            acb_mul(t, t, u, prec);

            acb_rising(u, b, n1, prec);
            acb_mul(t, t, u, prec);

            /* 1/(n+1)! */
            acb_rgamma(u, aa + 2, prec);
            acb_mul(t, t, u, prec);

            acb_set(res, t);

            acb_clear(n);
            acb_clear(n1);
            acb_clear(t);
            acb_clear(u);
            acb_clear(v);
            _acb_vec_clear(aa, 4);
        }
    }
    else
    {
        acb_ptr aa;
        int p, q;

        aa = _acb_vec_init(4);

        acb_set(aa + 2, c);

        if (acb_is_one(a))
        {
            p = q = 1;
            acb_set(aa, b);
        }
        else if (acb_is_one(b))
        {
            p = q = 1;
            acb_set(aa, a);
        }
        else
        {
            p = q = 2;
            acb_set(aa, a);
            acb_set(aa + 1, b);
            acb_one(aa + 3);
        }

        acb_hypgeom_pfq_direct(res, aa, p, aa + 2, q, z, -1, prec);

        if (regularized)
        {
            acb_rgamma(aa + 2, aa + 2, prec);
            acb_mul(res, res, aa + 2, prec);
        }

        _acb_vec_clear(aa, 4);
    }

    if (!acb_is_finite(res))
        acb_indeterminate(res);
}

