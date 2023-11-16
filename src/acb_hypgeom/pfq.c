/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

static int
_acb_vec_maybe_nonpositive_int(acb_srcptr b, slong q)
{
    slong i;

    for (i = 0; i < q; i++)
        if (!arb_is_positive(acb_realref(b + i)) && acb_contains_int(b + i))
            return 1;

    return 0;
}

void
acb_hypgeom_pfq(acb_t res, acb_srcptr a, slong p,
                           acb_srcptr b, slong q, const acb_t z, int regularized, slong prec)
{
    if (p == 0 && q == 0)
    {
        acb_exp(res, z, prec);
    }
    else if (p == 1 && q == 0)
    {
        acb_t t;
        acb_init(t);
        acb_neg(t, a);
        acb_sub_ui(res, z, 1, prec);
        acb_neg(res, res);
        acb_pow(res, res, t, prec);
        acb_clear(t);
    }
    else if (p == 0 && q == 1)
    {
        acb_hypgeom_0f1(res, b, z, regularized, prec);
    }
    else if (p == 1 && q == 1)
    {
        acb_hypgeom_m(res, a, b, z, regularized, prec);
    }
    else if (p == 2 && q == 1)
    {
        acb_hypgeom_2f1(res, a, a + 1, b, z, regularized, prec);
    }
    else if (regularized && _acb_vec_maybe_nonpositive_int(b, q))
    {
        /* todo: implement regularized sum without using polynomials */
        acb_poly_struct * tmp;
        slong i;

        tmp = flint_malloc(sizeof(acb_poly_struct) * (p + q + 2));
        for (i = 0; i < p + q + 2; i++)
            acb_poly_init(tmp + i);

        for (i = 0; i < p; i++)
            acb_poly_set_acb(tmp + i, a + i);
        for (i = 0; i < q; i++)
            acb_poly_set_acb(tmp + p + i, b + i);
        acb_poly_one(tmp + p + q);
        acb_poly_set_acb(tmp + p + q + 1, z);

        acb_hypgeom_pfq_series_direct(tmp, tmp, p, tmp + p, q + 1,
            tmp + p + q + 1, regularized, -1, 1, prec);

        acb_poly_get_coeff_acb(res, tmp, 0);

        for (i = 0; i < p + q + 2; i++)
            acb_poly_clear(tmp + i);
        flint_free(tmp);
    }
    else
    {
        acb_ptr tmp;
        slong i, j, alloc = 0;

        /* check if we can remove a '1' from the upper parameters */
        for (i = 0; i < p; i++)
        {
            if (acb_is_one(a + i))
            {
                alloc = p;
                tmp = _acb_vec_init(alloc);
                for (j = 0; j < p - 1; j++)
                    acb_set(tmp + 1 + j, a + j + (j >= i));
                acb_hypgeom_pfq_direct(tmp, tmp + 1, p - 1, b, q, z, -1, prec);
                break;
            }
        }

        if (alloc == 0)
        {
            alloc = q + 2;
            tmp = _acb_vec_init(alloc);

            for (j = 0; j < q; j++)
                acb_set(tmp + 1 + j, b + j);
            acb_one(tmp + 1 + q);
            acb_hypgeom_pfq_direct(tmp, a, p, tmp + 1, q + 1, z, -1, prec);
        }

        if (regularized && q > 0)
        {
            acb_t c, t;
            acb_init(c);
            acb_init(t);
            acb_rgamma(c, b, prec);

            for (i = 1; i < q; i++)
            {
                acb_rgamma(t, b + i, prec);
                acb_mul(c, c, t, prec);
            }

            acb_mul(tmp, tmp, c, prec);

            acb_clear(c);
            acb_clear(t);
        }

        acb_set(res, tmp);
        _acb_vec_clear(tmp, alloc);
    }

    if (!acb_is_finite(res))
        acb_indeterminate(res);
}

