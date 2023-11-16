/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "acb_poly.h"

#ifndef __compar_fn_t
#if defined(_MSC_VER)
typedef int(*__compar_fn_t) (const void *, const void *);
#else
typedef int(*__compar_fn_t) (__const void *, __const void *);
#endif
#endif

int arb_cmp_mid(const arb_t a, const arb_t b)
{
    return arf_cmp(arb_midref(a), arb_midref(b));
}

void _arb_vec_sort_mid(arb_ptr vec, slong len)
{
    qsort(vec, len, sizeof(arb_struct), (__compar_fn_t) arb_cmp_mid);
}

int
_acb_poly_validate_real_roots(acb_srcptr roots, acb_srcptr poly, slong len, slong prec)
{
    slong i, deg, num_real;
    arb_ptr real;
    int result;

    deg = len - 1;
    num_real = 0;
    result = 1;

    if (deg <= 1)
        return 1;

    real = _arb_vec_init(deg);

    /* pick out the candidate real roots */
    for (i = 0; i < deg; i++)
    {
        if (arb_contains_zero(acb_imagref(roots + i)))
        {
            arb_set(real + num_real, acb_realref(roots + i));
            num_real++;
        }
    }

    /* number of real roots must be even if the polynomial is even,
       and odd if the polynomial is odd (unless there are repeated roots...
       in which case the input is invalid) */
    if ((num_real % 2) != (deg % 2))
    {
        result = 0;
    }
    else if (num_real > 0)
    {
        int sign_neg_inf, sign_pos_inf, prev_sign;

        acb_t t;
        acb_init(t);

        /* by assumption that the roots are real and isolated, the lead
           coefficient really must be known to be either positive or negative */
        sign_pos_inf = arb_is_positive(acb_realref(poly + deg)) ? 1 : -1;
        sign_neg_inf = (deg % 2) ? -sign_pos_inf : sign_pos_inf;

        /* now we check that there's a sign change between each root */
        _arb_vec_sort_mid(real, num_real);

        prev_sign = sign_neg_inf;

        for (i = 0; i < num_real - 1; i++)
        {
            /* set t to the midpoint between the midpoints */
            arb_zero(acb_imagref(t));
            arf_add(arb_midref(acb_realref(t)),
                arb_midref(real + i), arb_midref(real + i + 1), prec, ARF_RND_DOWN);
            arf_mul_2exp_si(arb_midref(acb_realref(t)), arb_midref(acb_realref(t)), -1);
            mag_zero(arb_radref(acb_realref(t)));

            /* check that this point really is between both intervals (one interval
               could be much wider than the other */
            if (arb_lt(real + i, acb_realref(t)) && arb_lt(acb_realref(t), real + i + 1))
            {
                /* check sign change */
                _acb_poly_evaluate(t, poly, len, t, prec);

                if (prev_sign == 1)
                    result = arb_is_negative(acb_realref(t));
                else
                    result = arb_is_positive(acb_realref(t));

                if (!result)
                    break;

                prev_sign = -prev_sign;
            }
            else
            {
                result = 0;
                break;
            }
        }

        acb_clear(t);
    }

    _arb_vec_clear(real, deg);

    return result;
}

int
acb_poly_validate_real_roots(acb_srcptr roots, const acb_poly_t poly, slong prec)
{
    return _acb_poly_validate_real_roots(roots, poly->coeffs, poly->length, prec);
}

