/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_log(acb_t r, const acb_t z, slong prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)

    if (arb_is_zero(b))
    {
        if (arb_is_positive(a))
        {
            arb_log(acb_realref(r), a, prec);
            arb_zero(acb_imagref(r));
        }
        else if (arb_is_negative(a))
        {
            arb_neg(acb_realref(r), a);
            arb_log(acb_realref(r), acb_realref(r), prec);
            arb_const_pi(acb_imagref(r), prec);
        }
        else
        {
            acb_indeterminate(r);
        }
    }
    else if (arb_is_zero(a))
    {
        if (arb_is_positive(b))
        {
            arb_log(acb_realref(r), b, prec);
            arb_const_pi(acb_imagref(r), prec);
            arb_mul_2exp_si(acb_imagref(r), acb_imagref(r), -1);
        }
        else if (arb_is_negative(b))
        {
            arb_neg(acb_realref(r), b);
            arb_log(acb_realref(r), acb_realref(r), prec);
            arb_const_pi(acb_imagref(r), prec);
            arb_mul_2exp_si(acb_imagref(r), acb_imagref(r), -1);
            arb_neg(acb_imagref(r), acb_imagref(r));
        }
        else
        {
            acb_indeterminate(r);
        }
    }
    else
    {
        if (r != z)
        {
            arb_log_hypot(acb_realref(r), a, b, prec);
            if (arb_is_finite(acb_realref(r)))
                arb_atan2(acb_imagref(r), b, a, prec);
            else
                arb_indeterminate(acb_imagref(r));
        }
        else
        {
            arb_t t;
            arb_init(t);
            arb_log_hypot(t, a, b, prec);
            if (arb_is_finite(t))
                arb_atan2(acb_imagref(r), b, a, prec);
            else
                arb_indeterminate(acb_imagref(r));
            arb_swap(acb_realref(r), t);
            arb_clear(t);
        }
    }
#undef a
#undef b
}

void
acb_log_analytic(acb_t res, const acb_t z, int analytic, slong prec)
{
    if (analytic && arb_contains_zero(acb_imagref(z)) &&
        !arb_is_positive(acb_realref(z)))
    {
        acb_indeterminate(res);
    }
    else
    {
        acb_log(res, z, prec);
    }
}

