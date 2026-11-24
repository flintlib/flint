/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb.h"
#include "acb_hypgeom.h"

/* With complex parameters, the products of the factors a_i + k in
   rectangular splitting are balls that are multiplied by complex numbers
   at every step: the error boxes are rotated and grow each time (the
   wrapping effect), by a factor |cos(theta)| + |sin(theta)| for a
   rotation by theta = arg(a_i + k), about 1 + |Im(a_i)| / |Re(a_i) + k|.
   Over n terms this is a loss of about
   sum_i |Im(a_i)| log((|Re(a_i)| + n) / (|Re(a_i)| + 1)) nats, which is
   negligible for small parameters and thousands of bits for large ones.
   The choice is made if this exceeds a tenth of the precision.
   Binary splitting and fast multipoint evaluation, with few
   multiplications in sequence, lose much less. */
static int
_pfq_sum_wraps(acb_srcptr a, slong p, acb_srcptr b, slong q, slong n, slong prec)
{
    double loss = 0.0, re, im;
    slong i;

    for (i = 0; i < p + q; i++)
    {
        const acb_struct * x = (i < p) ? a + i : b + i - p;

        if (arb_is_zero(acb_imagref(x)))
            continue;

        re = fabs(arf_get_d(arb_midref(acb_realref(x)), ARF_RND_NEAR));
        im = fabs(arf_get_d(arb_midref(acb_imagref(x)), ARF_RND_NEAR));

        if (!(re < 1e300) || !(im < 1e300))
            return 0;

        loss += im * log((re + n) / (re + 1.0));
    }

    /* (in bits; the estimate is an upper bound, about twice what is
       observed) */
    return 0.5 * loss * 1.4426950408889634 > 0.1 * prec + 30.0;
}

void
acb_hypgeom_pfq_sum(acb_t s, acb_t t, acb_srcptr a, slong p,
    acb_srcptr b, slong q, const acb_t z, slong n, slong prec)
{
    if (n > 4 && prec >= 128
        && _acb_vec_bits(a, p) * p + _acb_vec_bits(b, q) * q + 10 < prec / 2)
    {
        if (prec >= 256 && acb_bits(z) < prec * 0.01)
            acb_hypgeom_pfq_sum_bs(s, t, a, p, b, q, z, n, prec);
        else if (_pfq_sum_wraps(a, p, b, q, n, prec))
        {
            if (prec >= 1500 && n >= 30 + 100000 / (prec - 1000))
                acb_hypgeom_pfq_sum_fme(s, t, a, p, b, q, z, n, prec);
            else
                acb_hypgeom_pfq_sum_bs(s, t, a, p, b, q, z, n, prec);
        }
        else
            acb_hypgeom_pfq_sum_rs(s, t, a, p, b, q, z, n, prec);
    }
    else if (prec >= 1500 && n >= 30 + 100000 / (prec - 1000))
    {
        acb_hypgeom_pfq_sum_fme(s, t, a, p, b, q, z, n, prec);
    }
    else
    {
        /* With inexact complex numbers, prefer binary splitting
           as it reduces the wrapping effect. */
        if (n > 8)
            acb_hypgeom_pfq_sum_bs(s, t, a, p, b, q, z, n, prec);
        else
            acb_hypgeom_pfq_sum_forward(s, t, a, p, b, q, z, n, prec);
    }
}

void
acb_hypgeom_pfq_sum_invz(acb_t s, acb_t t, acb_srcptr a, slong p,
    acb_srcptr b, slong q, const acb_t z, const acb_t zinv, slong n, slong prec)
{
    if (n > 4 && prec >= 128
        && _acb_vec_bits(a, p) * p + _acb_vec_bits(b, q) * q + 10 < prec / 2)
    {
        if (prec >= 256 && acb_bits(zinv) < prec * 0.01)
            acb_hypgeom_pfq_sum_bs_invz(s, t, a, p, b, q, zinv, n, prec);
        else if (_pfq_sum_wraps(a, p, b, q, n, prec))
        {
            if (prec >= 1500 && n >= 30 + 100000 / (prec - 1000))
                acb_hypgeom_pfq_sum_fme(s, t, a, p, b, q, z, n, prec);
            else
                acb_hypgeom_pfq_sum_bs_invz(s, t, a, p, b, q, zinv, n, prec);
        }
        else
            acb_hypgeom_pfq_sum_rs(s, t, a, p, b, q, z, n, prec);
    }
    else if (prec >= 1500 && n >= 30 + 100000 / (prec - 1000))
    {
        acb_hypgeom_pfq_sum_fme(s, t, a, p, b, q, z, n, prec);
    }
    else
    {
        /* With inexact complex numbers, binary splitting
           as it reduces the wrapping effect. */
        if (n > 8)
            acb_hypgeom_pfq_sum_bs_invz(s, t, a, p, b, q, zinv, n, prec);
        else
            acb_hypgeom_pfq_sum_forward(s, t, a, p, b, q, z, n, prec);
    }
}
