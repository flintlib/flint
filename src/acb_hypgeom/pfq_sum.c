/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_pfq_sum(acb_t s, acb_t t, acb_srcptr a, slong p,
    acb_srcptr b, slong q, const acb_t z, slong n, slong prec)
{
    if (n > 4 && prec >= 128
        && _acb_vec_bits(a, p) * p + _acb_vec_bits(b, q) * q + 10 < prec / 2)
    {
        if (prec >= 256 && acb_bits(z) < prec * 0.01)
            acb_hypgeom_pfq_sum_bs(s, t, a, p, b, q, z, n, prec);
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

