/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"

void
_acb_elliptic_rg(acb_t res, const acb_t x, const acb_t y, const acb_t z,
                    int flags, slong prec)
{
    acb_t a, b, c, t;
    slong wp;

    wp = prec + 10;

    acb_init(a);
    acb_init(b);
    acb_init(c);
    acb_init(t);

    acb_elliptic_rf(a, x, y, z, 0, wp);
    acb_mul(a, a, z, wp);

    acb_elliptic_rj(b, x, y, z, z, 0, wp);
    acb_sub(c, x, z, wp);
    acb_mul(b, b, c, wp);
    acb_sub(c, z, y, wp);
    acb_mul(b, b, c, wp);
    acb_div_ui(b, b, 3, wp);

    acb_sqrt(c, x, wp);
    acb_sqrt(t, y, wp);
    acb_mul(c, c, t, wp);
    acb_rsqrt(t, z, wp);
    acb_mul(c, c, t, wp);

    acb_add(res, a, b, wp);
    acb_add(res, res, c, prec);
    acb_mul_2exp_si(res, res, -1);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    acb_clear(t);
}

void
acb_elliptic_rg(acb_t res, const acb_t x, const acb_t y, const acb_t z,
                    int flags, slong prec)
{
    if (acb_is_zero(x) && acb_is_zero(y))
    {
        acb_sqrt(res, z, prec);
        acb_mul_2exp_si(res, res, -1);
    }
    else if (acb_is_zero(x) && acb_is_zero(z))
    {
        acb_sqrt(res, y, prec);
        acb_mul_2exp_si(res, res, -1);
    }
    else if (acb_is_zero(y) && acb_is_zero(z))
    {
        acb_sqrt(res, x, prec);
        acb_mul_2exp_si(res, res, -1);
    }
    else if (acb_contains_zero(z))
    {
        if (acb_contains_zero(y))
            _acb_elliptic_rg(res, y, z, x, flags, prec);
        else
            _acb_elliptic_rg(res, x, z, y, flags, prec);
    }
    else
    {
        _acb_elliptic_rg(res, x, y, z, flags, prec);
    }
}

