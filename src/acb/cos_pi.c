/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_cos_pi(acb_t r, const acb_t z, slong prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)
    if (arb_is_zero(b))
    {
        arb_cos_pi(acb_realref(r), a, prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(a))
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, prec);
        arb_mul(t, t, b, prec);
        arb_cosh(acb_realref(r), t, prec);
        arb_zero(acb_imagref(r));
        arb_clear(t);
    }
    else
    {
        arb_t sa, ca, sb, cb;

        arb_init(sa);
        arb_init(ca);
        arb_init(sb);
        arb_init(cb);

        arb_sin_cos_pi(sa, ca, a, prec);
        arb_const_pi(cb, prec);
        arb_mul(cb, cb, b, prec);
        arb_sinh_cosh(sb, cb, cb, prec);

        arb_mul(acb_realref(r), ca, cb, prec);
        arb_mul(acb_imagref(r), sa, sb, prec);
        arb_neg(acb_imagref(r), acb_imagref(r));

        arb_clear(sa);
        arb_clear(ca);
        arb_clear(sb);
        arb_clear(cb);
    }
#undef a
#undef b
}

