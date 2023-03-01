/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_sin_cos_pi(acb_t s, acb_t c, const acb_t z, slong prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)
    if (arb_is_zero(b))
    {
        arb_sin_cos_pi(acb_realref(s), acb_realref(c), a, prec);
        arb_zero(acb_imagref(s));
        arb_zero(acb_imagref(c));
    }
    else if (arb_is_zero(a))
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, prec);
        arb_mul(t, t, b, prec);
        arb_sinh_cosh(acb_imagref(s), acb_realref(c), t, prec);
        arb_zero(acb_realref(s));
        arb_zero(acb_imagref(c));
        arb_clear(t);
    }
    else
    {
        arb_t sa, ca, sb, cb;

        arb_init(sa);
        arb_init(ca);
        arb_init(sb);
        arb_init(cb);

        arb_const_pi(sb, prec);
        arb_mul(sb, sb, b, prec);

        arb_sin_cos_pi(sa, ca, a, prec);
        arb_sinh_cosh(sb, cb, sb, prec);

        arb_mul(acb_realref(s), sa, cb, prec);
        arb_mul(acb_imagref(s), sb, ca, prec);

        arb_mul(acb_realref(c), ca, cb, prec);
        arb_mul(acb_imagref(c), sa, sb, prec);
        arb_neg(acb_imagref(c), acb_imagref(c));

        arb_clear(sa);
        arb_clear(ca);
        arb_clear(sb);
        arb_clear(cb);
    }
#undef a
#undef b
}

