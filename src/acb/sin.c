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
acb_sin(acb_t r, const acb_t z, slong prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)
    if (arb_is_zero(b))
    {
        arb_sin(acb_realref(r), a, prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(a))
    {
        arb_sinh(acb_imagref(r), b, prec);
        arb_zero(acb_realref(r));
    }
    else
    {
        arb_t sa, ca, sb, cb;

        arb_init(sa);
        arb_init(ca);
        arb_init(sb);
        arb_init(cb);

        arb_sin_cos(sa, ca, a, prec);
        arb_sinh_cosh(sb, cb, b, prec);

        arb_mul(acb_realref(r), sa, cb, prec);
        arb_mul(acb_imagref(r), sb, ca, prec);

        arb_clear(sa);
        arb_clear(ca);
        arb_clear(sb);
        arb_clear(cb);
    }
#undef a
#undef b
}


