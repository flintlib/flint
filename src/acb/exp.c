/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb/impl.h"

void
acb_exp(acb_t r, const acb_t z, slong prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)
    if (arb_is_zero(b))
    {
        arb_exp(acb_realref(r), a, prec);
        arb_zero(acb_imagref(r));
    }
    else if (arb_is_zero(a))
    {
        arb_sin_cos(acb_imagref(r), acb_realref(r), b, prec);
    }
    else if (_acb_elem_small(r, NULL, z, ACB_ELEM_EXP, prec))
    {
        /* (the register paths) */
    }
    else
    {
        /* exp(a) (cos b + i sin b) with one temporary: the exponential
           first (a is then no longer needed), the sine and cosine into
           the parts of r (arb_sin_cos handles the aliasing of b) */
        arb_t t;

        arb_init(t);
        arb_exp(t, a, prec);
        arb_sin_cos(acb_imagref(r), acb_realref(r), b, prec);
        arb_mul(acb_realref(r), acb_realref(r), t, prec);
        arb_mul(acb_imagref(r), acb_imagref(r), t, prec);
        arb_clear(t);
    }
#undef a
#undef b
}
