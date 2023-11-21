/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_exp_invexp(acb_t r, acb_t s, const acb_t z, slong prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)
    if (arb_is_zero(b))
    {
        arb_exp_invexp(acb_realref(r), acb_realref(s), a, prec);
        arb_zero(acb_imagref(r));
        arb_zero(acb_imagref(s));
    }
    else if (arb_is_zero(a))
    {
        arb_sin_cos(acb_imagref(r), acb_realref(r), b, prec);
        acb_conj(s, r);
    }
    else
    {
        arb_t t, u, v, w;

        arb_init(t);
        arb_init(u);
        arb_init(v);
        arb_init(w);

        arb_exp_invexp(t, u, a, prec);
        arb_sin_cos(v, w, b, prec);

        arb_mul(acb_realref(r), t, w, prec);
        arb_mul(acb_imagref(r), t, v, prec);

        arb_mul(acb_realref(s), u, w, prec);
        arb_mul(acb_imagref(s), u, v, prec);
        arb_neg(acb_imagref(s), acb_imagref(s));

        arb_clear(t);
        arb_clear(u);
        arb_clear(v);
        arb_clear(w);
    }
#undef a
#undef b
}

