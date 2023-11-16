/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb.h"

int polylog_is_real(const acb_t s, const acb_t z);

void
arb_polylog(arb_t w, const arb_t s, const arb_t z, slong prec)
{
    acb_t ss, zz;
    acb_init(ss);
    acb_init(zz);
    acb_set_arb(ss, s);
    acb_set_arb(zz, z);
    if (polylog_is_real(ss, zz))
    {
        acb_polylog(zz, ss, zz, prec);
        arb_set(w, acb_realref(zz));
    }
    else
    {
        arb_indeterminate(w);
    }
    acb_clear(ss);
    acb_clear(zz);
}

void
arb_polylog_si(arb_t w, slong s, const arb_t z, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_si(t, s);
    arb_polylog(w, t, z, prec);
    arb_clear(t);
}
