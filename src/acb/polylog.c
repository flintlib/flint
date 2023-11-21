/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_poly.h"

void
acb_polylog(acb_t w, const acb_t s, const acb_t z, slong prec)
{
    acb_t t;
    acb_init(t);
    _acb_poly_polylog_cpx(t, s, z, 1, prec);
    acb_swap(w, t);
    acb_clear(t);
}

void
acb_polylog_si(acb_t w, slong s, const acb_t z, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_set_si(t, s);
    acb_polylog(w, t, z, prec);
    acb_clear(t);
}
