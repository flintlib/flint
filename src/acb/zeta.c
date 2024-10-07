/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_poly.h"
#include "acb_dirichlet.h"

/* FIXME: Remove this guard against warnings. Best thing would probably be to
 * implement an *-impl.h to keep track of local functions. */
#ifdef __GNUC__
# pragma GCC diagnostic ignored "-Wmissing-prototypes"
#endif

void
acb_zeta_si(acb_t z, slong s, slong prec)
{
    if (s >= 0)
    {
        arb_zeta_ui(acb_realref(z), s, prec);
    }
    else
    {
        arb_bernoulli_ui(acb_realref(z), 1-s, prec);
        arb_div_ui(acb_realref(z), acb_realref(z), 1-s, prec);
        arb_neg(acb_realref(z), acb_realref(z));
    }

    arb_zero(acb_imagref(z));
    return;
}

void
acb_hurwitz_zeta(acb_t z, const acb_t s, const acb_t a, slong prec)
{
    acb_dirichlet_hurwitz(z, s, a, prec);
}

void
acb_zeta(acb_t z, const acb_t s, slong prec)
{
    acb_dirichlet_zeta(z, s, prec);
}
