/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "arb.h"

void
arb_get_interval_mpfr(mpfr_t a, mpfr_t b, const arb_t x)
{
    arf_t r, t;
    arf_init(t);
    arf_init_set_mag_shallow(r, arb_radref(x));
    arf_sub(t, arb_midref(x), r, mpfr_get_prec(a), ARF_RND_FLOOR);
    arf_get_mpfr(a, t, MPFR_RNDD);
    arf_add(t, arb_midref(x), r, mpfr_get_prec(b), ARF_RND_CEIL);
    arf_get_mpfr(b, t, MPFR_RNDU);
    arf_clear(t);
}

