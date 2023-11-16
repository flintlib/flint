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
arb_set_interval_mpfr(arb_t x, const mpfr_t a, const mpfr_t b, slong prec)
{
    arf_t aa, bb;

    arf_init(aa);
    arf_init(bb);

    arf_set_mpfr(aa, a);
    arf_set_mpfr(bb, b);

    arb_set_interval_arf(x, aa, bb, prec);

    arf_clear(aa);
    arf_clear(bb);
}

