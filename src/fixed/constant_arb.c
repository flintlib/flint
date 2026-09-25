/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "arb.h"
#include "fixed.h"

/* res = f(x, n) as an arb ball accurate to prec bits: the fball
   constants give about FLINT_BITS (n - 1) bits, so two extra limbs
   cover the rounding to prec.  This is how the arb constants (which
   keep their own caches) call into this module. */
void
_fixed_constant_arb(arb_t res, fixed_constant_func f, slong prec)
{
    slong n = (prec + FLINT_BITS - 1) / FLINT_BITS + 2;
    fball_t x;

    fball_init(x);
    f(x, n);
    fball_get_arb(res, x);
    arb_set_round(res, res, prec);
    fball_clear(x);
}
