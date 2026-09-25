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
#include "mp_real.h"
#include "impl.h"

/* res = f(x, n) as an arb ball accurate to prec bits: the mp_real
   constants give about FLINT_BITS (n - 1) bits, so two extra limbs
   cover the rounding to prec.  This is how the arb constants (which
   keep their own caches) call into this module. */
void
_mp_real_const_arb(arb_t res, mp_real_const_func f, slong prec)
{
    slong n = (prec + FLINT_BITS - 1) / FLINT_BITS + 2;
    mp_real_t x;

    mp_real_init(x);
    f(x, n, 0);
    mp_real_get_arb(res, x);
    arb_set_round(res, res, prec);
    mp_real_clear(x);
}
