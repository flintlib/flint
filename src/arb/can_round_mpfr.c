/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "arb.h"

int mpfr_round_p(mp_srcptr, mp_size_t, mpfr_exp_t, mpfr_prec_t);

int
arb_can_round_arf(const arb_t x, slong prec, arf_rnd_t rnd)
{
    return arb_can_round_mpfr(x, prec, arf_rnd_to_mpfr(rnd));
}

int
arb_can_round_mpfr(const arb_t x, slong prec, mpfr_rnd_t rnd)
{
    if (!arb_is_finite(x))
    {
        return 0;
    }
    else if (mag_is_zero(arb_radref(x)))
    {
        return 1;
    }
    else if (arf_is_zero(arb_midref(x)))
    {
        return 0;
    }
    else
    {
        slong e, bits;
        mp_size_t n;
        mp_srcptr d;

        e = _fmpz_sub_small(ARF_EXPREF(arb_midref(x)), MAG_EXPREF(arb_radref(x)));

        if (e < prec)
            return 0;

        /* The relative exponent could be tiny (in which case _fmpz_sub_small
           has clamped it). Looking just past the end will be enough. */
        bits = arb_bits(x);
        e = FLINT_MIN(e, FLINT_MAX(bits, prec) + 10);

        ARF_GET_MPN_READONLY(d, n, arb_midref(x));

        return mpfr_round_p(d, n, e, prec + (rnd == MPFR_RNDN));
    }
}

