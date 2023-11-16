/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_euler_number_fmpz(arb_t res, const fmpz_t n, slong prec)
{
    if (fmpz_cmp_ui(n, UWORD_MAX) <= 0)
    {
        if (fmpz_sgn(n) >= 0)
            arb_euler_number_ui(res, fmpz_get_ui(n), prec);
        else
            arb_zero(res);
    }
    else if (fmpz_is_odd(n))
    {
        arb_zero(res);
    }
    else
    {
        arb_t t;
        fmpz_t m;
        slong wp;

        arb_init(t);
        fmpz_init(m);
        wp = prec + 2 * fmpz_bits(n);

        /* beta(n) ~= 1 */
        arf_one(arb_midref(res));
        mag_one(arb_radref(res));
        mag_mul_2exp_si(arb_radref(res), arb_radref(res), WORD_MIN);

        /* |E_n| = 2 n! beta(n+1) / (pi/2)^(n+1) */
        fmpz_add_ui(m, n, 1);
        arb_gamma_fmpz(t, m, wp);
        arb_mul(res, res, t, wp);

        arb_const_pi(t, wp);
        arb_mul_2exp_si(t, t, -1);
        arb_pow_fmpz(t, t, m, wp);

        arb_div(res, res, t, prec);
        arb_mul_2exp_si(res, res, 1);

        if (fmpz_fdiv_ui(n, 4) == 2)
            arb_neg(res, res);

        arb_clear(t);
        fmpz_clear(m);
    }
}

