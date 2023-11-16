/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

/* Todo: could special-case quadratics, in particular Gaussian rationals */
void
qqbar_get_acb(acb_t res, const qqbar_t x, slong prec)
{
    if (qqbar_degree(x) == 1)
    {
        arb_set_fmpz(acb_realref(res), QQBAR_COEFFS(x));
        arb_div_fmpz(acb_realref(res), acb_realref(res), QQBAR_COEFFS(x) + 1, prec);
        arb_neg(acb_realref(res), acb_realref(res));
        arb_zero(acb_imagref(res));
    }
    else
    {
        arb_t t;
        fmpz_t n;
        slong wp;
        int imag_zero, real_zero;

        imag_zero = (qqbar_sgn_im(x) == 0);
        real_zero = (qqbar_sgn_re(x) == 0);

        acb_set(res, QQBAR_ENCLOSURE(x));

        for (wp = prec + 30; ; wp *= 2)
        {
            _qqbar_enclosure_raw(res, QQBAR_POLY(x), res, wp);

            if (imag_zero)
                arb_zero(acb_imagref(res));
            if (real_zero)
                arb_zero(acb_realref(res));

            if (arb_rel_accuracy_bits(acb_realref(res)) > prec + 5 &&
                arb_rel_accuracy_bits(acb_imagref(res)) > prec + 5)
                break;
        }

        /* Detect exact (dyadic real/imag parts) */
        arb_init(t);
        fmpz_init(n);

        arb_mul_2exp_si(t, acb_realref(res), wp);
        if (!arb_is_exact(t) && arb_get_unique_fmpz(n, t))
        {
            qqbar_t u;
            qqbar_init(u);
            qqbar_set_fmpz(u, n);
            qqbar_mul_2exp_si(u, u, wp);
            qqbar_sub(u, x, u);
            if (qqbar_sgn_re(u) == 0)
            {
                arb_set_fmpz(acb_realref(res), n);
                arb_mul_2exp_si(acb_realref(res), acb_realref(res), wp);
            }
            qqbar_clear(u);
        }

        arb_mul_2exp_si(t, acb_imagref(res), wp);
        if (!arb_is_exact(t) && arb_get_unique_fmpz(n, t))
        {
            qqbar_t u;
            qqbar_init(u);
            qqbar_i(u);
            qqbar_mul_fmpz(u, u, n);
            qqbar_mul_2exp_si(u, u, wp);
            qqbar_sub(u, x, u);
            if (qqbar_sgn_im(u) == 0)
            {
                arb_set_fmpz(acb_imagref(res), n);
                arb_mul_2exp_si(acb_imagref(res), acb_imagref(res), wp);
            }
            qqbar_clear(u);
        }

        acb_set_round(res, res, prec);

        arb_clear(t);
        fmpz_clear(n);
    }
}

