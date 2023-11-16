/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec);

void
acb_hypgeom_gamma_stirling_sum_horner(acb_t s, const acb_t z, slong N, slong prec)
{
    acb_t b, t, zinv, w;
    mag_t zinv_mag;
    slong n, term_mag, term_prec;
    slong * term_mags;

    if (N <= 1)
    {
        acb_zero(s);
        return;
    }

    acb_init(b);
    acb_init(t);
    acb_init(zinv);
    acb_init(w);
    mag_init(zinv_mag);

    acb_inv(zinv, z, prec);
    acb_mul(w, zinv, zinv, prec);

    acb_get_mag(zinv_mag, zinv);
    term_mags = flint_malloc(sizeof(ulong) * N);

    _arb_hypgeom_gamma_stirling_term_bounds(term_mags, zinv_mag, N);

    acb_zero(s);

    for (n = N - 1; n >= 1; n--)
    {
        term_mag = term_mags[n];
        term_prec = prec + term_mag;
        term_prec = FLINT_MIN(term_prec, prec);
        term_prec = FLINT_MAX(term_prec, 10);

        if (prec - term_prec > 200)
        {
            acb_set_round(t, w, term_prec);
            acb_mul(s, s, t, term_prec);
        }
        else
            acb_mul(s, s, w, term_prec);

        arb_gamma_stirling_coeff(acb_realref(b), n, 0, term_prec);
        acb_add_arb(s, s, acb_realref(b), term_prec);
    }

    acb_mul(s, s, zinv, prec);

    flint_free(term_mags);

    acb_clear(t);
    acb_clear(b);
    acb_clear(zinv);
    acb_clear(w);
    mag_clear(zinv_mag);
}

