/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

void arb_gamma_stirling_coeff(arb_t b, ulong k, int digamma, slong prec);

void
arb_hypgeom_gamma_stirling_sum_horner(arb_t s, const arb_t z, slong N, slong prec)
{
    arb_t b, t, zinv, w;
    mag_t zinv_mag;
    slong n, term_mag, term_prec;
    slong * term_mags;

    if (N <= 1)
    {
        arb_zero(s);
        return;
    }

    arb_init(b);
    arb_init(t);
    arb_init(zinv);
    arb_init(w);
    mag_init(zinv_mag);

    arb_inv(zinv, z, prec);
    arb_mul(w, zinv, zinv, prec);

    arb_get_mag(zinv_mag, zinv);
    term_mags = flint_malloc(sizeof(ulong) * N);

    _arb_hypgeom_gamma_stirling_term_bounds(term_mags, zinv_mag, N);

    arb_zero(s);

    for (n = N - 1; n >= 1; n--)
    {
        term_mag = term_mags[n];
        term_prec = prec + term_mag;
        term_prec = FLINT_MIN(term_prec, prec);
        term_prec = FLINT_MAX(term_prec, 10);

        if (prec - term_prec > 200)
        {
            arb_set_round(t, w, term_prec);
            arb_mul(s, s, t, term_prec);
        }
        else
            arb_mul(s, s, w, term_prec);

        arb_gamma_stirling_coeff(b, n, 0, term_prec);
        arb_add(s, s, b, term_prec);
    }

    arb_mul(s, s, zinv, prec);

    flint_free(term_mags);

    arb_clear(t);
    arb_clear(b);
    arb_clear(zinv);
    arb_clear(w);
    mag_clear(zinv_mag);
}

