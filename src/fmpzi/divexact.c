/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpzi.h"

static void
_fmpzi_divexact(fmpzi_t q, const fmpzi_t x, const fmpzi_t y)
{
    fmpzi_t t, y_conj;
    fmpz_t v;
    mpz_t ytmp;

    fmpzi_init(t);
    fmpz_init(v);

    /* shallow conjugate */
    *fmpzi_realref(y_conj) = *fmpzi_realref(y);
    if (!COEFF_IS_MPZ(*fmpzi_imagref(y)))
    {
        *fmpzi_imagref(y_conj) = -*fmpzi_imagref(y);
    }
    else
    {
        *ytmp = *COEFF_TO_PTR(*fmpzi_imagref(y));
        mpz_neg(ytmp, ytmp);
        *fmpzi_imagref(y_conj) = PTR_TO_COEFF(ytmp);
    }

    fmpzi_mul(t, x, y_conj);

    fmpz_fmma(v, fmpzi_realref(y), fmpzi_realref(y),
                 fmpzi_imagref(y), fmpzi_imagref(y));

    fmpz_divexact(fmpzi_realref(q), fmpzi_realref(t), v);
    fmpz_divexact(fmpzi_imagref(q), fmpzi_imagref(t), v);

    fmpzi_clear(t);
    fmpz_clear(v);
}

void
fmpzi_divexact(fmpzi_t q, const fmpzi_t x, const fmpzi_t y)
{
    slong xbits, ybits, zbits, trunc;

    if (fmpz_is_zero(fmpzi_imagref(y)))
    {
        fmpz_divexact(fmpzi_imagref(q), fmpzi_imagref(x), fmpzi_realref(y));
        fmpz_divexact(fmpzi_realref(q), fmpzi_realref(x), fmpzi_realref(y));
        return;
    }

    if (fmpz_is_zero(fmpzi_realref(y)))
    {
        fmpz_divexact(fmpzi_realref(q), fmpzi_realref(x), fmpzi_imagref(y));
        fmpz_divexact(fmpzi_imagref(q), fmpzi_imagref(x), fmpzi_imagref(y));
        fmpzi_div_i(q, q);
        return;
    }

    xbits = fmpzi_bits(x);
    ybits = fmpzi_bits(y);

    if (ybits == 0)
    {
        flint_throw(FLINT_ERROR, "fmpzi_divexact: division by zero\n");
    }

    if (xbits == 0)
    {
        fmpzi_zero(q);
        return;
    }

    /* todo: special cases? */
    if (ybits == 1)
    {
    }

    zbits = xbits - ybits;

    if (zbits < 45)
    {
        double a, b, c, d, t, u, v, w, qa, qb;
        slong aexp, bexp, cexp, dexp;

        if (xbits < 500)
        {
            a = fmpz_get_d(fmpzi_realref(x));
            b = fmpz_get_d(fmpzi_imagref(x));
            c = fmpz_get_d(fmpzi_realref(y));
            d = fmpz_get_d(fmpzi_imagref(y));
        }
        else
        {
            a = fmpz_get_d_2exp(&aexp, fmpzi_realref(x));
            b = fmpz_get_d_2exp(&bexp, fmpzi_imagref(x));
            c = fmpz_get_d_2exp(&cexp, fmpzi_realref(y));
            d = fmpz_get_d_2exp(&dexp, fmpzi_imagref(y));

            a = ldexp(a, FLINT_MAX(aexp - xbits, -1024));
            b = ldexp(b, FLINT_MAX(bexp - xbits, -1024));
            c = ldexp(c, FLINT_MAX(cexp - xbits, -1024));
            d = ldexp(d, FLINT_MAX(dexp - xbits, -1024));
        }

        t = a * c + b * d;
        u = b * c - a * d;
        v = c * c + d * d;

        w = 0.5 / v;

        t = (2.0 * t + v) * w;
        u = (2.0 * u + v) * w;

        qa = floor(t);
        qb = floor(u);

        fmpz_set_d(fmpzi_realref(q), qa);
        fmpz_set_d(fmpzi_imagref(q), qb);

        return;
    }

    if (ybits > zbits * 1.25 + 256)
    {
        fmpzi_t t, u;

        fmpzi_init(t);
        fmpzi_init(u);

        trunc = ybits - zbits - 20;

        fmpz_tdiv_q_2exp(fmpzi_realref(t), fmpzi_realref(x), trunc);
        fmpz_tdiv_q_2exp(fmpzi_imagref(t), fmpzi_imagref(x), trunc);
        fmpz_tdiv_q_2exp(fmpzi_realref(u), fmpzi_realref(y), trunc);
        fmpz_tdiv_q_2exp(fmpzi_imagref(u), fmpzi_imagref(y), trunc);

        fmpzi_divrem_approx(q, NULL, t, u);

        fmpzi_clear(t);
        fmpzi_clear(u);
    }
    else
    {
        _fmpzi_divexact(q, x, y);
    }
}
