/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpzi.h"

void
fmpzi_divrem(fmpzi_t q, fmpzi_t r, const fmpzi_t x, const fmpzi_t y)
{
    fmpzi_t t, y_conj;
    fmpz_t v;
    mpz_t ytmp;
    slong xbits, ybits;

    xbits = fmpzi_bits(x);
    ybits = fmpzi_bits(y);

    if (ybits == 0)
    {
        flint_throw(FLINT_ERROR, "fmpzi_divrem: division by zero\n");
    }

    if (xbits == 0)
    {
        fmpzi_zero(q);
        if (r != NULL)
            fmpzi_zero(r);
        return;
    }

    if (xbits < ybits - 2)
    {
        if (r != NULL)
            fmpzi_set(r, x);
        fmpzi_zero(q);
        return;
    }

    if (q == x || q == y)
    {
        fmpzi_init(t);
        fmpzi_divrem(t, r, x, y);
        fmpzi_swap(q, t);
        fmpzi_clear(t);
        return;
    }

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

    fmpz_mul_2exp(fmpzi_realref(t), fmpzi_realref(t), 1);
    fmpz_mul_2exp(fmpzi_imagref(t), fmpzi_imagref(t), 1);

    fmpz_fmma(v, fmpzi_realref(y), fmpzi_realref(y),
                 fmpzi_imagref(y), fmpzi_imagref(y));

    fmpz_add(fmpzi_realref(t), fmpzi_realref(t), v);
    fmpz_add(fmpzi_imagref(t), fmpzi_imagref(t), v);

    fmpz_mul_2exp(v, v, 1);

    fmpz_fdiv_q(fmpzi_realref(q), fmpzi_realref(t), v);
    fmpz_fdiv_q(fmpzi_imagref(q), fmpzi_imagref(t), v);

    if (r != NULL)
    {
        fmpzi_mul(t, q, y);
        fmpzi_sub(r, x, t);
    }

    fmpzi_clear(t);
    fmpz_clear(v);
}
