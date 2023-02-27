/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void 
_fmpz_poly_evaluate_divconquer_fmpq(fmpz_t rnum, fmpz_t rden,
        const fmpz * poly, slong len, const fmpz_t xnum, const fmpz_t xden)
{
    slong c, h, i, k = 1;
    fmpz *ynum, *yden, *Tnum, *Tden, *tnum = rnum, *tden = rden, *unum, *uden;
    fmpz_t d;

    h = FLINT_BIT_COUNT(len - 1);  /* 2^{h-1} < len <= 2^h */
    ynum = _fmpz_vec_init(2 * h + 2); /* x^{2^0}, x^{2^1}, ..., x^{2^{h-1}} */
    yden = _fmpz_vec_init(2 * h + 2); /* x^{2^0}, x^{2^1}, ..., x^{2^{h-1}} */
    fmpz_init(d);

    Tnum = ynum + h;
    Tden = yden + h;

    unum = ynum + 2 * h + 1;
    uden = yden + 2 * h + 1;

    *ynum = *xnum;
    *yden = *xden;

    for (i = 1; i < h; i++)
    {
        fmpz_mul(ynum + i, ynum + (i - 1), ynum + (i - 1));
        fmpz_mul(yden + i, yden + (i - 1), yden + (i - 1));
    }

    for (i = 0; i < len - 1; )
    {
        /* t = poly[i] + y[0] * poly[i+1] */
        fmpz_mul(tnum, ynum + 0, poly + i + 1);
        fmpz_addmul(tnum, yden + 0, poly + i);
        fmpz_set(tden, yden + 0);

        i += 2;
        count_trailing_zeros(c, i);

        for (k = 1; k < c; k++)
        {
            /* t = T[k] + y[k] * t */
            fmpz_mul(unum, ynum + k, tnum);
            fmpz_mul(uden, yden + k, tden);
            fmpz_mul(tnum, unum, Tden + k);
            fmpz_addmul(tnum, uden, Tnum + k);
            fmpz_mul(tden, Tden + k, uden);
        }

        fmpz_swap(Tnum + k, tnum);
        fmpz_swap(Tden + k, tden);
    }

    if (len & WORD(1))
    {
        fmpz_set(tnum, poly + (len - 1));
        fmpz_one(tden);

        count_trailing_zeros(c, len + 1);

        for (k = 1; k < c; k++)
        {
            fmpz_mul(unum, ynum + k, tnum);
            fmpz_mul(uden, yden + k, tden);
            fmpz_mul(tnum, unum, Tden + k);
            fmpz_addmul(tnum, uden, Tnum + k);
            fmpz_mul(tden, Tden + k, uden);
        }

        fmpz_swap(Tnum + k, tnum);
        fmpz_swap(Tden + k, tden);
    }

    fmpz_swap(tnum, Tnum + k);
    fmpz_swap(tden, Tden + k);

    for ( ; k < h; k++)
    {
        if ((len - 1) & (WORD(1) << k))
        {
            fmpz_mul(unum, ynum + k, tnum);
            fmpz_mul(uden, yden + k, tden);
            fmpz_mul(tnum, unum, Tden + k);
            fmpz_addmul(tnum, uden, Tnum + k);
            fmpz_mul(tden, Tden + k, uden);
        }
    }

    fmpz_gcd(d, rnum, rden);
    fmpz_divexact(rnum, rnum, d);
    fmpz_divexact(rden, rden, d);

    *ynum = WORD(0);
    *yden = WORD(0);
    _fmpz_vec_clear(ynum, 2 * h + 2);
    _fmpz_vec_clear(yden, 2 * h + 2);
    fmpz_clear(d);
}

void
fmpz_poly_evaluate_divconquer_fmpq(fmpq_t res, const fmpz_poly_t f, 
                                   const fmpq_t a)
{
    if (fmpz_poly_is_zero(f))
    {
        fmpq_zero(res);
        return;
    }

    if (res == a)
    {
        fmpq_t t;
        fmpq_init(t);
        fmpz_poly_evaluate_divconquer_fmpq(t, f, a);
        fmpq_swap(res, t);
        fmpq_clear(t);
    }
    else
    {
        _fmpz_poly_evaluate_divconquer_fmpq(fmpq_numref(res), fmpq_denref(res),
            f->coeffs, f->length, fmpq_numref(a), fmpq_denref(a));
    }
}
