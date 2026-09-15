/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

/* b aliases neither q nor r; rounds to nearest, ties to even */
static void _fmpz_ndiv_qr(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    int c, rbsgn;

    fmpz_tdiv_qr(q, r, a, b);

    c = fmpz_cmp2abs(b, r);    /* compare |b| with 2 |r| */

    if (c > 0)
        return;

    /* |2r| >= |b|: round away from zero unless it is a tie with q even */
    if (c == 0 && fmpz_is_even(q))
        return;

    rbsgn = fmpz_sgn(r) * fmpz_sgn(b);   /* sign of a / b */

    if (rbsgn < 0)
    {
        fmpz_sub_ui(q, q, 1);
        fmpz_add(r, r, b);
    }
    else
    {
        fmpz_add_ui(q, q, 1);
        fmpz_sub(r, r, b);
    }
}

void
fmpz_ndiv_qr(fmpz_t q, fmpz_t r, const fmpz_t a, const fmpz_t b)
{
    slong A = *a;
    slong B = *b;

    if (fmpz_is_zero(b))
    {
        flint_throw(FLINT_DIVZERO, "Exception: division by zero in fmpz_ndiv_qr\n");
    }

    if (!COEFF_IS_MPZ(A) && !COEFF_IS_MPZ(B))
    {
        slong lquo, lrem;

        _fmpz_demote(q);
        _fmpz_demote(r);

        if (FLINT_ABS(*b) == 1) /* avoid overflow in case */
        {                       /* a = 2^(SMALL_FMPZ_BITCOUNT_MAX) */
            fmpz_set_si(q, A * FLINT_SGN(B));
            fmpz_zero(r);
            return;
        }

        *q = A / B;
        *r = A - B * *q;
        lquo = *q + FLINT_SGN(A) * FLINT_SGN(B);
        lrem = A - B * lquo;

        /* round away from zero if that gives a smaller remainder, or on a
           tie if the truncated quotient is odd (ties to even) */
        if (FLINT_ABS(lrem) < FLINT_ABS(*r)
            || (FLINT_ABS(lrem) == FLINT_ABS(*r) && (*q & 1)))
        {
            *q = lquo;
            *r = lrem;
        }
    }
    else
    {
        if (b == q)
        {
            fmpz_t t;
            fmpz_init(t);
            _fmpz_ndiv_qr(t, r, a, b);
            fmpz_swap(q, t);
            fmpz_clear(t);
        }
        else if (b == r)
        {
            fmpz_t t;
            fmpz_init(t);
            _fmpz_ndiv_qr(q, t, a, b);
            fmpz_swap(r, t);
            fmpz_clear(t);
        }
        else
        {
            _fmpz_ndiv_qr(q, r, a, b);
        }
    }
}
