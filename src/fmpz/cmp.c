/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "mpn_extras.h"
#include "fmpz.h"

static int flint_mpn_cmp2abs(mp_srcptr x, slong xn, mp_srcptr a, slong an)
{
    mp_limb_t xhi, ahi;

    FLINT_ASSERT(an >= 0);
    FLINT_ASSERT(xn >= 0);
    FLINT_ASSERT(xn == 0 || x[xn - 1] != 0);
    FLINT_ASSERT(an == 0 || a[an - 1] != 0);

    if (an > xn)
        return -1;

    if (an + 1 < xn)
        return 1;

    xhi = an == xn ? 0 : x[an];
    ahi = 0;

    while (an > 0)
    {
        ahi = MPN_LEFT_SHIFT_HI(ahi, a[an - 1], 1);
        if (xhi != ahi)
            return xhi > ahi ? 1 : -1;
        ahi = a[an - 1];
        xhi = x[an - 1];
        an--;
    }

    ahi = MPN_LEFT_SHIFT_HI(ahi, UWORD(0), 1);
    if (xhi != ahi)
        return xhi > ahi ? 1 : -1;

    return 0;
}

/* compare |a| and 2|b| */
int fmpz_cmp2abs(const fmpz_t a, const fmpz_t b)
{
    if (!COEFF_IS_MPZ(*b))
    {
        mp_limb_t ub = FLINT_ABS(*b);

        if (!COEFF_IS_MPZ(*a))
        {
            mp_limb_t ua = FLINT_ABS(*a);
            return ua < 2*ub ? -1 : ua > 2*ub ? 1 : 0;
        }
        else
        {
            return flint_mpn_cmp2abs(COEFF_TO_PTR(*a)->_mp_d,
                                     FLINT_ABS(COEFF_TO_PTR(*a)->_mp_size),
                                     &ub, ub != 0);
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(*a))
        {
            return -1;
        }
        else
        {
            return flint_mpn_cmp2abs(COEFF_TO_PTR(*a)->_mp_d,
                                     FLINT_ABS(COEFF_TO_PTR(*a)->_mp_size),
                                     COEFF_TO_PTR(*b)->_mp_d,
                                     FLINT_ABS(COEFF_TO_PTR(*b)->_mp_size));
        }
    }
}

int fmpz_cmpabs(const fmpz_t f, const fmpz_t g)
{
    if (f == g) return 0;  /* aliased inputs */

    if (!COEFF_IS_MPZ(*f))
    {
        if (!COEFF_IS_MPZ(*g))
        {
            mp_limb_t uf = FLINT_ABS(*f);
            mp_limb_t ug = FLINT_ABS(*g);

            return (uf < ug ? -1 : (uf > ug));
        }
        else return -1;
    }
    else
    {
        if (!COEFF_IS_MPZ(*g)) return 1;  /* f is large, so if g isn't... */
        else return mpz_cmpabs(COEFF_TO_PTR(*f), COEFF_TO_PTR(*g));
    }
}

int
fmpz_cmp(const fmpz_t f, const fmpz_t g)
{
    int sign;

    if (f == g)
        return 0;  /* aliased inputs */

    if (!COEFF_IS_MPZ(*f))
    {
        if (!COEFF_IS_MPZ(*g))
        {
            return (*f < *g ? -1 : *f > *g);
        }
        else  /* f is small, g is large */
        {
            sign = mpz_sgn(COEFF_TO_PTR(*g));
            return (*f >= 0 && sign < 0) ? 1 : -sign;
        }
    }
    else
    {
        if (!COEFF_IS_MPZ(*g))  /* f is large, and g is small */
        {
            sign = mpz_sgn(COEFF_TO_PTR(*f));
            return (*g >= 0 && sign < 0) ? -1 : sign;
        }
        else
            return mpz_cmp(COEFF_TO_PTR(*f), COEFF_TO_PTR(*g));
    }
}

int
fmpz_cmp_si(const fmpz_t f, slong g)
{
    fmpz c = *f;

    if (!COEFF_IS_MPZ(c))    /* f is small */
        return c < g ? -1 : c > g;
    else                     /* f is large */
        return flint_mpz_cmp_si(COEFF_TO_PTR(c), g);
}

int
fmpz_cmp_ui(const fmpz_t f, ulong g)
{
    fmpz c = *f;

    if (!COEFF_IS_MPZ(c))    /* f is small */
    {
        if (c < WORD(0) || g > COEFF_MAX)
            return -1;
        else
            return c < (slong) g ? -1 : c > (slong) g;
    }
    else                     /* f is large */
        return flint_mpz_cmp_ui(COEFF_TO_PTR(c), g);
}
