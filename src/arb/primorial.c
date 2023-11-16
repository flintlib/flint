/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "arb.h"

#define NUM_BASECASE 30

static int
basecase(arb_t res, n_primes_t primes, ulong a, ulong b, ulong nmax, slong prec)
{
    ulong n, p, pp;
    mp_limb_t prod[NUM_BASECASE];
    mp_limb_t top;
    mp_size_t nlimbs;
    mp_limb_t hi, lo;
    int inexact, more;
    slong shift;

    nlimbs = 0;
    pp = 1;
    more = 1;

    for (n = a; n < b; n++)
    {
        p = n_primes_next(primes);

        if (p > nmax)
        {
            more = 0;
            break;
        }

        umul_ppmm(hi, lo, pp, p);

        if (hi != 0)
        {
            if (nlimbs == 0)
            {
                prod[0] = lo;
                prod[1] = hi;
                pp = 1;
                nlimbs = 2;
            }
            else
            {
                prod[nlimbs] = top = mpn_mul_1(prod, prod, nlimbs, pp);
                nlimbs += (top != 0);
                pp = p;
            }
        }
        else
        {
            pp = lo;
        }
    }

    if (nlimbs == 0)
    {
        arb_set_ui(res, pp);
        arb_set_round(res, res, prec);
        return more;
    }

    if (pp != 1)
    {
        prod[nlimbs] = top = mpn_mul_1(prod, prod, nlimbs, pp);
        nlimbs += (top != 0);
    }

    inexact = _arf_set_round_mpn(arb_midref(res), &shift, prod, nlimbs, 0, prec, ARB_RND);
    fmpz_set_si(ARF_EXPREF(arb_midref(res)), nlimbs * FLINT_BITS + shift);

    if (inexact)
        arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
    else
        mag_zero(arb_radref(res));

    return more;
}

static int
bsplit(arb_t res, n_primes_t primes, ulong a, ulong b, ulong nmax, slong prec)
{
    if (b - a < NUM_BASECASE)
    {
        return basecase(res, primes, a, b, nmax, prec);
    }
    else
    {
        int more;
        more = bsplit(res, primes, a, a + (b - a) / 2, nmax, prec + 3);

        if (more)
        {
            arb_t t;
            arb_init(t);
            more = bsplit(t, primes, a + (b - a) / 2, b, nmax, prec + 3);
            arb_mul(res, res, t, prec);
            arb_clear(t);
        }
        else
        {
            arb_set_round(res, res, prec);
        }

        return more;
    }
}

void
arb_primorial_nth_ui(arb_t res, ulong n, slong prec)
{
    if (n < 10)
    {
        const unsigned int tab[] = { 1, 2, 6, 30, 210,
            2310, 30030, 510510, 9699690, 223092870 };

        arb_set_ui(res, tab[n]);
        arb_set_round(res, res, prec);
    }
    else if (FLINT_BITS == 32 && n >= 203280220)
    {
        arb_indeterminate(res);    /* p_n will not fit a ulong */
    }
    else
    {
        n_primes_t primes;
        n_primes_init(primes);
        bsplit(res, primes, 0, n, UWORD_MAX, prec);
        n_primes_clear(primes);
    }
}

void
arb_primorial_ui(arb_t res, ulong n, slong prec)
{
    if (n < 17)
    {
        const unsigned short tab[] = { 1, 1, 2, 6, 6, 30, 30,
            210, 210, 210, 210, 2310, 2310, 30030, 30030, 30030, 30030 };

        arb_set_ui(res, tab[n]);
        arb_set_round(res, res, prec);
    }
    else if (n >= UWORD_MAX / 2)
    {
        /* avoid potential slong/ulong issues and overflow in the multiply by
           two below. in any case, this would be quite big on 64-bit... */
        arb_indeterminate(res);
    }
    else
    {
        n_primes_t primes;
        n_primes_init(primes);
        /* crude upper bound for pi(n) just to balance the binary splitting;
           we stop the prime iterator at the actual correct count */
        bsplit(res, primes, 0, 2 * n / FLINT_BIT_COUNT(n) + 1, n, prec);
        n_primes_clear(primes);
    }
}

