/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <math.h>
#include "ulong_extras.h"

static int n_is_square_and_get_sqrt(ulong * s, ulong x)
{
    ulong sq = sqrt((double) x) + 0.5;

    *s = sq;
    return x == sq * sq;
}


#define r_shift(in, c) (((c) == FLINT_BITS) ? WORD(0) : ((in) >> (c)))

static ulong _ll_factor_SQUFOF(ulong n_hi, ulong n_lo, ulong max_iters)
{
    ulong n[2];
    ulong nsqrt[2];
    ulong rem[2];
    slong num, sqroot;

    ulong p, q;

    ulong l, l2, iq, pnext;
    ulong qarr[50];
    ulong qupto, qlast, t, r = 0;
    ulong i, j;

    n[0] = n_lo;
    n[1] = n_hi;

    if (n_hi) num = mpn_sqrtrem(nsqrt, rem, n, 2);
    else num = ((nsqrt[0] = n_sqrtrem(rem, n_lo)) != UWORD(0));

    sqroot = nsqrt[0];
    p = sqroot;
    q = rem[0];

    if ((q == 0) || (num == 0))
    {
        return sqroot;
    }

    l = 1 + 2*n_sqrt(2*p);
    l2 = l/2;
    qupto = 0;
    qlast = 1;

    for (i = 0; i < max_iters; i++)
    {
        iq = (sqroot + p)/q;
        pnext = iq*q - p;
        if (q <= l)
        {
            if ((q & UWORD(1)) == UWORD(0))
            {
                qarr[qupto] = q/2;
                qupto++;
                if (qupto >= UWORD(50)) return UWORD(0);
            } else if (q <= l2)
            {
                qarr[qupto] = q;
                qupto++;
                if (qupto >= UWORD(50)) return UWORD(0);
            }
        }

        t = qlast + iq*(p - pnext);
        qlast = q;
        q = t;
        p = pnext;
        if ((i & 1) == 1) continue;
        if (!n_is_square_and_get_sqrt(&r, q))
            continue;
        if (qupto == UWORD(0)) break;
        for (j = 0; j < qupto; j++)
            if (r == qarr[j]) goto cont;
        break;
cont: ;
      if (r == UWORD(1)) return UWORD(0);
    }

    if (i == max_iters) return UWORD(0);  /* taken too much time, give up */

    qlast = r;
    p = p + r*((sqroot - p)/r);

    umul_ppmm(rem[1], rem[0], p, p);
    sub_ddmmss(nsqrt[1], nsqrt[0], n[1], n[0], rem[1], rem[0]);
    if (nsqrt[1])
    {
        int norm;
        norm = flint_clz(qlast);
        udiv_qrnnd(q, rem[0], (nsqrt[1] << norm) + r_shift(nsqrt[0], FLINT_BITS - norm), nsqrt[0] << norm, qlast << norm);
        rem[0] >>= norm;
    }
    else
    {
        q = nsqrt[0]/qlast;
    }

    for (j = 0; j < max_iters; j++)
    {
        iq = (sqroot + p)/q;
        pnext = iq*q - p;
        if (p == pnext) break;
        t = qlast + iq*(p - pnext);
        qlast = q;
        q = t;
        p = pnext;
    }

    if (j == max_iters) return UWORD(0);  /* taken too much time, give up */

    if ((q & UWORD(1)) == UWORD(0)) q /= UWORD(2);

    return q;
}

/* Make sure multiplier does not overflow */
#define MAX_MULTIPLIER_BITS 16

ulong n_ll_factor_SQUFOF(ulong nhi, ulong nlo, ulong iters)
{
    ulong factor = _ll_factor_SQUFOF(nhi, nlo, iters);
    ulong multiplier;
    ulong quot, rem;
    ulong i;

    if (nhi >= UWORD(1) << (FLINT_BITS - MAX_MULTIPLIER_BITS))
        return 0;

    for (i = 1; (i < FLINT_NUM_PRIMES_SMALL) && !factor; i++)
    {
        ulong multn[2];
        multiplier = flint_primes_small[i];
        FLINT_ASSERT(multiplier < (UWORD(1) << MAX_MULTIPLIER_BITS));

        umul_ppmm(multn[1], multn[0], multiplier, nlo);
        multn[1] += multiplier * nhi;
        factor = _ll_factor_SQUFOF(multn[1], multn[0], iters);

        if (factor)
        {
            quot = factor/multiplier;
            rem = factor - quot*multiplier;
            if (!rem) factor = quot;
            /* The factor is trivial */
            if ((factor == UWORD(1)) || (factor == nlo && nhi == 0))
                factor = UWORD(0);
        }
    }

    if (i == FLINT_NUM_PRIMES_SMALL) return UWORD(0);

    return factor;
}

ulong n_factor_SQUFOF(ulong n, ulong iters)
{
    return n_ll_factor_SQUFOF(0, n, iters);
}

