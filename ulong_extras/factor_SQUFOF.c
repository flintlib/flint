/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t _ll_factor_SQUFOF(mp_limb_t n_hi, mp_limb_t n_lo, ulong max_iters)
{
    mp_limb_t n[2];
	 mp_limb_t sqrt[2];
	 mp_limb_t rem[2];
	 slong num, sqroot, p, q;

    mp_limb_t l, l2, iq, pnext;
    mp_limb_t qarr[50];
    mp_limb_t qupto, qlast, t, r = 0;
    ulong i, j;

	 n[0] = n_lo;
	 n[1] = n_hi;

    if (n_hi) num = mpn_sqrtrem(sqrt, rem, n, 2);
    else num = ((sqrt[0] = n_sqrtrem(rem, n_lo)) != UWORD(0));
	
    sqroot = sqrt[0];
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
        if (!n_is_square(q)) continue;
        r = n_sqrt(q);
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
    sub_ddmmss(sqrt[1], sqrt[0], n[1], n[0], rem[1], rem[0]);
	if (sqrt[1])
	{
        int norm;
        count_leading_zeros(norm, qlast);
        udiv_qrnnd(q, rem[0], (sqrt[1] << norm) + r_shift(sqrt[0], FLINT_BITS - norm), sqrt[0] << norm, qlast << norm); 
        rem[0] >>= norm;
    }
    else
    {
        q = sqrt[0]/qlast;
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

mp_limb_t n_factor_SQUFOF(mp_limb_t n, ulong iters)
{
    mp_limb_t factor = _ll_factor_SQUFOF(UWORD(0), n, iters);
    mp_limb_t multiplier;
    mp_limb_t quot, rem;
    ulong i;
    
    for (i = 1; (i < FLINT_NUM_PRIMES_SMALL) && !factor; i++)
    {
        mp_limb_t multn[2];
        multiplier = flint_primes_small[i];
        umul_ppmm(multn[1], multn[0], multiplier, n);
        factor = _ll_factor_SQUFOF(multn[1], multn[0], iters);

        if (factor) 
        {
            quot = factor/multiplier;
            rem = factor - quot*multiplier;
            if (!rem) factor = quot;
            if ((factor == UWORD(1)) || (factor == n)) factor = UWORD(0);
        }
    }

    if (i == FLINT_NUM_PRIMES_SMALL) return UWORD(0);

    return factor; 
}
