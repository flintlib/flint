/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "dlog.h"
#include <math.h>

#define vbs 0

/* TODO: tune the limit dlog -> index calculus */
void
dlog_vec_sieve_precomp(ulong *v, ulong nv, dlog_precomp_t pre,  ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    ulong smooth = 0, sievecount = 0, logcount = 0, missed = 0;
    ulong logcost;
#if 0
    ulong limcount;
#endif
    ulong k, p, pmax, logm1;
    n_primes_t iter;
    ulong X, aX, vaX;

    dlog_vec_fill(v, nv, DLOG_NOT_FOUND);
    v[1] = 0;

    logm1 = (na % 2) ? 0 : nmod_mul(na / 2, va, order);

    /* discrete log on first primes, then sieve */
    pmax = (nv < mod.n) ? nv : mod.n;
    logcost = pre->cost;

#if 0
    if (logcost < 15)
    {
        /* p1 = pmax; */
        limcount = mod.n;
    }
    else
    {
        limcount = ceil(pow((double)mod.n,1./2.3) * 40 / logcost);
    }
#endif

    /* take big power of gen */
    X = n_nextprime(3 * na / 2, 0) % na;
    aX = nmod_pow_ui(a, X, mod);
    vaX = nmod_mul(va, X % order.n, order);

    n_primes_init(iter);
    while ((p = n_primes_next(iter)) < pmax)
    {
        double cost;
        ulong m, vp;
        if (mod.n % p == 0)
            continue; /* won't be attained another time */
        cost = log(mod.n)/log(p);
        cost = pow(cost,cost);
        sievecount++;
        /* if (p < p1 || (wp = logp_sieve(w, nv, p, mod.n, logm1, order, logcost)) == NOT_FOUND) */
        /* if (smooth < limcount || (wp = logp_sieve_factor(w, nv, p, mod.n, a, na, va, logm1, order, logcost)) == NOT_FOUND)*/
        if (logcost < cost || (vp = dlog_vec_pindex_factorgcd(v, nv, p, mod, aX, na, vaX, logm1, order, cost)) == DLOG_NOT_FOUND)
        {
            if (logcost < cost)
                sievecount--;
            else
                missed++;
            logcount++;
            vp = nmod_mul(dlog_precomp(pre, p), va, order);
        }
        for (k = p, m = 1; k < nv; k += p, m++)
        {
            if (v[m] == DLOG_NOT_FOUND)
                continue;
            smooth++;
            v[k] = nmod_add(v[m],  vp, order);
        }
    }
#if vbs
    if (missed)
        flint_printf("[sieve: got %wu / %wu, n = %wu, cost %wu, logs %wu, sieve %wu missed %wu]\n",
                smooth, limcount, mod.n, logcost, logcount, sievecount, missed);
#endif
    n_primes_clear(iter);

    for (k = mod.n + 1; k < nv; k++)
        v[k] = v[k - mod.n];

}
