/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "dlog.h"

ulong
dlog_crt_init(dlog_crt_t t, ulong a, ulong mod, ulong n, ulong num)
{
    int k;
    n_factor_t fac;
    ulong * M, * u;
    ulong cost = 0;

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    t->num = fac.num;
    nmod_init(&t->mod,mod);
    nmod_init(&t->n, n);

    M = t->expo = flint_malloc(t->num * sizeof(ulong));
    u = t->crt_coeffs = flint_malloc(t->num * sizeof(ulong));
    t->pre = flint_malloc(t->num * sizeof(dlog_precomp_struct));

    for (k = 0; k < t->num; k++)
    {
        ulong p, e, mk;
        p = fac.p[k];
        e = fac.exp[k];
        if (0 && mod % p == 0)
        {
            flint_throw(FLINT_ERROR, "dlog_crt_init: modulus must be prime to order.\n");
        }
        mk = n_pow(p, e);
        M[k] = n / mk;
        u[k] = nmod_mul(M[k], n_invmod(M[k] % mk, mk), t->n);
        /* depends on the power */
#if 0
        flint_printf("[sub-crt -- init for size %wu mod %wu]\n", mk, mod);
#endif
        dlog_precomp_pe_init(t->pre + k, nmod_pow_ui(a, M[k], t->mod), mod, p, e, mk, num);
        cost += t->pre[k].cost;
    }
#if 0
    if (cost > 500)
    flint_printf("[crt init for size %wu mod %wu -> cost %wu]\n", n,mod,cost);
#endif
    return cost;
}
