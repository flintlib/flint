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

#ifdef __GNUC__
# define floor __builtin_floor
#else
# include <math.h>
#endif

static ulong
dlog_single(ulong b, ulong a, const nmod_t mod, ulong n)
{
    if (n < 50)
    {
        int k;
        ulong ak = 1;

        for (k=0; k < n; k++)
        {
            if (ak == b)
                return k;
            ak = nmod_mul(ak, a, mod);
        }

        flint_throw(FLINT_ERROR, "(dlog single): log(%wu,%wu) mod %wu not found (size %wu)\n", b, a, mod.n, n);
    }
    else
    {
        dlog_rho_t t;
        dlog_rho_init(t, a, mod.n, n);
        return dlog_rho(t, b);
    }
}

/* solve log knowing equation  e = f * log(b) [n] */
static ulong
dlog_quotient(const dlog_rho_t t, ulong e, ulong f, ulong g, ulong b)
{
    ulong r, b_ar, an;
    nmod_t n = t->n;

    if (g == n.n)
    {
        flint_throw(FLINT_ERROR, "(dlog quotient): trivial relation e = %wu, f = %wu mod %wu\n", e, f, n.n);
    }

    nmod_init(&n, n.n / g);
    e = e / g;
    f = f / g;
    r = nmod_div(e, f, n);
    an = nmod_pow_ui(t->a, n.n, t->mod);
    b_ar = nmod_div(b, nmod_pow_ui(t->a, r, t->mod), t->mod);

    return r + n.n * dlog_single(b_ar, an, t->mod, g);
}

#define RWALK 20
ulong
dlog_rho(const dlog_rho_t t, ulong b)
{
    int j, k, l;
    ulong m[RWALK], n[RWALK], ab[RWALK];
    ulong x[2], e[2], f[2], g;
    flint_rand_t state;

    flint_randinit(state);

    do {

        for (k = 0; k < RWALK; k++)
        {
            m[k] = 1 + n_randint(state, t->n.n - 1);
            n[k] = 1 + n_randint(state, t->n.n - 1);
            ab[k] = nmod_mul(nmod_pow_ui(t->a, m[k], t->mod), nmod_pow_ui(b, n[k], t->mod), t->mod);
        }

        /* x[l] = a^e[l] * b^f[l] */
        x[0] = x[1] = 1;
        e[0] = e[1] = 0;
        f[0] = f[1] = 0;

        do {

            for(j = 0; j < 3; j++)
            {
                l = (j > 0);
                k = floor( (double) RWALK * x[l] / t->mod.n );
                x[l] = nmod_mul(x[l], ab[k], t->mod);
                e[l] = nmod_add(e[l], m[k], t->n);
                f[l] = nmod_add(f[l], n[k], t->n);
            }

        } while (x[0] != x[1]);

    } while (e[0] == e[1] && f[0] == f[1]);

    flint_randclear(state);

    /* e = f * log(b) */
    e[0] = nmod_sub(e[0], e[1], t->n);
    f[0] = nmod_sub(f[1], f[0], t->n);

    if (!t->nisprime && (g = n_gcd(f[0], t->n.n)) > 1)
        return dlog_quotient(t, e[0], f[0], g, b);
    else
        return nmod_div(e[0], f[0], t->n);
}
