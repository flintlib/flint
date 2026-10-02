/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Atkin's information at a prime l where E has no rational l-isogeny.

    Frobenius then has conjugate eigenvalues lambda, lambda^l in F_{l^2}
    outside F_l, and their ratio gamma = lambda^(l-1) has norm one and
    exact order r, the degree of the irreducible factors of the modular
    polynomial at j (see _atkin_order in elkies.c). Since t = lambda +
    lambda^l and q = lambda^(l+1),

        t^2 / q = gamma + 1/gamma + 2,

    so t mod l is one of the square roots of q (gamma + 1/gamma + 2) for
    gamma of order r in the norm one subgroup of F_{l^2}^*, which is cyclic
    of order l + 1. gamma and 1/gamma give the same value, so there are at
    most phi(r) candidates, often far fewer than l.

    F_{l^2} is F_l[w]/(w^2 - n) for the least non-residue n; an element of
    norm one has inverse its conjugate, so gamma + 1/gamma is twice its
    rational part.
*/

#include "ulong_extras.h"
#include "nmod.h"
#include "gr_ec.h"
#include "impl.h"

typedef struct
{
    ulong a, b;         /* a + b w */
}
_fl2_struct;

static _fl2_struct
_fl2_mul(_fl2_struct x, _fl2_struct y, ulong n, nmod_t mod)
{
    _fl2_struct z;
    ulong t;

    t = nmod_mul(nmod_mul(x.b, y.b, mod), n, mod);
    z.a = nmod_add(nmod_mul(x.a, y.a, mod), t, mod);
    z.b = nmod_add(nmod_mul(x.a, y.b, mod), nmod_mul(x.b, y.a, mod), mod);

    return z;
}

static _fl2_struct
_fl2_pow(_fl2_struct x, ulong e, ulong n, nmod_t mod)
{
    _fl2_struct z;

    z.a = 1;
    z.b = 0;

    while (e != 0)
    {
        if (e & 1)
            z = _fl2_mul(z, x, n, mod);

        x = _fl2_mul(x, x, n, mod);
        e >>= 1;
    }

    return z;
}

static int
_fl2_is_one(_fl2_struct x)
{
    return x.a == 1 && x.b == 0;
}

/*
    The candidates for t mod l given r, into T (room for l + 1 entries),
    sorted and without repetition. Returns their number.
*/
slong
_gr_ec_atkin_candidates(ulong * T, ulong l, ulong r, ulong ql)
{
    nmod_t mod;
    n_factor_t fac;
    _fl2_struct x, gen, h, g;
    ulong n, k, e, v, s, i;
    slong len = 0, a, b;

    if (l < 3 || r < 2 || (l + 1) % r != 0 || ql % l == 0)
        return 0;

    nmod_init(&mod, l);

    for (n = 2; n_jacobi(n, l) != -1; n++)
        ;

    /* a generator of the norm one subgroup: x^(l-1) for a suitable x */
    n_factor_init(&fac);
    n_factor(&fac, l + 1, 1);

    for (k = 0; ; k++)
    {
        int ok = 1;

        x.a = k;
        x.b = 1;
        gen = _fl2_pow(x, l - 1, n, mod);

        for (i = 0; i < (ulong) fac.num && ok; i++)
            if (_fl2_is_one(_fl2_pow(gen, (l + 1) / fac.p[i], n, mod)))
                ok = 0;

        if (ok)
            break;
    }

    /* the elements of order r are h^e, gcd(e, r) = 1; e and r - e pair up */
    h = _fl2_pow(gen, (l + 1) / r, n, mod);

    for (e = 1; 2 * e <= r; e++)
    {
        if (n_gcd(e, r) != 1)
            continue;

        g = _fl2_pow(h, e, n, mod);

        /* v = q (gamma + 1/gamma + 2) = q (2 a + 2) */
        v = nmod_add(g.a, 1, mod);
        v = nmod_add(v, v, mod);
        v = nmod_mul(v, ql % l, mod);

        if (v == 0)
            T[len++] = 0;
        else if (n_jacobi(v, l) == 1)
        {
            s = n_sqrtmod(v, l);
            T[len++] = s;
            T[len++] = l - s;
        }
    }

    /* sort and remove repetitions; the lists are short */
    for (a = 1; a < len; a++)
        for (b = a; b > 0 && T[b - 1] > T[b]; b--)
        {
            ulong t = T[b];
            T[b] = T[b - 1];
            T[b - 1] = t;
        }

    for (a = 0, b = 0; a < len; a++)
        if (b == 0 || T[b - 1] != T[a])
            T[b++] = T[a];

    return b;
}
