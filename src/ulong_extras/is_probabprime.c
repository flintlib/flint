/*
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Dana Jacobsen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"

typedef struct
{
    ulong x;
    ulong y;
}
n_pair_t;

/* Obsolete function; n_is_prime is faster than n_is_probabprime_BPSW,
   so just make n_is_probabprime a wrapper for n_is_prime. */
int n_is_probabprime(ulong n)
{
    return n_is_prime(n);
}

int
n_is_probabprime_BPSW(ulong n)
{
    if (n <= UWORD(1))
        return 0;

    if ((n & UWORD(1)) == UWORD(0))
    {
        if (n == UWORD(2))
            return 1;
        return 0;
    }

    if (((n % 10) == 3) || ((n % 10) == 7))
    {
        if (n_is_probabprime_fermat(n, 2) == 0)
            return 0;

        return n_is_probabprime_fibonacci(n);
    }
    else
    {
        ulong d;

        d = n - UWORD(1);
        while ((d & UWORD(1)) == UWORD(0))
            d >>= 1;

        if (FLINT_BIT_COUNT(n) <= FLINT_D_BITS)
        {
            double npre = n_precompute_inverse(n);
            if (n_is_strong_probabprime_precomp(n, npre, WORD(2), d) == 0)
                return 0;
        }
        else
        {
            ulong ninv = n_preinvert_limb(n);
            if (n_is_strong_probabprime2_preinv(n, ninv, WORD(2), d) == 0)
                return 0;
        }

        return (n_is_probabprime_lucas(n) == 1);
    }
}

int
n_is_probabprime_fermat(ulong n, ulong i)
{
    if (FLINT_BIT_COUNT(n) <= FLINT_D_BITS)
        return (n_powmod(i, n - 1, n) == UWORD(1));
    else
        return n_powmod2_ui_preinv(i, n - 1, n, n_preinvert_limb(n)) == UWORD(1);
}

static n_pair_t
fchain_precomp(ulong m, ulong n, double npre)
{
    n_pair_t current = {0, 0}, old;
    int length;
    ulong power, xy;

    old.x = UWORD(2);
    old.y = n - UWORD(3);

    length = FLINT_BIT_COUNT(m);
    power = (UWORD(1) << (length - 1));

    for (; length > 0; length--)
    {
        xy = n_mulmod_precomp(old.x, old.y, n, npre);

        xy = n_addmod(xy, UWORD(3), n);

        if (m & power)
        {
            current.y =
                n_submod(n_mulmod_precomp(old.y, old.y, n, npre), UWORD(2), n);
            current.x = xy;
        }
        else
        {
            current.x =
                n_submod(n_mulmod_precomp(old.x, old.x, n, npre), UWORD(2), n);
            current.y = xy;
        }

        power >>= 1;
        old = current;
    }

    return current;
}

/* the operands are reduced: nmod_mul rather than n_mulmod2_preinv */
static n_pair_t
fchain2_nmod(ulong m, nmod_t mod)
{
    const ulong n = mod.n;
    n_pair_t current = {0, 0}, old;
    int length;
    ulong power, xy;

    old.x = UWORD(2);
    old.y = n - UWORD(3);

    length = FLINT_BIT_COUNT(m);
    power = (UWORD(1) << (length - 1));

    for (; length > 0; length--)
    {
        xy = nmod_mul(old.x, old.y, mod);

        xy = n_addmod(xy, UWORD(3), n);

        if (m & power)
        {
            current.y =
                n_submod(nmod_mul(old.y, old.y, mod), UWORD(2), n);
            current.x = xy;
        }
        else
        {
            current.x =
                n_submod(nmod_mul(old.x, old.x, mod), UWORD(2), n);
            current.y = xy;
        }

        power >>= 1;
        old = current;
    }

    return current;
}

int
n_is_probabprime_fibonacci(ulong n)
{
    ulong m;
    n_pair_t V;

    if ((ulong) FLINT_ABS((slong) n) <= UWORD(3))
    {
        if (n >= UWORD(2))
            return 1;
        return 0;
    }

    m = (n - n_jacobi(WORD(5), n)) / 2;  /* cannot overflow
                                       as (5/n) = 0 for n = 2^64-1 */

    if (FLINT_BIT_COUNT(n) <= FLINT_D_BITS)
    {
        double npre = n_precompute_inverse(n);

        V = fchain_precomp(m, n, npre);
        return (n_mulmod_precomp(n - UWORD(3), V.x, n, npre) ==
                n_mulmod_precomp(UWORD(2), V.y, n, npre));
    }
    else
    {
        nmod_t mod;
        nmod_init(&mod, n);

        V = fchain2_nmod(m, mod);
        return (nmod_mul(n - UWORD(3), V.x, mod) ==
                nmod_mul(UWORD(2), V.y, mod));
    }
}

static n_pair_t
lchain_precomp(ulong m, ulong a, ulong n, double npre)
{
    n_pair_t current = {0, 0}, old;
    int length, i;
    ulong power, xy, xx, yy;

    old.x = UWORD(2);
    old.y = a;

    length = FLINT_BIT_COUNT(m);
    power = (UWORD(1) << (length - 1));

    for (i = 0; i < length; i++)
    {
        xy = n_submod(n_mulmod_precomp(old.x, old.y, n, npre), a, n);

        if (m & power)
        {
            yy = n_submod(n_mulmod_precomp(old.y, old.y, n, npre), UWORD(2), n);
            current.x = xy;
            current.y = yy;
        }
        else
        {
            xx = n_submod(n_mulmod_precomp(old.x, old.x, n, npre), UWORD(2), n);
            current.x = xx;
            current.y = xy;
        }

        power >>= 1;
        old = current;
    }

    return current;
}

/* the operands are reduced: nmod_mul rather than n_mulmod2_preinv */
static n_pair_t
lchain2_nmod(ulong m, ulong a, nmod_t mod)
{
    const ulong n = mod.n;
    n_pair_t current = {0, 0}, old;
    int length, i;
    ulong power, xy, xx, yy;

    old.x = UWORD(2);
    old.y = a;

    length = FLINT_BIT_COUNT(m);
    power = (UWORD(1) << (length - 1));

    for (i = 0; i < length; i++)
    {
        xy = n_submod(nmod_mul(old.x, old.y, mod), a, n);

        if (m & power)
        {
            yy = n_submod(nmod_mul(old.y, old.y, mod), UWORD(2), n);
            current.x = xy;
            current.y = yy;
        }
        else
        {
            xx = n_submod(nmod_mul(old.x, old.x, mod), UWORD(2), n);
            current.x = xx;
            current.y = xy;
        }

        power >>= 1;
        old = current;
    }

    return current;
}

int
n_is_probabprime_lucas(ulong n)
{
    int i;
    slong D, Q;
    ulong A;
    ulong left, right;
    n_pair_t V;

    D = 0;
    Q = 0;

    if (((n % 2) == 0) || (FLINT_ABS((slong) n) <= 2))
    {
        return (n == UWORD(2));
    }

    for (i = 0; i < 100; i++)
    {
        D = 5 + 2 * i;
        if (n_gcd(D, n % D) != UWORD(1))
        {
            if (n == (ulong) D)
                continue;
            else
                return 0;
        }
        if (i % 2 == 1)
            D = -D;
        if (n_jacobi(D, n) == -1)
            break;
    }

    if (i == 100)
    {
        return (n_is_square(n) ? -1 : 1);
    }

    Q = (1 - D) / 4;
    if (Q < 0)
    {
        if (n < UWORD(52))
        {
            while (Q < 0)
                Q += n;
            A = n_submod(n_invmod(Q, n), UWORD(2), n);
        }
        else
            A = n_submod(n_invmod(Q + n, n), UWORD(2), n);
    }
    else
    {
        if (n < UWORD(52))
        {
            while ((ulong) Q >= n)
                Q -= n;
            A = n_submod(n_invmod(Q, n), UWORD(2), n);
        }
        else
            A = n_submod(n_invmod(Q, n), UWORD(2), n);
    }

    if (FLINT_BIT_COUNT(n) <= FLINT_D_BITS)
    {
        double npre = n_precompute_inverse(n);
        V = lchain_precomp(n + 1, A, n, npre);

        left = n_mulmod_precomp(A, V.x, n, npre);
        right = n_mulmod_precomp(2, V.y, n, npre);
    }
    else
    {
        nmod_t mod;
        nmod_init(&mod, n);
        V = lchain2_nmod(n + 1, A, mod);

        left = nmod_mul(A, V.x, mod);
        right = nmod_mul(2, V.y, mod);
    }

    return (left == right);
}
