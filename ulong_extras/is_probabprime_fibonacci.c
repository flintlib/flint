/*
    Copyright (C) 2008 Peter Shrimpton
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

n_pair_t
fchain_precomp(mp_limb_t m, mp_limb_t n, double npre)
{
    n_pair_t current = {0, 0}, old;
    int length;
    mp_limb_t power, xy;

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

n_pair_t
fchain2_preinv(mp_limb_t m, mp_limb_t n, mp_limb_t ninv)
{
    n_pair_t current = {0, 0}, old;
    int length;
    mp_limb_t power, xy;

    old.x = UWORD(2);
    old.y = n - UWORD(3);

    length = FLINT_BIT_COUNT(m);
    power = (UWORD(1) << (length - 1));

    for (; length > 0; length--)
    {
        xy = n_mulmod2_preinv(old.x, old.y, n, ninv);

        xy = n_addmod(xy, UWORD(3), n);

        if (m & power)
        {
            current.y =
                n_submod(n_mulmod2_preinv(old.y, old.y, n, ninv), UWORD(2), n);
            current.x = xy;
        }
        else
        {
            current.x =
                n_submod(n_mulmod2_preinv(old.x, old.x, n, ninv), UWORD(2), n);
            current.y = xy;
        }

        power >>= 1;
        old = current;
    }

    return current;
}

int
n_is_probabprime_fibonacci(mp_limb_t n)
{
    mp_limb_t m;
    n_pair_t V;

    if (FLINT_ABS((mp_limb_signed_t) n) <= UWORD(3))
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
        mp_limb_t ninv = n_preinvert_limb(n);

        V = fchain2_preinv(m, n, ninv);
        return (n_mulmod2_preinv(n - UWORD(3), V.x, n, ninv) ==
                n_mulmod2_preinv(UWORD(2), V.y, n, ninv));
    }
}
