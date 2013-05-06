/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart
   
******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

n_pair_t
fchain_precomp(mp_limb_t m, mp_limb_t n, double npre)
{
    n_pair_t current = {0, 0}, old;
    int length;
    mp_limb_t power, xy;

    old.x = 2UL;
    old.y = n - 3UL;

    length = FLINT_BIT_COUNT(m);
    power = (1UL << (length - 1));

    for (; length > 0; length--)
    {
        xy = n_mulmod_precomp(old.x, old.y, n, npre);

        xy = n_addmod(xy, 3UL, n);

        if (m & power)
        {
            current.y =
                n_submod(n_mulmod_precomp(old.y, old.y, n, npre), 2UL, n);
            current.x = xy;
        }
        else
        {
            current.x =
                n_submod(n_mulmod_precomp(old.x, old.x, n, npre), 2UL, n);
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

    old.x = 2UL;
    old.y = n - 3UL;

    length = FLINT_BIT_COUNT(m);
    power = (1UL << (length - 1));

    for (; length > 0; length--)
    {
        xy = n_mulmod2_preinv(old.x, old.y, n, ninv);

        xy = n_addmod(xy, 3UL, n);

        if (m & power)
        {
            current.y =
                n_submod(n_mulmod2_preinv(old.y, old.y, n, ninv), 2UL, n);
            current.x = xy;
        }
        else
        {
            current.x =
                n_submod(n_mulmod2_preinv(old.x, old.x, n, ninv), 2UL, n);
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

    if (FLINT_ABS((mp_limb_signed_t) n) <= 3UL)
    {
        if (n >= 2UL)
            return 1;
        return 0;
    }

    m = (n - n_jacobi(5L, n)) / 2;  /* cannot overflow 
                                       as (5/n) = 0 for n = 2^64-1 */

    if (FLINT_BIT_COUNT(n) <= FLINT_D_BITS)
    {
        double npre = n_precompute_inverse(n);

        V = fchain_precomp(m, n, npre);
        return (n_mulmod_precomp(n - 3UL, V.x, n, npre) ==
                n_mulmod_precomp(2UL, V.y, n, npre));
    }
    else
    {
        mp_limb_t ninv = n_preinvert_limb(n);

        V = fchain2_preinv(m, n, ninv);
        return (n_mulmod2_preinv(n - 3UL, V.x, n, ninv) ==
                n_mulmod2_preinv(2UL, V.y, n, ninv));
    }
}
