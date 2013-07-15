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
lchain_precomp(mp_limb_t m, mp_limb_t a, mp_limb_t n, double npre)
{
    n_pair_t current = {0, 0}, old;
    int length, i;
    mp_limb_t power, xy, xx, yy;

    old.x = 2UL;
    old.y = a;

    length = FLINT_BIT_COUNT(m);
    power = (1UL << (length - 1));

    for (i = 0; i < length; i++)
    {
        xy = n_submod(n_mulmod_precomp(old.x, old.y, n, npre), a, n);

        if (m & power)
        {
            yy = n_submod(n_mulmod_precomp(old.y, old.y, n, npre), 2UL, n);
            current.x = xy;
            current.y = yy;
        }
        else
        {
            xx = n_submod(n_mulmod_precomp(old.x, old.x, n, npre), 2UL, n);
            current.x = xx;
            current.y = xy;
        }

        power >>= 1;
        old = current;
    }

    return current;
}

n_pair_t
lchain2_preinv(mp_limb_t m, mp_limb_t a, mp_limb_t n, mp_limb_t ninv)
{
    n_pair_t current = {0, 0}, old;
    int length, i;
    mp_limb_t power, xy, xx, yy;

    old.x = 2UL;
    old.y = a;

    length = FLINT_BIT_COUNT(m);
    power = (1UL << (length - 1));

    for (i = 0; i < length; i++)
    {
        xy = n_submod(n_mulmod2_preinv(old.x, old.y, n, ninv), a, n);

        if (m & power)
        {
            yy = n_submod(n_mulmod2_preinv(old.y, old.y, n, ninv), 2UL, n);
            current.x = xy;
            current.y = yy;
        }
        else
        {
            xx = n_submod(n_mulmod2_preinv(old.x, old.x, n, ninv), 2UL, n);
            current.x = xx;
            current.y = xy;
        }

        power >>= 1;
        old = current;
    }

    return current;
}

int
n_is_probabprime_lucas(mp_limb_t n)
{
    int i, D, Q;
    mp_limb_t A;
    mp_limb_t left, right;
    n_pair_t V;

    D = 0;
    Q = 0;

    if (((n % 2) == 0) || (FLINT_ABS((mp_limb_signed_t) n) <= 2))
    {
        return (n == 2UL);
    }

    for (i = 0; i < 100; i++)
    {
        D = 5 + 2 * i;
        if (n_gcd(D, n % D) != 1UL)
        {
            if (n == D)
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
        if (n < 52UL)
        {
            while (Q < 0)
                Q += n;
            A = n_submod(n_invmod(Q, n), 2UL, n);
        }
        else
            A = n_submod(n_invmod(Q + n, n), 2UL, n);
    }
    else
    {
        if (n < 52UL)
        {
            while (Q >= n)
                Q -= n;
            A = n_submod(n_invmod(Q, n), 2UL, n);
        }
        else
            A = n_submod(n_invmod(Q, n), 2UL, n);
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
        mp_limb_t ninv = n_preinvert_limb(n);
        V = lchain2_preinv(n + 1, A, n, ninv);

        left = n_mulmod_precomp(A, V.x, n, ninv);
        right = n_mulmod_precomp(2, V.y, n, ninv);
    }

    return (left == right);
}
