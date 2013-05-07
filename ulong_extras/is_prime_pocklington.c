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
#undef ulong /* prevent clash with standard library */
#include <stdlib.h>
#define ulong unsigned long
#include "flint.h"
#include "ulong_extras.h"

int
n_is_prime_pocklington(mp_limb_t n, ulong iterations)
{
    int i, j, pass;
    mp_limb_t n1, cofactor, b, c = 0, ninv, limit;
    n_factor_t factors;

    if (n % 2 == 0)
    {
        return (n == 2UL);
    }

    n1 = n - 1;

    n_factor_init(&factors);

    limit = n_sqrt(n1);
    cofactor = n_factor_partial(&factors, n1, limit, 1);

    if (cofactor != 1) /* check that cofactor is coprime to factors found */
    {
        for (i = 0; i < factors.num; i++)
        {
            if (factors.p[i] > flint_primes[FLINT_FACTOR_TRIAL_PRIMES - 1])
            {
                while (cofactor >= factors.p[i] && (cofactor % factors.p[i]) == 0)
                {
                    factors.exp[i]++;
                    cofactor /= factors.p[i];
                }
            }
        }
    }

    ninv = n_preinvert_limb(n);

    c = 1;
    for (i = factors.num - 1; i >= 0; i--)
    {
        mp_limb_t exp = n1 / factors.p[i];
        pass = 0;

        for (j = 2; j < iterations && pass == 0; j++)
        {
            b = n_powmod2_preinv(j, exp, n, ninv);
            if (n_powmod2_ui_preinv(b, factors.p[i], n, ninv) != 1UL)
                return 0;

            b = n_submod(b, 1UL, n);
            if (b != 0UL)
            {
                c = n_mulmod2_preinv(c, b, n, ninv);
                pass = 1;
            }

            if (c == 0)
                return 0;
        }

        if (j == iterations)
            return -1;
    }

    return (n_gcd(n, c) == 1UL);
}
