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

    Copyright (C) 2015 Nitin Kumar

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <flint.h>
#include <ulong_extras.h>

/* computes primitive root for modulo n where n has following form
   2, 4, p ^ k, 2 * (p ^ k), p is an odd prime*/

mp_limb_t
n_primitive_root_prefactor(mp_limb_t n, n_factor_t * factors, mp_limb_t * phi)
{
    slong i;
    int found;
    mp_limb_t result, a, pm1;
    double pinv;

    pm1 = *phi;
    pinv = n_precompute_inverse(n);

    for (a = 2; a < n; a++)
    {
        if (n_gcd(n, a) == 1)
        {
            found = 1;
            for (i = 0; i < factors->num; i++)
            {
                result = n_powmod_precomp(a, pm1 / factors->p[i], n, pinv);
                if (result == 1)
                {
                    found = 0;
                    break;
                }
            }
            if (found)
            {
                return a;
            }
        }
    }
    flint_printf("Exception (n_primitive_root_prime_prefactor).  root not found.\n");
    abort();
}

mp_limb_t n_primitive_root(mp_limb_t n)
{
    if (n_is_prime(n)) return n_primitive_root_prime(n);
    else
    {
        if (n == 4) return 3;
        else if (n % 4 == 0)
        {
            flint_printf("Exception (n_primitive_root_prime_prefactor).  root not found.\n");
            return 0;
        }
        else
        {
            mp_limb_t a, m, root, phi;
            int k;
            n_factor_t factors;
            n_factor_init(&factors);
            if (n % 2 == 0)
            {
                m = 2;
                n /= 2;
            }
            else
            {
                m = 1;
            }
            k = n_is_perfect_power(&root, n);
            if (n_is_prime(root) != 1)
            {
                flint_printf("Exception (n_primitive_root_prime_prefactor).  root not found.\n");
                return 0;
            }
            else
            {
                phi = n;
                /* calculating  phi(n), using formula */
                phi = (phi / root) * (root - 1);
                /* factoring (root - 1), which is a factor of phi(n) */
                n_factor(&factors, root - 1, 1);
                n *= m;
                /* adding extra factor if exist */
                if (k > 1)
                {
                    factors.p[factors.num] = root;
                    factors.exp[factors.num] = 1;
                    factors.num++;
                }
                a = n_primitive_root_prefactor(n, &factors, &phi);
                return a;
            }
        }
    }
}
