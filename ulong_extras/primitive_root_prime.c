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

    Copyright (C) 2013 Mike Hansen
 
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_primitive_root_prime_prefactor(mp_limb_t p, n_factor_t * factors)
{
    slong i;
    int found;
    mp_limb_t result, a, pm1;
    double pinv;

    if (p == 2)
    {
        return 1;
    }

    pm1 = p - 1;
    pinv = n_precompute_inverse(p);

    for (a = 2; a < p; a++)
    {
        found = 1;
        for (i = 0; i < factors->num; i++)
        {
            result = n_powmod_precomp(a, pm1 / factors->p[i], p, pinv);
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
    flint_printf("Exception (n_primitive_root_prime_prefactor).  root not found.\n");
    flint_abort();
}

mp_limb_t n_primitive_root_prime(mp_limb_t p)
{
    mp_limb_t a;
    n_factor_t factors;

    n_factor_init(&factors);
    n_factor(&factors, p - 1, 1);
    
    a = n_primitive_root_prime_prefactor(p, &factors);

    return a;
}
