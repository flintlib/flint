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

    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int
n_factor_ecm_stage_I(mp_limb_t *f, const mp_limb_t *prime_array, mp_limb_t num,
                     mp_limb_t B1, mp_limb_t n, n_ecm_t n_ecm_inf)
{
    mp_limb_t times;
    int i, j, p;

    for (i = 0; i < num; i++)
    {
        p = n_flog(B1, prime_array[i]);
        times = prime_array[i];

        for (j = 1; j <= p; j ++)
        {
            n_factor_ecm_mul_montgomery_ladder(&(n_ecm_inf->x), &(n_ecm_inf->z),
                                               n_ecm_inf->x, n_ecm_inf->z, 
                                               times, n, n_ecm_inf);
        }
        
        *f = n_gcd(n_ecm_inf->z, n);

        if ((*f > n_ecm_inf->one) && (*f < n))
        {
            /* Found factor in stage I */
            return 1;
        }
    }

    return 0;
}
