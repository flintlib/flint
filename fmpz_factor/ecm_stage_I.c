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
#include "fmpz.h"
#include "mpn_extras.h"

/* Implementation of the stage I of ECM */

int
fmpz_factor_ecm_stage_I(mp_ptr f, const mp_limb_t *prime_array, mp_limb_t num,
                        mp_limb_t B1, mp_ptr n, ecm_t ecm_inf)
{
    mp_limb_t times;
    mp_size_t sz, gcdlimbs;
    int i, j, p;

    for (i = 0; i < num; i++)
    {
        p = n_flog(B1, prime_array[i]);
        times = prime_array[i];

        for (j = 1; j <= p; j ++)
        {
            fmpz_factor_ecm_mul_montgomery_ladder(ecm_inf->x, ecm_inf->z,
                                                   ecm_inf->x, ecm_inf->z,
                                                   times, n, ecm_inf);

        }

        sz = ecm_inf->n_size;
        MPN_NORM(ecm_inf->z, sz);

        if (sz == 0)
            return 0;

        gcdlimbs = flint_mpn_gcd_full(f, n, ecm_inf->n_size, ecm_inf->z, sz);

        /* condition one -> gcd = n_ecm->one
           condition two -> gcd = n
           if neither is true, factor found */

        if ((((gcdlimbs == 1) && f[0] == ecm_inf->one[0]) || 
            ((gcdlimbs == ecm_inf->n_size) && mpn_cmp(f, n, ecm_inf->n_size) == 0)) == 0)
        {
            /* Found factor in stage I */
            return gcdlimbs;
        }
    }

    return 0;
}
