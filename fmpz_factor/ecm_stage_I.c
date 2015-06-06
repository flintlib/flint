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

/* Implementation of the stage I of ECM */

int
fmpz_factor_ecm_stage_I(fmpz_t f, fmpz_t x, fmpz_t z, 
                       const mp_limb_t *prime_array, mp_limb_t num, 
                       mp_limb_t B1, fmpz_t a24, fmpz_t n)
{
    fmpz_t times;
    fmpz_init(times);

    int i, j, p, ret;

    ret = 0;

    for (i = 0; i < num; i++)
    {
        p = n_flog(B1, prime_array[i]);
        fmpz_set_ui(times, prime_array[i]);

        for (j = 1; j <= p; j ++)
            fmpz_factor_ecm_mul_montgomery_ladder(x, z, x, z, times, a24, n);

        fmpz_gcd(f, z, n);
        if (!fmpz_is_one(f) && fmpz_cmp(f, n))
        {
            /* Found factor after stage I */
            ret = 1;
            goto cleanup;
        }
    }

    cleanup:
    fmpz_clear(times);

    return ret;
}
