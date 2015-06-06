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
fmpz_factor_ecm_stage1(fmpz_t f, flint_rand_t state, fmpz_t x, fmpz_t z, int curves, 
           const mp_limb_t *prime_array, mp_limb_t num, fmpz_t B1, fmpz_t n)
{
    int i, j, p, ok, ret;
    fmpz_t a, sig, nm6, a24, times;
    fmpz_init(a);
    fmpz_init(sig);
    fmpz_init(nm6);
    fmpz_init(a24);
    fmpz_init(times);

    fmpz_sub_ui(nm6, n, 8);
    ret = 0;

    for (j = 0; j < curves; j++)
    {
        fmpz_set_ui(z, 1);
        fmpz_randm(sig, state, nm6);
        fmpz_add_ui(sig, sig, 7);

        if (fmpz_factor_ecm_select_curve(f, x, a, sig, n)) 
        {
            /* Found factor while selecting curve,
               very very lucky :) */
            return 1;
        }

        fmpz_add_ui(a24, a, 2);
        fmpz_fdiv_q_2exp(a24, a24, 2);

        for (i = 0; i < num; i++)
        {
            p = n_flog(fmpz_get_ui(B1), prime_array[i]);
            fmpz_set_ui(times, prime_array[i]);

            for (ok = 1; ok <= p; ok ++)
                fmpz_factor_ecm_mul_montgomery_ladder(x, z, x, z, times, a24, n);

            fmpz_gcd(f, z, n);
            if (!fmpz_is_one(f) && fmpz_cmp(f, n))
            {
                /* Found factor after stage I */
                ret = 1;
                goto cleanup;
            }
        }

        if(ecm_stage2(f, x, z, fmpz_get_ui(B1), 100*fmpz_get_ui(B1), a24, n))
        {
            /* Found factor after stage II */
            ret = 1;
            goto cleanup;
        }
    }
    cleanup:

    fmpz_clear(a);
    fmpz_clear(sig);
    fmpz_clear(nm6);
    fmpz_clear(a24);
    fmpz_clear(times);

    return ret;
}