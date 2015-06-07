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

/* Outer wrapper for ECM 
   makes calls to stage I and stage II */

int
fmpz_factor_ecm(fmpz_t f, mp_limb_t curves, flint_rand_t state, fmpz_t n)
{
    mp_limb_t B1, B2, num;
    int j, ret;

    B1 = 1000;

    fmpz_t x, z, a, sig, nm8, a24;
    fmpz_init(x);
    fmpz_init(z);


    fmpz_init(a);
    fmpz_init(sig);
    fmpz_init(nm8);
    fmpz_init(a24);
    fmpz_sub_ui(nm8, n, 8);
    ret = 0;

    /* STAGE I PRECOMPUTATIONS */

    num = n_prime_pi(B1);   /* number of primes under B1 */

    /* compute list of primes under B1 for stage I */
    const mp_limb_t *prime_array = flint_malloc(num * sizeof(mp_limb_t));
    prime_array = n_primes_arr_readonly(num);   


    for (j = 0; j < curves; j++)
    {
        fmpz_set_ui(z, 1);
        fmpz_randm(sig, state, nm8);
        fmpz_add_ui(sig, sig, 7);

        if (fmpz_factor_ecm_select_curve(f, x, a, sig, n)) 
        {
            /* Found factor while selecting curve,
               very very lucky :) */
            ret = 1;
            goto cleanup;
        }        

        fmpz_add_ui(a24, a, 2);
        fmpz_fdiv_q_2exp(a24, a24, 2);

        if (fmpz_factor_ecm_stage_I(f, x, z, prime_array, num, B1, a24, n))
        {
            /* Found factor after stage I */
            ret = 1;
            goto cleanup;
        }

        if(fmpz_factor_ecm_stage_II(f, x, z, B1, 100*B1, a24, n))
        {
            /* Found factor after stage II */
            ret = 1;
            goto cleanup;
        }             
    }

    cleanup:

    fmpz_clear(a);
    fmpz_clear(sig);
    fmpz_clear(nm8);
    fmpz_clear(a24);

    return ret;
}
