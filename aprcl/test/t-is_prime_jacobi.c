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

    Copyright (C) 2015 Vladimir Glazachev

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "aprcl.h"

#include <time.h>

int main(void)
{
    int i;
    FLINT_TEST_INIT(state);
   
    flint_printf("is_prime_jacobi....");
    fflush(stdout);

    /* Test _is_prime_jacobi_check_pk() */
    {
        ulong p, q, v, k, p_pow;
        unity_zp j;
        fmpz_t n, u;

        q = 19;
        p = 3;
        k = 2;
        p_pow = n_pow(p, k);

        fmpz_init(u);
        fmpz_init_set_ui(n, 31);

        fmpz_tdiv_q_ui(u, n, p_pow);
        v = fmpz_tdiv_ui(n, p_pow);

        unity_zp_init(j, p, k, n);
        jacobi_pq(j, q, p);

        if (_is_prime_jacobi_check_pk(j, u, v) < 0)
        {
            flint_printf("FAIL\n");
            flint_printf("_is_prime_jacobi_check_pk() wrong answer");
            abort();
        }

        unity_zp_clear(j);
        fmpz_clear(n);
        fmpz_clear(u);
    }

    /* Test _is_prime_jacobi_check_22() */
    {
        ulong p, q, v, k, p_pow;
        unity_zp j;
        fmpz_t n, u;

        q = 13;
        p = 2;
        k = 2;
        p_pow = n_pow(p, k);

        fmpz_init(u);
        fmpz_init_set_ui(n, 31);

        fmpz_tdiv_q_ui(u, n, p_pow);
        v = fmpz_tdiv_ui(n, p_pow);

        unity_zp_init(j, p, k, n);
        jacobi_pq(j, q, p);

        if (_is_prime_jacobi_check_22(j, u, v, q) < 0)
        {
            flint_printf("FAIL\n");
            flint_printf("_is_prime_jacobi_check_pk() wrong answer");
            abort();
        }

        unity_zp_clear(j);
        fmpz_clear(n);
        fmpz_clear(u);

    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("NOT FINISHED\n");
    return 0;
}

