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
            flint_printf("in function _is_prime_jacobi_check_pk() wrong answer\n");
            abort();
        }

        unity_zp_clear(j);
        fmpz_clear(n);
        fmpz_clear(u);
    }

    /* Test _is_prime_jacobi_check_21() */
    {
        ulong p, q, v, k, p_pow;
        fmpz_t n;

        q = 7;
        p = 2;
        k = 1;
        p_pow = n_pow(p, k);

        fmpz_init_set_ui(n, 101);

        if (_is_prime_jacobi_check_21(q, n) < 0)
        {
            flint_printf("FAIL\n");
            flint_printf("in function _is_prime_jacobi_check_21() wrong answer\n");
            abort();
        }

        fmpz_clear(n);
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
        fmpz_init_set_ui(n, 101);

        fmpz_tdiv_q_ui(u, n, p_pow);
        v = fmpz_tdiv_ui(n, p_pow);

        unity_zp_init(j, p, k, n);
        jacobi_pq(j, q, p);

        if (_is_prime_jacobi_check_22(j, u, v, q) < 0)
        {
            flint_printf("FAIL\n");
            flint_printf("in function _is_prime_jacobi_check_22() wrong answer\n");
            abort();
        }

        unity_zp_clear(j);
        fmpz_clear(n);
        fmpz_clear(u);
    }


    /* Test _is_prime_jacobi_check_2k() */
    {
        ulong p, q, v, k, p_pow;
        unity_zp j, j2_1, j2_2;
        fmpz_t n, u;

        q = 41;
        p = 2;
        k = 3;
        p_pow = n_pow(p, k);

        fmpz_init(u);
        fmpz_init_set_ui(n, 101);

        fmpz_tdiv_q_ui(u, n, p_pow);
        v = fmpz_tdiv_ui(n, p_pow);

        unity_zp_init(j, p, k, n);
        unity_zp_init(j2_1, p, k, n);
        unity_zp_init(j2_2, p, k, n);

        jacobi_pq(j, q, p);
        jacobi_2q_one(j2_1, q);
        jacobi_2q_two(j2_2, q);

        flint_printf("\nresult = %w\n", _is_prime_jacobi_check_2k(j, j2_1, j2_2, u, v));

        if (_is_prime_jacobi_check_2k(j, j2_1, j2_2, u, v) < 0)
        {
            flint_printf("FAIL\n");
            flint_printf("in function _is_prime_jacobi_check_2k() wrong answer\n");
            abort();
        }

        unity_zp_clear(j);
        unity_zp_clear(j2_1);
        unity_zp_clear(j2_2);
        fmpz_clear(n);
        fmpz_clear(u);
    }

    /* Test is_prime_jacobi. */
    {
        fmpz_t n;
        aprcl_config config;

        fmpz_init(n);
        fmpz_set_str(n, "521419622856657689423872613771", 10);
        aprcl_config_init_min_R(config, n, 720);

        flint_printf("\ngocheck:\n");

        flint_printf("STATUS = %i\n", _is_prime_jacobi(n, config));

        aprcl_config_clear(config);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("NOT FINISHED\n");
    return 0;
}

