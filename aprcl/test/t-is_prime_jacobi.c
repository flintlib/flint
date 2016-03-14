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
        unity_zp_jacobi_sum_pq(j, q, p);

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
        ulong q;
        fmpz_t n;

        q = 7;

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
        fmpz_init_set_ui(n, 659);

        fmpz_tdiv_q_ui(u, n, p_pow);
        v = fmpz_tdiv_ui(n, p_pow);

        unity_zp_init(j, p, k, n);
        unity_zp_jacobi_sum_pq(j, q, p);

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

        unity_zp_jacobi_sum_pq(j, q, p);
        unity_zp_jacobi_sum_2q_one(j2_1, q);
        unity_zp_jacobi_sum_2q_two(j2_2, q);

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
        for (i = 0; i < 10000; i++)
        {
            int pbprime, cycloprime;
            fmpz_t n;
            fmpz_init(n);

            fmpz_randtest_unsigned(n, state, 1000);
            while (fmpz_cmp_ui(n, 100) <= 0)
                fmpz_randtest_unsigned(n, state, 1000);

            pbprime = fmpz_is_probabprime(n);
            cycloprime = is_prime_jacobi(n);
        
            if (pbprime != cycloprime)
            {
                flint_printf("FAIL\n");
                flint_printf("Testing number = ");
                fmpz_print(n);
                flint_printf("\nis_probabprime = %i, is_prime_gauss = %i\n", pbprime, cycloprime);
                abort();
            }

            fmpz_clear(n);
        }

        {
            int result;
            fmpz_t n;
            fmpz_init(n);
            result = 1;

            /* Test big primes. */
            fmpz_set_str(n, "40206835204840513073", 10);
            if (is_prime_jacobi(n) == 0)
                result = 0;

            /* 521419622856657689423872613771 % 4 == 3 */
            fmpz_set_str(n, "521419622856657689423872613771", 10);
            if (is_prime_jacobi(n) == 0)
                result = 0;

            /* 5991810554633396517767024967580894321153 % 4 == 1 */
            fmpz_set_str(n, "5991810554633396517767024967580894321153", 10);
            if (is_prime_jacobi(n) == 0)
                result = 0;

            fmpz_set_str(n, "2074722246773485207821695222107608587480996474721117292752992589912196684750549658310084416732550077", 10);
            if (is_prime_jacobi(n) == 0)
                result = 0;

            /* prime power */
            fmpz_pow_ui(n, n, 7);
            if (is_prime_jacobi(n) == 1)
                result = 0;

            /* Test big composite. */
            /* 1500450271 * 5915587277 */
            fmpz_set_str(n, "8876044532898802067", 10);
            if (is_prime_jacobi(n) == 1)
                result = 0;

            /* 5915587277 * 54673257461630679457 */
            fmpz_set_str(n, "323424426232167763068694468589", 10);
            if (is_prime_jacobi(n) == 1)
                result = 0;

            /* 48112959837082048697 * 66405897020462343733 */
            fmpz_set_str(n, "3194984256290911228520362769161858765901", 10);
            if (is_prime_jacobi(n) == 1)
                result = 0;

            if (result == 0)
            {
                flint_printf("FAIL\n");
                abort();
            }

            fmpz_clear(n);
        }

    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

