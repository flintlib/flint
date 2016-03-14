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

int main(void)
{
    int i, j;
    FLINT_TEST_INIT(state);
   
    flint_printf("unity_zp_jacobi_sum....");
    fflush(stdout);

    /*
        test computation of
        j(\chi_{p, q}, \chi_{p, q})
    */
    for (i = 0; i < 100; i++)
    {
        ulong ind, q, p, k;
        fmpz_t n;
        unity_zp f, g;
        n_factor_t q_factors;
        mp_ptr table;

        n_factor_init(&q_factors);

        q = n_randprime(state, 2 + n_randint(state, 4), 0);
        while (q < 3)
            q = n_randprime(state, 2 + n_randint(state, 4), 0);

        n_factor(&q_factors, q - 1, 1);
        ind = n_randint(state, q_factors.num);
        p = q_factors.p[ind];
        k = q_factors.exp[ind];

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        table = f_table(q);

        unity_zp_init(f, p, k, n);
        unity_zp_init(g, p, k, n);

        for (j = 1; j < q - 1; j++)
        {
            unity_zp_coeff_inc(g, j + table[j]);
        }

        _unity_zp_reduce_cyclotomic(g);
        unity_zp_jacobi_sum_pq(f, q, p);

        if (unity_zp_equal(f, g) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_clear(n);
        unity_zp_clear(f);
        unity_zp_clear(g);
        _nmod_vec_clear(table);
    }

    /*
        test computation of
        j(\chi_{2, q}^{2^{k - 3}}, \chi_{2, q}^{3 * 2^{k - 3}}))^2
    */
    for (i = 0; i < 100; i++)
    {
        ulong q, p, k;
        fmpz_t n;
        unity_zp f, g;
        n_factor_t q_factors;
        mp_ptr table;

        n_factor_init(&q_factors);

        q = n_randprime(state, 2 + n_randint(state, 4), 0);
        while (q < 3)
            q = n_randprime(state, 2 + n_randint(state, 4), 0);

        n_factor(&q_factors, q - 1, 1);
        p = 2;
        k = q_factors.exp[0];

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        table = f_table(q);

        unity_zp_init(f, p, k, n);
        unity_zp_init(g, p, k, n);

        for (j = 1; j < q - 1; j++)
        {
            unity_zp_coeff_inc(g, 2 * j + table[j]);
        }

        _unity_zp_reduce_cyclotomic(g);
        unity_zp_jacobi_sum_2q_one(f, q);

        if (unity_zp_equal(f, g) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_clear(n);
        unity_zp_clear(f);
        unity_zp_clear(g);
        _nmod_vec_clear(table);
    }

    /*
        test computation of
        j(\chi_{2, q}, \chi_{2, q}, \chi_{2, q})
    */
    for (i = 0; i < 100; i++)
    {
        ulong a, b, q, p, k;
        fmpz_t n;
        unity_zp f, g;
        n_factor_t q_factors;
        mp_ptr table;

        n_factor_init(&q_factors);

        q = n_randprime(state, 2 + n_randint(state, 6), 0);
        while (q < 3)
            q = n_randprime(state, 2 + n_randint(state, 6), 0);

        n_factor(&q_factors, q - 1, 1);
        p = 2;
        k = q_factors.exp[0];

        if (k < 3) continue;

        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);
        while (fmpz_equal_ui(n, 0) != 0)
            fmpz_randtest_unsigned(n, state, 200);

        table = f_table(q);

        unity_zp_init(f, p, k, n);
        unity_zp_init(g, p, k, n);

        b = n_pow(2, k - 3);
        a = 3 * b;

        for (j = 1; j < q - 1; j++)
        {
            unity_zp_coeff_inc(g, a * j + b * table[j]);
        }

        _unity_zp_reduce_cyclotomic(g);
        unity_zp_jacobi_sum_2q_two(f, q);

        if (unity_zp_equal(f, g) == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }

        fmpz_clear(n);
        unity_zp_clear(f);
        unity_zp_clear(g);
        _nmod_vec_clear(table);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

