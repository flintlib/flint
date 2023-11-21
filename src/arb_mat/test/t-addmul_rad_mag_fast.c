/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"

#define BLOCK_SIZE 32

static void
fallback(arb_mat_t C, mag_srcptr A, mag_srcptr B, slong ar, slong ac, slong bc)
{
    slong i, j, k, ii, jj, kk;

    for (ii = 0; ii < ar; ii += BLOCK_SIZE)
    {
        for (jj = 0; jj < bc; jj += BLOCK_SIZE)
        {
            for (kk = 0; kk < ac; kk += BLOCK_SIZE)
            {
                for (i = ii; i < FLINT_MIN(ii + BLOCK_SIZE, ar); i++)
                {
                    for (j = jj; j < FLINT_MIN(jj + BLOCK_SIZE, bc); j++)
                    {
                        for (k = kk; k < FLINT_MIN(kk + BLOCK_SIZE, ac); k++)
                        {
                            mag_fast_addmul(arb_radref(arb_mat_entry(C, i, j)),
                                A + i * ac + k, B + j * ac + k);
                        }
                    }
                }
            }
        }
    }
}

TEST_FUNCTION_START(arb_mat_addmul_rad_mag_fast, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_mat_t C, D;
        slong m, n, p, density, i, j, off;
        mag_ptr A, B;
        mag_t lo, hi;

        m = n_randint(state, 40);
        n = n_randint(state, 40);
        p = n_randint(state, 40);

        A = _mag_vec_init(m * n);
        B = _mag_vec_init(p * n);
        arb_mat_init(C, m, p);
        arb_mat_init(D, m, p);
        mag_init(lo);
        mag_init(hi);

        density = 1 + n_randint(state, 100);

        off = n_randint(state, 10000);
        for (i = 0; i < m * n; i++)
        {
            if (n_randint(state, 100) < density)
            {
                mag_randtest(A + i, state, 8 + n_randint(state, 10));
                mag_mul_2exp_si(A + i, A + i, off);
            }
        }

        off = n_randint(state, 10000);
        for (i = 0; i < p * n; i++)
        {
            if (n_randint(state, 100) < density)
            {
                mag_randtest(B + i, state, 8 + n_randint(state, 10));
                mag_mul_2exp_si(B + i, B + i, off);
            }
        }

        fallback(C, A, B, m, n, p);
        _arb_mat_addmul_rad_mag_fast(D, A, B, m, n, p);

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < p; j++)
            {
                mag_set_d(lo, 0.9999);
                mag_set_d(hi, 1.0001);
                mag_mul(lo, lo, arb_radref(arb_mat_entry(C, i, j)));
                mag_mul(hi, hi, arb_radref(arb_mat_entry(C, i, j)));

                if (mag_cmp(arb_radref(arb_mat_entry(D, i, j)), lo) < 0 ||
                    mag_cmp(arb_radref(arb_mat_entry(D, i, j)), hi) > 0)
                {
                    flint_printf("FAIL\n");
                    flint_printf("m = %wd, n = %wd, p = %wd\n", m, n, p);
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    mag_printd(arb_radref(arb_mat_entry(C, i, j)), 10);
                    flint_printf("\n");
                    mag_printd(arb_radref(arb_mat_entry(D, i, j)), 10);
                    flint_printf("\n");
                    flint_abort();
                }
            }
        }

        _mag_vec_clear(A, m * n);
        _mag_vec_clear(B, p * n);
        arb_mat_clear(C);
        arb_mat_clear(D);
        mag_clear(lo);
        mag_clear(hi);
    }

    TEST_FUNCTION_END(state);
}
