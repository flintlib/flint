/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_theta.h"

static int
sp2gz_comes_from_g1(const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t x, y;
    slong k, l;
    int res = 0;

    fmpz_mat_init(x, 2, 2);
    fmpz_mat_init(y, 2 * g, 2 * g);

    for (k = 0; k < 2; k++)
    {
        for (l = 0; l < 2; l++)
        {
            fmpz_set(fmpz_mat_entry(x, k, l), fmpz_mat_entry(mat, k * g, l * g));
        }
    }

    if (sp2gz_is_correct(x))
    {
        sp2gz_embed(y, x);
        res = fmpz_mat_equal(mat, y);
    }

    fmpz_mat_clear(x);
    fmpz_mat_clear(y);
    return res;
}

static int
sp2gz_is_allowed_in_dec(const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t alpha, beta, gamma, x, y;
    slong r;
    int res;

    if (g == 1 || sp2gz_comes_from_g1(mat))
    {
        return 1;
    }

    fmpz_mat_window_init(alpha, mat, 0, 0, g, g);
    fmpz_mat_window_init(beta, mat, 0, g, g, 2 * g);
    fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
    fmpz_mat_init(x, 2 * g, 2 * g);

    if (!fmpz_mat_is_zero(gamma))
    {
        r = fmpz_mat_rank(gamma);

        fmpz_mat_init(y, 2 * r, 2 * r);
        sp2gz_j(y);
        sp2gz_embed(x, y);
        fmpz_mat_clear(y);
    }
    else if (!fmpz_mat_is_zero(beta))
    {
        sp2gz_trig(x, beta);
    }
    else
    {
        sp2gz_block_diag(x, alpha);
    }

    res = fmpz_mat_equal(mat, x);
    fmpz_mat_window_clear(alpha);
    fmpz_mat_window_clear(beta);
    fmpz_mat_window_clear(gamma);
    fmpz_mat_clear(x);
    return res;
}

TEST_FUNCTION_START(acb_theta_sp2gz_decompose, state)
{
    slong iter;

    /* Test: decomposition consists of elementary matrices and product is the
       original matrix */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 2 + n_randint(state, 5);
        slong bits = n_randint(state, 20);
        fmpz_mat_t m, x;
        fmpz_mat_struct * dec = NULL;
        slong nb_dec = 0;
        slong k;

        fmpz_mat_init(m, 2 * g, 2 * g);
        fmpz_mat_init(x, 2 * g, 2 * g);

        sp2gz_randtest(m, state, bits);
        dec = sp2gz_decompose(&nb_dec, m);

        for (k = 0; k < nb_dec; k++)
        {
            if (!sp2gz_is_allowed_in_dec(&dec[k]))
            {
                flint_printf("FAIL (not elementary)\n");
                fmpz_mat_print_pretty(&dec[k]);
                flint_printf("\n");
                flint_abort();
            }
        }

        fmpz_mat_one(x);
        for (k = 0; k < nb_dec; k++)
        {
            fmpz_mat_mul(x, x, &dec[k]);
        }
        if (!fmpz_mat_equal(m, x))
        {
            flint_printf("FAIL (product)\n");
            fmpz_mat_print_pretty(x);
            flint_printf("\n");
            fmpz_mat_print_pretty(m);
            flint_printf("\ndecomposition in %wd matrices:\n", nb_dec);
            for (k = 0; k < nb_dec; k++)
            {
                fmpz_mat_print_pretty(&dec[k]);
                flint_printf("\n");
            }
            flint_abort();
        }

        fmpz_mat_clear(m);
        fmpz_mat_clear(x);
        for (k = 0; k < nb_dec; k++)
        {
            fmpz_mat_clear(&dec[k]);
        }
        flint_free(dec);
    }

    TEST_FUNCTION_END(state);
}

