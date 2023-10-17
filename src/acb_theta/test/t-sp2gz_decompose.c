/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("sp2gz_decompose....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: decomposition consist of block-diagonal, trigonal or J matrices,
       and product is the original matrix */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        slong bits = n_randint(state, 100);
        fmpz_mat_t m, x, alpha, beta, gamma;
        fmpz_mat_struct* dec = NULL;
        slong nb_dec = 0;
        slong k;

        fmpz_mat_init(m, 2 * g, 2 * g);
        fmpz_mat_init(x, 2 * g, 2 * g);

        sp2gz_randtest(m, state, bits);
        dec = sp2gz_decompose(&nb_dec, m);

        for (k = 0; k < nb_dec; k++)
        {
            fmpz_mat_window_init(alpha, &dec[k], 0, 0, g, g);
            fmpz_mat_window_init(beta, &dec[k], 0, g, g, 2 * g);
            fmpz_mat_window_init(gamma, &dec[k], g, 0, 2 * g, g);

            if (!fmpz_mat_is_zero(gamma))
            {
                sp2gz_j(x);
            }
            else if (!fmpz_mat_is_zero(beta))
            {
                sp2gz_trig(x, beta);
            }
            else
            {
                sp2gz_block_diag(x, alpha);
            }

            if (!fmpz_mat_equal(&dec[k], x))
            {
                flint_printf("FAIL (not elementary)\n");
                fmpz_mat_print_pretty(&dec[k]);
                flint_printf("\n");
                flint_abort();
            }

            fmpz_mat_window_clear(alpha);
            fmpz_mat_window_clear(beta);
            fmpz_mat_window_clear(gamma);
        }

        fmpz_mat_one(x);
        for (k = 0; k < nb_dec; k++)
        {
            fmpz_mat_mul(x, &dec[k], x);
        }
        if (!fmpz_mat_equal(m, x))
        {
            flint_printf("FAIL (product)\n");
            fmpz_mat_print_pretty(x);
            flint_printf("\n");
            fmpz_mat_print_pretty(m);
            flint_printf("\n");
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

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

