/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("spd_radius....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: random matrix within radius is still positive definite */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 10);
        slong prec = 200 + n_randint(state, 500);
        slong mag_bits = 1 + n_randint(state, 5);
        arb_mat_t A, B, cho;
        arf_t rad, c;
        slong k, j;
        int res;

        arb_mat_init(A, g, g);
        arb_mat_init(B, g, g);
        arb_mat_init(cho, g, g);
        arf_init(rad);
        arf_init(c);

        arb_mat_randtest_spd(A, state, prec, mag_bits);
        arb_mat_spd_radius(rad, A, prec);

        if (!arf_is_finite(rad) || arf_cmp_si(rad, 0) <= 0)
        {
            flint_printf("FAIL (not positive)\n");
            flint_printf("g = %wd, prec = %wd, mag_bits = %wd\n", g, prec, mag_bits);
            arb_mat_printd(A, 5);
            arf_printd(rad, 10);
            flint_printf("\n");
            flint_abort();
        }

        for (k = 0; k < g; k++)
        {
            for (j = 0; j <= k; j++)
            {
                /* Get c between -rad and rad */
                arf_urandom(c, state, prec, ARF_RND_FLOOR);
                arf_mul_2exp_si(c, c, 1);
                arf_sub_si(c, c, 1, prec, ARF_RND_DOWN);
                arf_mul(c, c, rad, prec, ARF_RND_DOWN);

                arb_add_arf(arb_mat_entry(B, k, j), arb_mat_entry(A, k, j),
                    c, prec);
                arb_set(arb_mat_entry(B, j, k), arb_mat_entry(B, k, j));
            }
        }

        res = arb_mat_cho(cho, B, prec);
        if (!res)
        {
            flint_printf("FAIL (cholesky)\n");
            flint_printf("g = %wd, prec = %wd, mag_bits = %wd\n", g, prec, mag_bits);
            arb_mat_printd(A, 5);
            flint_printf("radius: ");
            arf_printd(rad, 10);
            flint_printf("\n");
            flint_printf("Deformed matrix:\n");
            arb_mat_printd(B, 5);
            flint_abort();
        }

        arb_mat_clear(A);
        arb_mat_clear(B);
        arb_mat_clear(cho);
        arf_clear(rad);
        arf_clear(c);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
