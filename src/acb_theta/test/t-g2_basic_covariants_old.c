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

    flint_printf("g2_basic_covariants_old....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: basic_covariants and covariants_old agree */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong prec = 100 + n_randint(state, 100);
        slong bits = 2;
        slong nb = ACB_THETA_G2_BASIC_NB;
        acb_poly_struct* r;
        acb_poly_struct* test;
        acb_poly_t f;
        slong k;

        r = flint_malloc(nb * sizeof(acb_poly_struct));
        test = flint_malloc(nb * sizeof(acb_poly_struct));
        for (k = 0; k < nb; k++)
        {
            acb_poly_init(&r[k]);
            acb_poly_init(&test[k]);
        }
        acb_poly_init(f);

        acb_poly_randtest(f, state, 7, prec, bits);

        acb_theta_g2_basic_covariants(r, f, prec);
        acb_theta_g2_basic_covariants_old(test, f, prec);

        for (k = 0; k < nb; k++)
        {
            if (!acb_poly_overlaps(&r[k], &test[k]))
            {
                flint_printf("FAIL (k = %wd)\n", k);
                acb_poly_printd(&r[k], 5);
                flint_printf("\n");
                acb_poly_printd(&test[k], 5);
                flint_printf("\n");
                flint_abort();
            }
        }

        for (k = 0; k < nb; k++)
        {
            acb_poly_clear(&r[k]);
            acb_poly_clear(&test[k]);
        }
        flint_free(r);
        flint_free(test);
        acb_poly_clear(f);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
