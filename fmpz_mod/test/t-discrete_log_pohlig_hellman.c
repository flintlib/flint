/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

int
main(void)
{
    slong i, j, k;
    FLINT_TEST_INIT(state);

    flint_printf("discrete_log_pohlig_hellman....");
    fflush(stdout);
   
    flint_randinit(state);

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_t p, pm1;
        fmpz_mod_ctx_t fpctx;
        fmpz_mod_discrete_log_pohlig_hellman_t L;

        fmpz_init_set_ui(p, 2);
        fmpz_init_set_ui(pm1, 1);
        fmpz_mod_ctx_init(fpctx, p);
        fmpz_mod_discrete_log_pohlig_hellman_init(L);

        for (j = 0; j < 10; j++)
        {
            double score;

            fmpz_set_ui(p, n_randtest_prime(state, 1));
            fmpz_sub_ui(pm1, p, 1);
            fmpz_mod_ctx_set_mod(fpctx, p);
            score = fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(L, p);
            if (score > 10000)
            {
                continue;
            }
            
            for (k = 0; k < 10; k++)
            {
                fmpz_t x, y, xr;
                const fmpz * alpha = fmpz_mod_discrete_log_pohlig_hellman_primitive_root(L);
                fmpz_init(x);
                fmpz_init(y);
                fmpz_init(xr);

                fmpz_randm(x, state, pm1);
                fmpz_mod_pow_fmpz(y, alpha, x, fpctx);
                fmpz_mod_discrete_log_pohlig_hellman_run(L, xr, y);
                if (!fmpz_equal(x, xr))
                {
                    printf("FAIL\n");
                    flint_printf("i = %wd, j = %wd, k = %wd\n", i, j, k);
                    flint_abort();
                }

                fmpz_clear(x);
                fmpz_clear(y);
                fmpz_clear(xr);
            }
        }
        fmpz_mod_discrete_log_pohlig_hellman_clear(L);
        fmpz_mod_ctx_clear(fpctx);
        fmpz_clear(p);
        fmpz_clear(pm1);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
