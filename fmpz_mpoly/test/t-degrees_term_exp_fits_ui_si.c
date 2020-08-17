/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i, j, k;
    int result;
    FLINT_TEST_INIT(state);

    flint_printf("degrees_term_exp_fits_ui_si....");
    fflush(stdout);

    /* basic corner cases */
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        const char * vars[] = {"x","y","z"};

        fmpz_mpoly_ctx_init(ctx, 3, ORD_DEGLEX);
        fmpz_mpoly_init(f, ctx);

        if (FLINT_BITS == 64)
        {
            fmpz_mpoly_set_str_pretty(f, "x^9223372036854775807", vars, ctx);
            if (!fmpz_mpoly_term_exp_fits_si(f, 0, ctx)
                || !fmpz_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 1\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "(x*y)^9223372036854775807", vars, ctx);
            if (!fmpz_mpoly_term_exp_fits_si(f, 0, ctx)
                || !fmpz_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 2\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "x^9223372036854775808", vars, ctx);
            if (fmpz_mpoly_term_exp_fits_si(f, 0, ctx)
                || fmpz_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 3\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "(x*y)^9223372036854775808", vars, ctx);
            if (fmpz_mpoly_term_exp_fits_si(f, 0, ctx)
                || fmpz_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 4\n");
                flint_abort();
            }


            fmpz_mpoly_set_str_pretty(f, "x^18446744073709551615", vars, ctx);
            if (!fmpz_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 1\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "(x*y)^18446744073709551615", vars, ctx);
            if (!fmpz_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 2\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "x^18446744073709551616", vars, ctx);
            if (fmpz_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 3\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "(x*y)^18446744073709551616", vars, ctx);
            if (fmpz_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 4\n");
                flint_abort();
            }

        } else if (FLINT_BITS == 32)
        {
            fmpz_mpoly_set_str_pretty(f, "x^2147483647", vars, ctx);
            if (!fmpz_mpoly_term_exp_fits_si(f, 0, ctx)
                || !fmpz_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 1\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "x^2147483647*y^2147483647", vars, ctx);
            if (!fmpz_mpoly_term_exp_fits_si(f, 0, ctx)
                || !fmpz_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 2\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "x^2147483648", vars, ctx);
            if (fmpz_mpoly_term_exp_fits_si(f, 0, ctx)
                || fmpz_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 3\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "x^2147483648*y^2147483648", vars, ctx);
            if (fmpz_mpoly_term_exp_fits_si(f, 0, ctx)
                || fmpz_mpoly_degrees_fit_si(f, ctx))
            {
                printf("FAIL\nsi test 4\n");
                flint_abort();
            }


            fmpz_mpoly_set_str_pretty(f, "x^4294967295", vars, ctx);
            if (!fmpz_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 1\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "x^4294967295*y^4294967295", vars, ctx);
            if (!fmpz_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 2\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "x^4294967296", vars, ctx);
            if (fmpz_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 3\n");
                flint_abort();
            }

            fmpz_mpoly_set_str_pretty(f, "x^4294967296*y^4294967296", vars, ctx);
            if (fmpz_mpoly_term_exp_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 4\n");
                flint_abort();
            }

        } else
        {
            printf("FAIL\nFLINT_BITS is not 64 or 32\n");
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check fit */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        fmpz ** exp, ** deg;
        slong nvars, len;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);

        nvars = fmpz_mpoly_ctx_nvars(ctx);

        exp = (fmpz **) flint_malloc(nvars*sizeof(fmpz *));
        deg = (fmpz **) flint_malloc(nvars*sizeof(fmpz *));
        for (k = 0; k < nvars; k++)
        {
            exp[k] = (fmpz *) flint_malloc(sizeof(fmpz)); 
            deg[k] = (fmpz *) flint_malloc(sizeof(fmpz)); 
            fmpz_init(exp[k]);
            fmpz_init(deg[k]);
            fmpz_set_si(deg[k], -WORD(1));
        }

        len = n_randint(state, 100);
        exp_bits = n_randint(state, FLINT_BITS + 10) + 1;
        coeff_bits = n_randint(state, 30);

        fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < fmpz_mpoly_length(f, ctx); j++)
        {
            fmpz_mpoly_get_term_exp_fmpz(exp, f, j, ctx);
            for (k = 0; k < nvars; k++)
                if (fmpz_cmp(deg[k], exp[k]) < 0)
                    fmpz_set(deg[k], exp[k]);

            result = 1;
            for (k = 0; k < nvars; k++)
                result = result && fmpz_fits_si(exp[k]);
            if (result != fmpz_mpoly_term_exp_fits_si(f, j, ctx))
            {
                flint_printf("FAIL\nCheck monomial_fit_si\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }

            result = 1;
            for (k = 0; k < nvars; k++)
                result = result && fmpz_abs_fits_ui(exp[k]);
            if (result != fmpz_mpoly_term_exp_fits_ui(f, j, ctx))
            {
                flint_printf("FAIL\nCheck monomial_fit_ui\ni = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        result = 1;
        fmpz_mpoly_degrees_fmpz(exp, f, ctx);
        for (k = 0; k < nvars; k++)
        {
            if (!fmpz_equal(exp[k], deg[k]))
            {
                flint_printf("FAIL\nCheck degree computation\ni = %wd\n", i);
                flint_abort();
            }            
            result = result && fmpz_fits_si(exp[k]);
        }
        if (result != fmpz_mpoly_degrees_fit_si(f, ctx))
        {
            flint_printf("FAIL\nCheck degrees_fit_si\ni = %wd\n", i);
            flint_abort();
        }            

        for (k = 0; k < nvars; k++)
        {
            fmpz_clear(deg[k]);
            flint_free(deg[k]); 
            fmpz_clear(exp[k]);
            flint_free(exp[k]); 
        }
        flint_free(deg);
        flint_free(exp);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

