/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpq_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    FLINT_TEST_INIT(state);

    flint_printf("degrees/monomial_fits_ui/si....");
    fflush(stdout);

    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        const char * vars[] = {"x","y","z"};

        fmpq_mpoly_ctx_init(ctx, 3, ORD_DEGLEX);
        fmpq_mpoly_init(f, ctx);

        if (FLINT_BITS == 64)
        {
            fmpq_mpoly_set_str_pretty(f, "x^9223372036854775807", vars, ctx);
            if (!fmpq_mpoly_monomial_fits_si(f, 0, ctx)
                || !fmpq_mpoly_degrees_fit_si(f, 0, ctx))
            {
                printf("FAIL\nsi test 1\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "(x*y)^9223372036854775807", vars, ctx);
            if (!fmpq_mpoly_monomial_fits_si(f, 0, ctx)
                || !fmpq_mpoly_degrees_fit_si(f, 0, ctx))
            {
                printf("FAIL\nsi test 2\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^9223372036854775808", vars, ctx);
            if (fmpq_mpoly_monomial_fits_si(f, 0, ctx)
                || fmpq_mpoly_degrees_fit_si(f, 0, ctx))
            {
                printf("FAIL\nsi test 3\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "(x*y)^9223372036854775808", vars, ctx);
            if (fmpq_mpoly_monomial_fits_si(f, 0, ctx)
                || fmpq_mpoly_degrees_fit_si(f, 0, ctx))
            {
                printf("FAIL\nsi test 4\n");
                flint_abort();
            }


            fmpq_mpoly_set_str_pretty(f, "x^18446744073709551615", vars, ctx);
            if (!fmpq_mpoly_monomial_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 1\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "(x*y)^18446744073709551615", vars, ctx);
            if (!fmpq_mpoly_monomial_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 2\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^18446744073709551616", vars, ctx);
            if (fmpq_mpoly_monomial_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 3\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "(x*y)^18446744073709551616", vars, ctx);
            if (fmpq_mpoly_monomial_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 4\n");
                flint_abort();
            }

        } else if (FLINT_BITS == 32)
        {
            fmpq_mpoly_set_str_pretty(f, "x^2147483647", vars, ctx);
            if (!fmpq_mpoly_monomial_fits_si(f, 0, ctx)
                || !fmpq_mpoly_degrees_fit_si(f, 0, ctx))
            {
                printf("FAIL\nsi test 1\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^2147483647*y^2147483647", vars, ctx);
            if (!fmpq_mpoly_monomial_fits_si(f, 0, ctx)
                || !fmpq_mpoly_degrees_fit_si(f, 0, ctx))
            {
                printf("FAIL\nsi test 2\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^2147483648", vars, ctx);
            if (fmpq_mpoly_monomial_fits_si(f, 0, ctx)
                || fmpq_mpoly_degrees_fit_si(f, 0, ctx))
            {
                printf("FAIL\nsi test 3\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^2147483648*y^2147483648", vars, ctx);
            if (fmpq_mpoly_monomial_fits_si(f, 0, ctx)
                || fmpq_mpoly_degrees_fit_si(f, 0, ctx))
            {
                printf("FAIL\nsi test 4\n");
                flint_abort();
            }


            fmpq_mpoly_set_str_pretty(f, "x^4294967295", vars, ctx);
            if (!fmpq_mpoly_monomial_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 1\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^4294967295*y^4294967295", vars, ctx);
            if (!fmpq_mpoly_monomial_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 2\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^4294967296", vars, ctx);
            if (fmpq_mpoly_monomial_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 3\n");
                flint_abort();
            }

            fmpq_mpoly_set_str_pretty(f, "x^4294967296*y^4294967296", vars, ctx);
            if (fmpq_mpoly_monomial_fits_ui(f, 0, ctx))
            {
                printf("FAIL\nui test 4\n");
                flint_abort();
            }

        } else
        {
            printf("FAIL\nFLINT_BITS is not 64 or 32\n");
            flint_abort();
        }




        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

