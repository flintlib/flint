/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_factor.h>
#include <flint/profiler.h>

int
main(int argc, char * argv[])
{
    fmpz_mpoly_t f;
    fmpz_mpoly_ctx_t mctx;
    fmpz_mpoly_factor_t fac;
    slong i;
    int num_threads = 1;
    int timing = 0;
    int suppress = 0;
    const char * vars[] = { "x", "y", "z" };

    if (argc < 2)
    {
        flint_printf("usage: factor_polynomial [-threads t] [-timing] [-suppress] expr\n");
        flint_printf("example: build/examples/factor_polynomial -timing \"(1+x+2*y+3*z)^3+1\"\n");
        return 1;
    }

    fmpz_mpoly_ctx_init(mctx, 3, ORD_LEX);
    fmpz_mpoly_init(f, mctx);
    fmpz_mpoly_factor_init(fac, mctx);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-threads"))
        {
            num_threads = atoi(argv[i+1]);
            flint_set_num_threads(num_threads);
            i++;
        }
        else if (!strcmp(argv[i], "-timing"))
        {
            timing = 1;
        }
        else if (!strcmp(argv[i], "-suppress"))
        {
            suppress = 1;
        }
        else
        {
            if (fmpz_mpoly_set_str_pretty(f, argv[i], vars, mctx) != 0)
            {
                flint_printf("unable to parse polynomial\n");
                return 1;
            }
        }
    }

    if (timing)
    {
        TIMEIT_START
        fmpz_mpoly_factor(fac, f, mctx);
        TIMEIT_STOP
        SHOW_MEMORY_USAGE
    }
    else
    {
        fmpz_mpoly_factor(fac, f, mctx);
    }

    if (!suppress)
    {
        fmpz_mpoly_print_pretty(f, vars, mctx);
        flint_printf(" =\n");
        fmpz_mpoly_factor_print_pretty(fac, vars, mctx);
        flint_printf("\n");
    }
    else
    {
        flint_printf("%wd irreducible factors\n", fac->num);
    }

    fmpz_mpoly_factor_clear(fac, mctx);
    fmpz_mpoly_clear(f, mctx);
    fmpz_mpoly_ctx_clear(mctx);

    flint_cleanup_master();
    return 0;
}
