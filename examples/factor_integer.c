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
#include "flint.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "fmpz_mpoly.h"
#include "profiler.h"

int
main(int argc, char * argv[])
{
    fmpz_t n;
    fmpz_factor_t fac;
    slong i;
    int num_threads = 1;
    int timing = 0;

    if (argc < 2)
    {
        flint_printf("usage: factor_integer [-threads t] [-timing] n\n");
        flint_printf("n can be given as an expression (no spaces)\n");
        return 1;
    }

    fmpz_init(n);

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
        else
        {
            /* allow input like 2^64+1 */
            /* we don't have an expression parser for fmpz yet, so use
               fmpz_mpoly with 0 variables  */
            {
                fmpz_mpoly_ctx_t mctx;
                fmpz_mpoly_t f;

                fmpz_mpoly_ctx_init(mctx, 0, ORD_LEX);
                fmpz_mpoly_init(f, mctx);

                if (fmpz_mpoly_set_str_pretty(f, argv[i], NULL, mctx) != 0)
                {
                    flint_printf("unable to parse integer\n");
                    return 1;
                }

                fmpz_mpoly_get_fmpz(n, f, mctx);

                fmpz_mpoly_clear(f, mctx);
                fmpz_mpoly_ctx_clear(mctx);
            }
        }
    }

    fmpz_factor_init(fac);

    if (timing)
    {
        TIMEIT_START
        fmpz_factor(fac, n);
        TIMEIT_STOP
        SHOW_MEMORY_USAGE
    }
    else
    {
        fmpz_factor(fac, n);
    }

    fmpz_print(n);
    flint_printf(" =\n");
    if (fac->sign != 1 || fac->num == 0)
    {
        flint_printf("%d", fac->sign);
        if (fac->num > 0)
            flint_printf(" * ");
    }
    for (i = 0; i < fac->num; i++)
    {
        fmpz_print(fac->p + i);
        if (fac->exp[i] >= 2)
            flint_printf("^%lu", fac->exp[i]);
        if (i < fac->num - 1)
            flint_printf(" * ");
    }
    flint_printf("\n");

    fmpz_factor_clear(fac);
    fmpz_clear(n);

    flint_cleanup_master();
    return 0;
}
