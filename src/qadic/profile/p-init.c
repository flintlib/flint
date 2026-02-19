/*
    Copyright (C) 2026 Alexey Orlov

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Benchmarks for the context creation, when no conway polynomial is available.

    Currently we test only "small" primes, to ensure that q-adic context creation
      does not "pessimize" in this case (nmod_poly routines should be used)
 */

#include <time.h>
#include <unistd.h>

#include "fmpz.h"
#include "qadic.h"

int
main(void)
{
    const ulong p_small[] = {
        2, 3, 5, 7, 11,
    };
    const ulong ps_degree[] = {
        123, 628, 1024, 2048, 2755, 2867,
    };

    fmpz_t p;
    qadic_ctx_t ctx;

    clock_t c0, c1;

    c0 = clock();
    for (int i = 0, in = sizeof(p_small)/sizeof(p_small[0]) ; i < in ; ++i)
    {
        fmpz_init_set_ui(p, p_small[i]);
        for (int j = 0, jn = sizeof(ps_degree)/sizeof(ps_degree[0]) ; j < jn ; ++j)
        {
            qadic_ctx_init(ctx, p, ps_degree[j], 0, 1, "x", PADIC_SERIES);
            qadic_ctx_clear(ctx);
        }
        fmpz_clear(p);
    }
    c1 = clock();

    flint_printf("%.4Lf s\n", (long double)(c1 - c0) / (long double)CLOCKS_PER_SEC);

    return 0;
}
