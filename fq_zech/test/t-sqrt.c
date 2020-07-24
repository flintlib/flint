/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>

#include "fq_zech.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int j, i, result;
    fq_zech_ctx_t ctx;

    FLINT_TEST_INIT(state);
    
    flint_printf("sqrt... ");
    fflush(stdout);

    for (j = 0; j < 10; j++)
    {
        fq_zech_ctx_randtest(ctx, state);

        /* Check aliasing: a = a * a */
        for (i = 0; i < 200; i++)
        {
            fq_zech_t a, c;

            fq_zech_init(a, ctx);
            fq_zech_init(c, ctx);

            fq_zech_randtest(a, state, ctx);

            fq_zech_sqr(c, a, ctx);

	    fq_zech_sqrt(a, c, ctx);
            fq_zech_sqrt(c, c, ctx);

            result = (fq_zech_equal(a, c, ctx));
            if (!result)
            {
                flint_printf("FAIL (aliasing):\n\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("c = "), fq_zech_print_pretty(c, ctx), flint_printf("\n");
                abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(c, ctx);
        }

        /* Check sqrt(a^2) = a */
        for (i = 0; i < 200; i++)
        {
            fq_zech_t a, b, c, d;

            fq_zech_init(a, ctx);
            fq_zech_init(b, ctx);
            fq_zech_init(c, ctx);
            fq_zech_init(d, ctx);

            fq_zech_randtest(a, state, ctx);

            fq_zech_sqr(b, a, ctx);

            fq_zech_sqrt(c, b, ctx);
            fq_zech_sqr(d, c, ctx);

            result = (fq_zech_equal(d, b, ctx));
            if (!result)
            {
                flint_printf("FAIL (sqrt(a^2) == a):\n\n");
                flint_printf("a = "), fq_zech_print_pretty(a, ctx), flint_printf("\n");
                flint_printf("b = "), fq_zech_print_pretty(b, ctx), flint_printf("\n");
                flint_printf("c = "), fq_zech_print_pretty(c, ctx), flint_printf("\n");
                flint_printf("d = "), fq_zech_print_pretty(c, ctx), flint_printf("\n");
                abort();
            }

            fq_zech_clear(a, ctx);
            fq_zech_clear(b, ctx);
            fq_zech_clear(c, ctx);
            fq_zech_clear(d, ctx);
        }

        fq_zech_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");

    return EXIT_SUCCESS;
}
