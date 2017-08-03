/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "templates.h"
#include "fq_nmod_poly.h"

#include <stdio.h>
#include <stdlib.h>

#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("set_mod_poly... ");
    fflush(stdout);

    /* Check litfed polynomials by evaluating at random points */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        slong len;
        fq_nmod_ctx_t ctx;
        fq_nmod_poly_t a;
	nmod_poly_t b;
	ulong p;
	fq_nmod_t r, s;

        len = n_randint(state, 15) + 1;
        fq_nmod_ctx_randtest(ctx, state);
        fq_nmod_poly_init(a, ctx);
        nmod_poly_init(b, ctx->p);
	fq_nmod_init(r, ctx); fq_nmod_init(s, ctx);

        nmod_poly_randtest(b, state, len);
	p = n_randint(state, 10);
	
        fq_nmod_poly_set_mod_poly(a, b, ctx);
	fq_nmod_set_ui(r, p, ctx);
	fq_nmod_poly_evaluate_fq_nmod(r, a, r, ctx);
	p = nmod_poly_evaluate_nmod(b, p);
	fq_nmod_set_ui(s, p, ctx);

        result = fq_nmod_equal(r, s, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), fq_nmod_ctx_print(ctx),
                flint_printf("\n");
            flint_printf("b = "), nmod_poly_print_pretty(b, "X"),
                flint_printf("\n");
            flint_printf("p = %u\n", p);
            abort();
        }

        fq_nmod_clear(r, ctx); fq_nmod_clear(s, ctx);
	nmod_poly_clear(b);
	fq_nmod_poly_clear(a, ctx);
        fq_nmod_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
