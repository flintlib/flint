/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "templates.h"
#include "fq_zech_poly.h"

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
        fq_zech_ctx_t ctx;
        fq_zech_poly_t a;
	nmod_poly_t b;
	ulong p;
	fq_zech_t r, s;

        len = n_randint(state, 15) + 1;
        fq_zech_ctx_randtest(ctx, state);
        fq_zech_poly_init(a, ctx);
        nmod_poly_init(b, ctx->p);
	fq_zech_init(r, ctx); fq_zech_init(s, ctx);

        nmod_poly_randtest(b, state, len);
	p = n_randint(state, 10);
	
        fq_zech_poly_set_mod_poly(a, b, ctx);
	fq_zech_set_ui(r, p, ctx);
	fq_zech_poly_evaluate_fq_zech(r, a, r, ctx);
	p = nmod_poly_evaluate_nmod(b, p);
	fq_zech_set_ui(s, p, ctx);

        result = fq_zech_equal(r, s, ctx);
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("CTX\n"), fq_zech_ctx_print(ctx),
                flint_printf("\n");
            flint_printf("b = "), nmod_poly_print_pretty(b, "X"),
                flint_printf("\n");
            flint_printf("p = %u\n", p);
            abort();
        }

        fq_zech_clear(r, ctx); fq_zech_clear(s, ctx);
	nmod_poly_clear(b);
	fq_zech_poly_clear(a, ctx);
        fq_zech_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
