/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, poly_factor_iterated_frobenius_preinv, state)
{
    int i, j;

    for (j = 0; j < 20 * flint_test_multiplier(); j++)
    {
        int result;
        fmpz_t q;
        slong n;

        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, poly_t) v, vinv, *h1, *h2;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        fmpz_init(q);
        TEMPLATE(T, ctx_order) (q, ctx);

        TEMPLATE(T, poly_init) (v, ctx);
        TEMPLATE(T, poly_init) (vinv, ctx);

        TEMPLATE(T, poly_randtest_monic) (v, state, n_randint(state, 20) + 1,
                                          ctx);

        TEMPLATE(T, poly_reverse) (vinv, v, v->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (vinv, vinv, v->length, ctx);

        n = n_randint(state, 5) + 2;
        if (!(h1 = flint_malloc((2 * n) * sizeof(TEMPLATE(T, poly_struct)))))
        {
            flint_printf("Exception (t-fq_poly_iterated_frobenius_preinv):\n");
            flint_printf("Not enough memory.\n");
            fflush(stdout);
            flint_abort();
        }
        h2 = h1 + n;

        for (i = 0; i < 2 * n; i++)
            TEMPLATE(T, poly_init) (h1[i], ctx);

        TEMPLATE(T, poly_iterated_frobenius_preinv) (h1, n, v, vinv, ctx);

        TEMPLATE(T, poly_gen) (h2[0], ctx);
        for (i = 1; i < n; i++)
            TEMPLATE(T, poly_powmod_fmpz_sliding_preinv) (h2[i], h2[i - 1], q,
                                                          0, v, vinv, ctx);

        result = 1;
        for (i = 0; i < n; i++)
            result = result && TEMPLATE(T, poly_equal) (h1[i], h2[i], ctx);

        if (!result)
        {
            flint_printf("FAIL (composition):\n");
            flint_printf("v:\n");
            TEMPLATE(T, poly_print) (v, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, poly_clear) (v, ctx);
        TEMPLATE(T, poly_clear) (vinv, ctx);

        for (i = 0; i < 2 * n; i++)
            TEMPLATE(T, poly_clear) (h1[i], ctx);

        flint_free(h1);
        TEMPLATE(T, ctx_clear) (ctx);
        fmpz_clear(q);
    }

    TEST_FUNCTION_END(state);
}
#endif
