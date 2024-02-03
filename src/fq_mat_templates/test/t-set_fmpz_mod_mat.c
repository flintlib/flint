/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod.h"

TEST_TEMPLATE_FUNCTION_START(T, mat_set_fmpz_mod_mat, state)
{
    int i, result;

    /* Check conversion of identity matrix */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) a;
        fmpz_mod_mat_t m;
        slong r, c;
        fmpz_mod_ctx_t pctx;

        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);

        r = n_randint(state, 10);
        c = n_randint(state, 10);
        TEMPLATE(T, mat_init)(a, r, c, ctx);

        TEMPLATE(T, mat_randtest)(a, state, ctx);

#if defined(FQ_NMOD_MAT_H) || defined(FQ_ZECH_MAT_H)
        fmpz_mod_ctx_init_ui(pctx, TEMPLATE(T, ctx_prime)(ctx));
#else
        fmpz_mod_ctx_init(pctx, TEMPLATE(T, ctx_prime)(ctx));
#endif

        fmpz_mod_mat_init(m, r, c, pctx);
        fmpz_mod_mat_one(m, pctx);

        TEMPLATE(T, mat_set_fmpz_mod_mat)(a, m, ctx);

        result = (TEMPLATE(T, mat_is_one)(a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "), TEMPLATE(T, mat_print)(a, ctx), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mat_clear(m, pctx);
        fmpz_mod_ctx_clear(pctx);

        TEMPLATE(T, mat_clear)(a, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
