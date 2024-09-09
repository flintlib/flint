/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_default.h"

TEST_FUNCTION_START(fq_default_init, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_t fq;
        fmpz_t p;

        fmpz_init(p);
        fmpz_set_ui(p, 3);

        fq_default_ctx_init(ctx, p, 3, "x");
        fq_default_init(fq, ctx);
        fq_default_clear(fq, ctx);
        fq_default_randtest(fq, state, ctx);
        fq_default_ctx_clear(ctx);

        fq_default_ctx_init(ctx, p, 16, "x");
        fq_default_init(fq, ctx);
        fq_default_randtest(fq, state, ctx);
        fq_default_clear(fq, ctx);
        fq_default_ctx_clear(ctx);

        fmpz_set_str(p, "73786976294838206473", 10);
        fq_default_ctx_init(ctx, p, 1, "x");
        fq_default_init(fq, ctx);
        fq_default_randtest(fq, state, ctx);
        fq_default_clear(fq, ctx);
        fq_default_ctx_clear(ctx);

        fmpz_clear(p);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fq_default_ctx_t ctx;
        fq_default_t x, y;
        int type = 1 + n_randint(state, 5);

        fmpz_init(p);
        fmpz_set_ui(p, 7);
        fq_default_ctx_init_type(ctx, p, 1, "x", type);
        FLINT_TEST(fq_default_ctx_type(ctx) == type);
        fq_default_init(x, ctx);
        fq_default_init(y, ctx);
        fq_default_randtest(x, state, ctx);
        fq_default_sqr(y, x, ctx);

        if (type == FQ_DEFAULT_FQ)
            fq_sqr((fq_struct *) x, (const fq_struct *) x, (const fq_ctx_struct *) fq_default_ctx_inner(ctx));
        else if (type == FQ_DEFAULT_FQ_NMOD)
            fq_nmod_sqr((fq_nmod_struct *) x, (const fq_nmod_struct *) x, (const fq_nmod_ctx_struct *) fq_default_ctx_inner(ctx));
        else if (type == FQ_DEFAULT_FQ_ZECH)
            fq_zech_sqr((fq_zech_struct *) x, (const fq_zech_struct *) x, (const fq_zech_ctx_struct *) fq_default_ctx_inner(ctx));
        else if (type == FQ_DEFAULT_FMPZ_MOD)
            fmpz_mod_mul((fmpz *) x, (const fmpz *) x, (const fmpz *) x, (const fmpz_mod_ctx_struct *) fq_default_ctx_inner(ctx));
        else
            ((ulong *) x)[0] = nmod_mul(((ulong *) x)[0], ((ulong *) x)[0], *((const nmod_t *) fq_default_ctx_inner(ctx)));

        FLINT_TEST(fq_default_equal(x, y, ctx));
        fq_default_clear(x, ctx);
        fq_default_clear(y, ctx);
        fq_default_ctx_clear(ctx);
        fmpz_clear(p);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_t fq;

        fq_default_ctx_init_randtest(ctx, state);
        fq_default_init(fq, ctx);
        fq_default_randtest(fq, state, ctx);
        fq_default_clear(fq, ctx);
        fq_default_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
