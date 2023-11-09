/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_default.h"

TEST_FUNCTION_START(fq_default_ctx_init_modulus_nmod, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_t fq;
        ulong p;
        slong len;
        nmod_poly_t mod;

        p = 3;

        nmod_poly_init(mod, p);

        len = n_randint(state, 16) + 2;
        nmod_poly_randtest_irreducible(mod, state, len);

        fq_default_ctx_init_modulus_nmod(ctx, mod, "x");

        fq_default_init(fq, ctx);

        fq_default_randtest(fq, state, ctx);

        fq_default_clear(fq, ctx);

        fq_default_ctx_clear(ctx);
    }

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_t fq;
        ulong p;
        slong len;
        nmod_poly_t mod;

        p = 3;

        nmod_poly_init(mod, p);

        len = n_randint(state, 16) + 2;
        nmod_poly_randtest_irreducible(mod, state, len);

        fq_default_ctx_init_modulus_nmod_type(ctx, mod, "x", 3);

        fq_default_init(fq, ctx);

        fq_default_randtest(fq, state, ctx);

        fq_default_clear(fq, ctx);

        fq_default_ctx_clear(ctx);
    }

    {
        fq_default_ctx_t ctx;
        fq_default_t fq;
        ulong p;
        int result;
        nmod_poly_t mod;
        fmpz_mod_ctx_t mod_ctx;
        fmpz_mod_poly_t mod2, mod3;
        fmpz_t pp;

        p = 3;

        nmod_poly_init(mod, p);

        nmod_poly_fit_length(mod, 2);

        nmod_poly_set_coeff_ui(mod, 0, 2);
        nmod_poly_set_coeff_ui(mod, 1, 1);

        fq_default_ctx_init_modulus_nmod(ctx, mod, "x");

        fmpz_init(pp);
        fmpz_set_ui(pp, 3);
        fmpz_mod_ctx_init(mod_ctx, pp);
        fmpz_mod_poly_init(mod2, mod_ctx);
        fmpz_mod_poly_init(mod3, mod_ctx);
        fmpz_mod_poly_set_coeff_ui(mod3, 0, 2, mod_ctx);
        fmpz_mod_poly_set_coeff_ui(mod3, 1, 1, mod_ctx);
        fq_default_ctx_modulus(mod2, ctx);
        result = fmpz_mod_poly_equal(mod2, mod3, mod_ctx);

        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(mod2, mod_ctx); flint_printf("\n");
            fmpz_mod_poly_print(mod3, mod_ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fq_default_clear(fq, ctx);
        fq_default_ctx_clear(ctx);
        fmpz_mod_poly_clear(mod2, mod_ctx);
        fmpz_mod_poly_clear(mod3, mod_ctx);
        fmpz_mod_ctx_clear(mod_ctx);
        fmpz_clear(pp);
        nmod_poly_clear(mod);
    }

    {
        fq_default_ctx_t ctx;
        fq_default_t fq;
        ulong p;
        int result;
        nmod_poly_t mod;
        fmpz_mod_ctx_t mod_ctx;
        fmpz_mod_poly_t mod2, mod3;
        fmpz_t pp;

        p = 3;

        nmod_poly_init(mod, p);

        nmod_poly_fit_length(mod, 2);

        nmod_poly_set_coeff_ui(mod, 0, 2);
        nmod_poly_set_coeff_ui(mod, 1, 1);

        fq_default_ctx_init_modulus_nmod_type(ctx, mod, "x", FQ_DEFAULT_FMPZ_MOD);

        fmpz_init(pp);
        fmpz_set_ui(pp, 3);
        fmpz_mod_ctx_init(mod_ctx, pp);
        fmpz_mod_poly_init(mod2, mod_ctx);
        fmpz_mod_poly_init(mod3, mod_ctx);
        fmpz_mod_poly_set_coeff_ui(mod3, 0, 2, mod_ctx);
        fmpz_mod_poly_set_coeff_ui(mod3, 1, 1, mod_ctx);
        fq_default_ctx_modulus(mod2, ctx);
        result = fmpz_mod_poly_equal(mod2, mod3, mod_ctx);

        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_mod_poly_print(mod2, mod_ctx); flint_printf("\n");
            fmpz_mod_poly_print(mod3, mod_ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fq_default_clear(fq, ctx);
        fq_default_ctx_clear(ctx);
        fmpz_mod_poly_clear(mod2, mod_ctx);
        fmpz_mod_poly_clear(mod3, mod_ctx);
        fmpz_mod_ctx_clear(mod_ctx);
        fmpz_clear(pp);
        nmod_poly_clear(mod);
    }

    TEST_FUNCTION_END(state);
}
