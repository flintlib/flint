/*
    Copyright (C) 2026 Rubén Muñoz--Bertrand

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "padic.h"
#include "padic_nmod.h"

TEST_FUNCTION_START(padic_nmod_neg, state)
{
    /* Check that negation coincides with the padic_t one */
    for (int i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong p = n_randprime(state, FLINT_BITS / 4, 1);
        ulong max_val = FLINT_BITS / n_clog(p, 2);

        fmpz_t p_fmpz;
        fmpz_init(p_fmpz);
        fmpz_set_ui(p_fmpz, p);

        padic_ctx_t ctx_padic;
        gr_ctx_t ctx_nmod;

        padic_ctx_init(ctx_padic, p_fmpz, 0, max_val, PADIC_VAL_UNIT);
        padic_nmod_ctx_init(ctx_nmod, p, max_val);

        padic_t a, res;
        padic_nmod_t b, res_float;

        padic_init2(a, max_val / 3);
        padic_init2(res, max_val);
        padic_nmod_init(b, ctx_nmod);
        padic_nmod_init(res_float, ctx_nmod);

        padic_randtest_not_zero(a, state, ctx_padic);
        padic_nmod_set_ui(b, fmpz_get_ui(padic_unit(a)), ctx_nmod);
        b->v = padic_get_val(a);
        _padic_nmod_canonicalise(b, ctx_nmod);

        padic_neg(res, a, ctx_padic);
        padic_nmod_neg(res_float, b, ctx_nmod);

        if (!fmpz_equal_ui(padic_unit(a), b->u))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = ");
            padic_print(a, ctx_padic);
            flint_printf("\nb = ");
            padic_nmod_println(b, ctx_nmod);
            flint_printf("p = %wu\n", p);
            flint_printf("precision = %wu\n", max_val);
            flint_printf("p^precision = %wu\n",
                         PADIC_NMOD_CTX_PN_MOD(ctx_nmod).n);
            fflush(stdout);
            flint_abort();
        }
        if ((fmpz_fdiv_ui(padic_unit(res), PADIC_NMOD_CTX_PN_MOD(ctx_nmod).n)
             != res_float->u) || (padic_get_val(res) != res_float->v))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = ");
            padic_print(a, ctx_padic);
            flint_printf("\nb = ");
            padic_nmod_println(b, ctx_nmod);
            flint_printf("res = ");
            padic_print(res, ctx_padic);
            flint_printf("\nres_float = ");
            padic_nmod_println(res_float, ctx_nmod);
            flint_printf("p = %wu\n", p);
            flint_printf("precision = %wu\n", max_val);
            flint_printf("p^precision = %wu\n",
                         PADIC_NMOD_CTX_PN_MOD(ctx_nmod).n);
            fflush(stdout);
            flint_abort();
        }

        padic_clear(a);
        padic_clear(res);

        padic_ctx_clear(ctx_padic);
        padic_nmod_ctx_clear(ctx_nmod);
    }

    TEST_FUNCTION_END(state);
}
