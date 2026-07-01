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

TEST_FUNCTION_START(padic_nmod_div, state)
{
    /* Check that division coincides with the padic_t one */
    for (int i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong p = n_randprime(state, FLINT_BITS / 4, 1);
        ulong max_val = FLINT_BITS / n_clog(p, 2);

        fmpz_t p_fmpz;
        fmpz_init(p_fmpz);
        fmpz_set_ui(p_fmpz, p);

        padic_ctx_t ctx_padic;
        gr_ctx_t ctx_nmod;

        padic_ctx_init(ctx_padic, p_fmpz, 0, max_val + 1, PADIC_VAL_UNIT);
        padic_nmod_ctx_init(ctx_nmod, p, max_val);

        padic_t a, b, a_temp, b_temp, res;
        padic_nmod_t c, d, res_float;

        padic_init2(a, max_val / 3);
        padic_init2(b, max_val / 3);
        padic_init2(a_temp, max_val);
        padic_init2(b_temp, max_val);
        padic_init2(res, max_val);
        padic_nmod_init(c, ctx_nmod);
        padic_nmod_init(d, ctx_nmod);
        padic_nmod_init(res_float, ctx_nmod);

        padic_randtest(a, state, ctx_padic);
        padic_randtest_not_zero(b, state, ctx_padic);
        padic_nmod_set_ui(c, fmpz_get_ui(padic_unit(a)), ctx_nmod);
        c->v = padic_get_val(a);
        _padic_nmod_canonicalise(c, ctx_nmod);
        padic_nmod_set_ui(d, fmpz_get_ui(padic_unit(b)), ctx_nmod);
        d->v = padic_get_val(b);
        _padic_nmod_canonicalise(d, ctx_nmod);
        padic_set(a_temp, a, ctx_padic);
        padic_set(b_temp, b, ctx_padic);

        padic_div(res, a_temp, b_temp, ctx_padic);
        padic_nmod_div(res_float, c, d, ctx_nmod);

        if (!fmpz_equal_ui(padic_unit(a), c->u))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = ");
            padic_print(a, ctx_padic);
            flint_printf("\nc = ");
            padic_nmod_println(c, ctx_nmod);
            flint_printf("p = %wd\n", p);
            flint_printf("precision = %wd\n", max_val);
            flint_printf("p^precision = %wd\n",
                         PADIC_NMOD_CTX_PN_MOD(ctx_nmod).n);
            fflush(stdout);
            flint_abort();
        }
        if (!fmpz_equal_ui(padic_unit(b), d->u))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("b = ");
            padic_print(b, ctx_padic);
            flint_printf("\nd = ");
            padic_nmod_println(d, ctx_nmod);
            flint_printf("p = %wu\n", p);
            flint_printf("precision = %wu\n", max_val);
            flint_printf("p^precision = %wu\n",
                         PADIC_NMOD_CTX_PN_MOD(ctx_nmod).n);
            fflush(stdout);
            flint_abort();
        }

        ulong inv_prec;

        if (c->v < d->v)
            inv_prec = max_val;
        else if (c->v > d->v + (signed) max_val)
            inv_prec = 0;
        else
            inv_prec = max_val - (c->v - d->v);

        if ((fmpz_fdiv_ui(padic_unit(res), PADIC_NMOD_CTX_PN_MOD(ctx_nmod).n)
             != res_float->u % fmpz_get_ui(ctx_padic->pow + inv_prec))
            || ((padic_get_val(res) != res_float->v)
                && (!padic_is_zero(res)
                    || (padic_nmod_is_zero(res_float, ctx_nmod) != T_TRUE))))
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = ");
            padic_print(a, ctx_padic);
            flint_printf("\nb = ");
            padic_print(b, ctx_padic);
            flint_printf("\nc = ");
            padic_nmod_println(c, ctx_nmod);
            flint_printf("d = ");
            padic_nmod_println(d, ctx_nmod);
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
        padic_clear(b);
        padic_clear(res);

        padic_ctx_clear(ctx_padic);
        padic_nmod_ctx_clear(ctx_nmod);
    }

    TEST_FUNCTION_END(state);
}
