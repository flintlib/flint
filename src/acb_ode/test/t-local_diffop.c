/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Code (claude-opus-5-5) */

#include "test_helpers.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"

static slong
valuation(const gr_poly_t f, gr_ctx_t Cst)
{
    slong v = 0;
    while (v < f->length
           && gr_is_zero(gr_poly_coeff_srcptr(f, v, Cst), Cst) == T_TRUE)
        v++;
    return v;
}

TEST_FUNCTION_START(acb_ode_local_diffop, state)
{
    gr_ctx_t ZZ, QQ, ZZi, QQbar;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_fmpzi(ZZi);
    gr_ctx_init_complex_qqbar(QQbar);

    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        const char * fail = NULL;
        gr_ctx_t Pol, Dop, Pol1, Euler;

        gr_ctx_struct * rings[] = { ZZ, QQ, ZZi, QQbar };
        gr_ctx_struct * Cst = rings[n_randint(state, 4)];

        slong maxlen = (Cst == QQbar) ? 2 : 5;

        gr_ctx_init_gr_poly(Pol, Cst);
        gr_ore_poly_ctx_init(Dop, Pol, 0, n_randint(state, 2)
                                          ? ORE_ALGEBRA_EULER_DERIVATIVE
                                          : ORE_ALGEBRA_DERIVATIVE);
        acb_ode_local_diffop_ctx_init(Pol1, Euler, Dop, Cst);

        gr_ore_poly_t dop, ldop;
        gr_poly_t f, g, lg;
        gr_ore_poly_init(dop, Dop);
        gr_ore_poly_init(ldop, Euler);
        gr_poly_init(f, Cst);
        gr_poly_init(g, Cst);
        gr_poly_init(lg, Cst);

        status |= gr_ore_poly_randtest(dop, state,
                                       n_randint(state, maxlen + 1), Dop);
        if (n_randint(state, 2))
        {
            /* multiply by a power of x */
            gr_poly_t xk;
            gr_poly_init(xk, Cst);
            status |= gr_poly_set_coeff_si(xk, n_randint(state, 4), 1, Cst);
            status |= gr_ore_poly_other_mul(dop, xk, Pol, dop, Dop);
            gr_poly_clear(xk, Cst);
        }
        GR_MUST_SUCCEED(status);

        status |= acb_ode_local_diffop(ldop, dop, Dop, Euler);

        if (status != GR_SUCCESS)
        {
            fail = "unexpected failure";
            goto end;
        }

        int val0 = 0;
        for (slong i = 0; i < ldop->length; i++)
        {
            gr_poly_struct * c = gr_ore_poly_coeff_ptr(ldop, i, Euler);
            val0 |= (c->length > 0 && valuation(c, Cst) == 0);
        }
        if (dop->length > 0 && !val0)
        {
            fail = "valuation";
            goto end;
        }

        /* ldop(f) = x^s·dop(f) with s independent of f */
        slong s = WORD_MIN;
        for (int j = 0; j < 2; j++)
        {
            status |= gr_poly_randtest(f, state,
                                       n_randint(state, maxlen + 1), Cst);
            status |= gr_ore_poly_apply(g, dop, f, Dop);
            status |= gr_ore_poly_apply(lg, ldop, f, Euler);
            GR_MUST_SUCCEED(status);

            if ((g->length == 0) != (lg->length == 0))
            {
                fail = "apply (zero)";
                goto end;
            }
            if (g->length == 0)
                continue;

            slong v = valuation(g, Cst);
            slong lv = valuation(lg, Cst);
            if (s == WORD_MIN)
                s = lv - v;

            status |= gr_poly_shift_left(g, g, lv, Cst);
            status |= gr_poly_shift_left(lg, lg, v, Cst);
            if (lv - v != s || gr_poly_equal(g, lg, Cst) != T_TRUE)
            {
                fail = "apply";
                goto end;
            }
        }

end:
        if (fail != NULL)
        {
            flint_printf("FAIL (%s)\n\n", fail);
            flint_printf("Dop = %{gr_ctx}\n", Dop);
            flint_printf("dop = %{gr}\n", dop, Dop);
            flint_printf("ldop = %{gr}\n", ldop, Euler);
            flint_printf("f = %{gr}\n", f, Pol);
            flint_abort();
        }

        gr_poly_clear(lg, Cst);
        gr_poly_clear(g, Cst);
        gr_poly_clear(f, Cst);
        gr_ore_poly_clear(ldop, Euler);
        gr_ore_poly_clear(dop, Dop);
        gr_ctx_clear(Euler);
        gr_ctx_clear(Pol1);
        gr_ctx_clear(Dop);
        gr_ctx_clear(Pol);
    }

    gr_ctx_clear(QQbar);
    gr_ctx_clear(ZZi);
    gr_ctx_clear(QQ);
    gr_ctx_clear(ZZ);

    TEST_FUNCTION_END(state);
}
