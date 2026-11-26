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
#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"

TEST_FUNCTION_START(acb_ode_lc_roots, state)
{
    gr_ctx_t ZZ, QQ, QQbar;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_complex_qqbar(QQbar);

    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t Rnd, CC, Pol, Dop;

        slong prec = 2 + n_randint(state, 200);
        gr_ctx_init_complex_acb(CC, prec);
        gr_ctx_init_random(Rnd, state);

        gr_ctx_struct * rings[] = { ZZ, QQ, QQbar, CC, Rnd };
        slong r = n_randint(state, 5);
        gr_ctx_struct * Cst = rings[r];
        int expect_success = (r < 3);

        gr_ctx_init_gr_poly(Pol, Cst);
        gr_ore_poly_ctx_init(Dop, Pol, 0, n_randint(state, 2)
                                          ? ORE_ALGEBRA_EULER_DERIVATIVE
                                          : ORE_ALGEBRA_DERIVATIVE);

        gr_ore_poly_t dop;
        acb_poly_t lc, ref;
        gr_ore_poly_init(dop, Dop);
        acb_poly_init(lc);
        acb_poly_init(ref);

        slong len = n_randint(state, (Cst == QQbar) ? 3 : 4);
        GR_MUST_SUCCEED(gr_ore_poly_randtest(dop, state, len, Dop));
        int include_zero = n_randint(state, 2);

        if (dop->length == 0)
        {
            status = acb_ode_lc_roots(NULL, dop, Dop, include_zero, CC);
            FLINT_TEST(status == GR_DOMAIN);
            goto cleanup;
        }

        gr_poly_struct * dlc = gr_ore_poly_coeff_ptr(dop, dop->length - 1, Dop);
        slong deg = dlc->length - 1;
        acb_ptr lcrt = _acb_vec_init(deg);

        status |= acb_ode_lc_roots(lcrt, dop, Dop, include_zero, CC);

        if (expect_success && status != GR_SUCCESS)
        {
            flint_printf("FAIL (unexpected failure)\n\n");
            flint_printf("dop = %{gr}\n", dop, Dop);
            flint_abort();
        }

        status |= gr_poly_set_gr_poly_other((gr_poly_struct *) lc, dlc, Cst,
                                            CC);

        if (status == GR_SUCCESS)
        {
            slong val = 0;
            if (!include_zero)
                while (val < deg && acb_is_zero(lc->coeffs + val))
                    val++;

            if (val <= deg
                && gr_is_zero(gr_poly_coeff_ptr(dlc, val, Cst), Cst) == T_FALSE)
            {
                acb_poly_product_roots(ref, lcrt, deg - val, prec);
                acb_poly_shift_left(ref, ref, val);
                acb_poly_scalar_mul(ref, ref, lc->coeffs + deg, prec);

                if (!acb_poly_overlaps(ref, lc))
                {
                    flint_printf("FAIL\n\n");
                    flint_printf("dop = %{gr}\n", dop, Dop);
                    flint_printf("include_zero = %d\n", include_zero);
                    flint_printf("lcrt = %{acb*}\n", lcrt, deg - val);
                    flint_printf("lc = %{acb_poly}\n", lc);
                    flint_printf("ref = %{acb_poly}\n", ref);
                    flint_abort();
                }
            }
        }

        _acb_vec_clear(lcrt, deg);

cleanup:
        acb_poly_clear(ref);
        acb_poly_clear(lc);
        gr_ore_poly_clear(dop, Dop);
        gr_ctx_clear(Dop);
        gr_ctx_clear(Pol);
        gr_ctx_clear(Rnd);
        gr_ctx_clear(CC);
    }

    gr_ctx_clear(QQbar);
    gr_ctx_clear(QQ);
    gr_ctx_clear(ZZ);

    TEST_FUNCTION_END(state);
}
