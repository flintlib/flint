/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Generated using Claude Opus 4.8 */

#include "test_helpers.h"
#include "acb.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "acb_types.h"
#include "gr.h"
#include "gr_ore_poly.h"

TEST_FUNCTION_START(acb_ode_exponents, state)
{
    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t CC, Pol;
        gr_ore_poly_ctx_t Dop;
        gr_ore_poly_t dop;
        acb_ode_exponents_t expos_ref, expos;
        acb_poly_t ind_ref, ind;
        acb_ptr lcroots;
        slong prec, dop_len, clen, lcdeg, disp;
        int status;

        prec = 2 + n_randint(state, 100);

        gr_ctx_init_complex_acb(CC, prec);
        gr_ctx_init_gr_poly(Pol, CC);
        gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);

        gr_ore_poly_init(dop, Dop);
        acb_ode_exponents_init(expos_ref);
        acb_ode_exponents_init(expos);
        acb_poly_init(ind_ref);
        acb_poly_init(ind);

        dop_len = n_randint(state, 5);
        clen = 1 + n_randint(state, 5);
        lcdeg = n_randint(state, clen);
        disp = n_randint(state, 10);

        lcroots = _acb_vec_init(FLINT_MAX(lcdeg, 1));

        gr_ore_poly_fit_length(dop, dop_len, Dop);
        acb_ode_randtest_acb(dop->coeffs, lcroots, expos_ref, state, disp,
                             lcdeg, clen, dop_len, prec);
        _gr_ore_poly_set_length(dop, dop_len, Dop);

        status = acb_ode_exponents(expos, dop, Dop, CC);

        if (status == GR_SUCCESS)
        {
            acb_ode_indicial_polynomial_from_exponents(ind_ref, expos_ref,
                                                       prec);
            acb_ode_indicial_polynomial_from_exponents(ind, expos, prec);

            if (acb_ode_exponents_length(expos)
                    != acb_ode_exponents_length(expos_ref)
                || !acb_poly_overlaps(ind, ind_ref))
            {
                flint_printf("FAIL\n\n");
                flint_printf("dop = %{gr}\n", dop, Dop);
                flint_printf("expos_ref = ");
                acb_ode_exponents_println(expos_ref);
                flint_printf("expos     = ");
                acb_ode_exponents_println(expos);
                flint_printf("ind_ref = %{acb_poly}\n", ind_ref);
                flint_printf("ind     = %{acb_poly}\n", ind);
                flint_abort();
            }
        }

        _acb_vec_clear(lcroots, FLINT_MAX(lcdeg, 1));
        acb_poly_clear(ind);
        acb_poly_clear(ind_ref);
        acb_ode_exponents_clear(expos);
        acb_ode_exponents_clear(expos_ref);
        gr_ore_poly_clear(dop, Dop);
        gr_ctx_clear(Dop);
        gr_ctx_clear(Pol);
        gr_ctx_clear(CC);
    }

    TEST_FUNCTION_END(state);
}
