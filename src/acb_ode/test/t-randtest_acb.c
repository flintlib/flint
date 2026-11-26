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

TEST_FUNCTION_START(acb_ode_randtest_acb, state)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        acb_ode_exponents_t expos;
        acb_poly_t ind, ref;
        acb_t c, c0;
        const char * fail = NULL;

        slong prec = 2 + n_randint(state, 200);
        slong len = 1 + n_randint(state, 6);
        slong lcdeg = n_randint(state, 5);
        slong clen = 1 + n_randint(state, 6);
        slong disp = n_randint(state, 20);

        acb_poly_struct * dop = _acb_poly_vec_init(len);
        acb_ptr lcroots = _acb_vec_init(lcdeg);
        acb_ode_exponents_init(expos);
        acb_poly_init(ind);
        acb_poly_init(ref);
        acb_init(c);
        acb_init(c0);

        acb_ode_randtest_acb(dop, lcroots, expos, state, disp, lcdeg, clen,
                             len, prec);

        acb_poly_struct * lc = dop + len - 1;

        if (acb_ode_exponents_length(expos) != len - 1)
        {
            fail = "number of exponents";
            goto end;
        }

        for (slong g = 0; g < expos->ngroups; g++)
        {
            acb_ode_group_struct * grp = expos->groups + g;
            if (grp->shifts[grp->nshifts - 1].n > disp)
            {
                fail = "dispersion";
                goto end;
            }
        }

        if (acb_poly_length(lc) != lcdeg + 1)
        {
            fail = "degree of the leading coefficient";
            goto end;
        }

        for (slong i = 0; i < len - 1; i++)
        {
            if (acb_poly_length(dop + i) > clen)
            {
                fail = "length of the coefficients";
                goto end;
            }
        }

        acb_poly_product_roots(ref, lcroots, lcdeg, prec);
        acb_poly_scalar_mul(ref, ref, lc->coeffs + lcdeg, prec);
        if (!acb_poly_overlaps(ref, lc))
        {
            fail = "roots of the leading coefficient";
            goto end;
        }

        acb_ode_indicial_polynomial_from_exponents(ind, expos, prec);
        for (slong i = 0; i < len; i++)
        {
            acb_mul(c, lc->coeffs, ind->coeffs + i, prec);
            acb_poly_get_coeff_acb(c0, dop + i, 0);
            if (!acb_overlaps(c, c0))
            {
                fail = "indicial polynomial";
                goto end;
            }
        }

end:
        if (fail != NULL)
        {
            flint_printf("FAIL (%s)\n\n", fail);
            flint_printf("len = %wd, lcdeg = %wd, clen = %wd, disp = %wd, "
                         "prec = %wd\n", len, lcdeg, clen, disp, prec);
            for (slong i = 0; i < len; i++)
                flint_printf("dop[%wd] = %{acb_poly}\n", i, dop + i);
            flint_printf("lcroots = %{acb*}\n", lcroots, lcdeg);
            flint_printf("expos = ");
            acb_ode_exponents_println(expos);
            flint_abort();
        }

        acb_clear(c0);
        acb_clear(c);
        acb_poly_clear(ref);
        acb_poly_clear(ind);
        acb_ode_exponents_clear(expos);
        _acb_vec_clear(lcroots, lcdeg);
        _acb_poly_vec_clear(dop, len);
    }

    TEST_FUNCTION_END(state);
}
