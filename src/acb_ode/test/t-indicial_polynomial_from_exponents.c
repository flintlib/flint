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

TEST_FUNCTION_START(acb_ode_indicial_polynomial_from_exponents, state)
{
    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        acb_ode_exponents_t expos;
        acb_poly_t ind, der;
        acb_t e, val, one;
        slong len, prec, total;

        acb_ode_exponents_init(expos);
        acb_poly_init(ind);
        acb_poly_init(der);
        acb_init(e);
        acb_init(val);
        acb_init(one);
        acb_one(one);

        len = n_randint(state, 16);
        prec = 100;

        acb_ode_exponents_randtest(expos, state, len, 20, prec, 6);
        acb_ode_indicial_polynomial_from_exponents(ind, expos, prec);

        total = acb_ode_exponents_length(expos);

        if (acb_poly_degree(ind) != total)
        {
            flint_printf("FAIL (degree)\n\n");
            flint_printf("expos = "); acb_ode_exponents_println(expos);
            flint_printf("ind = %{acb_poly}\n", ind);
            flint_printf("degree = %wd, expected %wd\n",
                         acb_poly_degree(ind), total);
            flint_abort();
        }

        if (!acb_contains(acb_poly_get_coeff_ptr(ind, total), one))
        {
            flint_printf("FAIL (not monic)\n\n");
            flint_printf("expos = "); acb_ode_exponents_println(expos);
            flint_printf("leading coeff = %{acb}\n",
                         acb_poly_get_coeff_ptr(ind, total));
            flint_abort();
        }

        /* every exponent e = leader + n is a root of multiplicity mult: ind and
           its first mult - 1 derivatives vanish at e */
        for (slong g = 0; g < expos->ngroups; g++)
        {
            acb_ode_group_struct * grp = expos->groups + g;
            for (slong s = 0; s < grp->nshifts; s++)
            {
                slong n = grp->shifts[s].n;
                slong mult = grp->shifts[s].mult;

                acb_add_si(e, grp->leader, n, prec);

                acb_poly_set(der, ind);
                for (slong k = 0; k < mult; k++)
                {
                    acb_poly_evaluate(val, der, e, prec);
                    if (!acb_contains_zero(val))
                    {
                        flint_printf("FAIL (root)\n\n");
                        flint_printf("expos = ");
                        acb_ode_exponents_println(expos);
                        flint_printf("ind = %{acb_poly}\n", ind);
                        flint_printf("e = %{acb}\n", e);
                        flint_printf("derivative order k = %wd (mult = %wd)\n",
                                     k, mult);
                        flint_printf("value = %{acb}\n", val);
                        flint_abort();
                    }
                    acb_poly_derivative(der, der, prec);
                }
            }
        }

        acb_clear(one);
        acb_clear(val);
        acb_clear(e);
        acb_poly_clear(der);
        acb_poly_clear(ind);
        acb_ode_exponents_clear(expos);
    }

    TEST_FUNCTION_END(state);
}
