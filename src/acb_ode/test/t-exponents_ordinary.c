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

TEST_FUNCTION_START(acb_ode_exponents_ordinary, state)
{
    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        acb_ode_exponents_t expos;
        acb_poly_t ind, ref;
        acb_ptr roots;
        slong order, prec;

        acb_ode_exponents_init(expos);
        acb_poly_init(ind);
        acb_poly_init(ref);

        order = n_randint(state, 8);
        prec = 30 + n_randint(state, 100);

        if (n_randint(state, 2))
            /* test overwriting an existing object */
            acb_ode_exponents_randtest(expos, state, 1 + n_randint(state, 6),
                                       10, prec, 3);
        acb_ode_exponents_ordinary(expos, order);

        FLINT_TEST(expos->ngroups == 1);
        FLINT_TEST(acb_ode_exponents_length(expos) == order);

        acb_ode_group_struct * grp = expos->groups + 0;
        FLINT_TEST(acb_is_zero(grp->leader));
        FLINT_TEST(acb_ode_group_length(grp) == order);
        FLINT_TEST(acb_ode_group_multiplicity(grp, -1) == 0);
        for (slong n = 0; n < order; n++)
        {
            FLINT_TEST(acb_ode_group_multiplicity(grp, n) == 1);
            FLINT_TEST(acb_ode_group_nlogs(grp, n) == n + 1);
        }
        FLINT_TEST(acb_ode_group_multiplicity(grp, order) == 0);

        /* the indicial polynomial is prod_{k=0}^{order-1} (x - k) */
        roots = _acb_vec_init(order);
        for (slong k = 0; k < order; k++)
            acb_set_si(roots + k, k);
        acb_poly_product_roots(ref, roots, order, prec);

        acb_ode_indicial_polynomial_from_exponents(ind, expos, prec);

        if (!acb_poly_overlaps(ind, ref))
        {
            flint_printf("FAIL (ordinary indicial polynomial)\n\n");
            flint_printf("order = %wd\n", order);
            flint_printf("ind = %{acb_poly}\n", ind);
            flint_printf("ref = %{acb_poly}\n", ref);
            flint_abort();
        }

        _acb_vec_clear(roots, order);
        acb_poly_clear(ref);
        acb_poly_clear(ind);
        acb_ode_exponents_clear(expos);
    }

    TEST_FUNCTION_END(state);
}
