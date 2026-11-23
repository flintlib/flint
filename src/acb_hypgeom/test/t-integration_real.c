/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

/* acb_hypgeom_{1f1,u,2f1}_integration with real input take the integral on
   the real line first: the result overlaps the one along the paths, and in
   the domain of the real integral it is as accurate as that one */
TEST_FUNCTION_START(acb_hypgeom_integration_real, state)
{
    slong iter;

    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z, r1, r2;
        slong prec, acc1, acc2;
        double size;
        int which, regularized, domain;

        acb_init(a); acb_init(b); acb_init(c); acb_init(z);
        acb_init(r1); acb_init(r2);

        prec = 30 + n_randint(state, 200);
        which = n_randint(state, 3);
        regularized = n_randint(state, 2);
        size = (n_randint(state, 2) == 0) ? 20.0 : 1000.0;

        /* parameters which are mostly in the domain of the real integral,
           sometimes not */
        acb_set_d(a, size * (n_randint(state, 1000) / 1000.0) + ((n_randint(state, 4) == 0) ? -size / 2 : 1.0));
        acb_set_d(b, size * (n_randint(state, 1000) / 1000.0) + 1.0);
        acb_add(b, b, a, prec);
        acb_add_ui(b, b, 1, prec);
        acb_set_d(c, size * (n_randint(state, 1000) / 1000.0) + 1.0);
        acb_add(c, c, b, prec);
        if (which == 2)
            acb_set_d(z, 1.8 * (n_randint(state, 1000) / 1000.0) - 1.0);
        else
            acb_set_d(z, size * (n_randint(state, 2000) / 1000.0 - 1.0) * ((which == 1) ? 0.5 : 1.0) + ((which == 1) ? size : 0.0));

        if (n_randint(state, 4) == 0)
        {
            mag_set_ui_2exp_si(arb_radref(acb_realref(a)), 1, -prec + 8);
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -prec + 8);
        }

        if (which == 0)
        {
            acb_hypgeom_1f1_integration(r1, a, b, z, regularized, prec);
            _acb_hypgeom_1f1_integration_paths(r2, a, b, z, regularized, prec);
        }
        else if (which == 1)
        {
            acb_hypgeom_u_integration(r1, a, b, z, prec);
            _acb_hypgeom_u_integration_paths(r2, a, b, z, prec);
        }
        else
        {
            acb_hypgeom_2f1_integration(r1, a, b, c, z, regularized, prec);
            _acb_hypgeom_2f1_integration_paths(r2, a, b, c, z, regularized, prec);
        }

        if (!acb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap (%d)\n\n", which);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
            flint_abort();
        }

        if (!acb_is_real(r1) && acb_is_finite(r1))
        {
            flint_printf("FAIL: not real (%d)\n\n", which);
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
            flint_abort();
        }

        /* the domain of the real integral (U: a >= 1 and z > 0) */
        domain = arf_cmp_si(arb_midref(acb_realref(a)), 1) >= 0 && arb_is_exact(acb_realref(a))
            && (which == 1 ? arb_is_positive(acb_realref(z))
                : (which == 0 ? 1 : arf_cmp_si(arb_midref(acb_realref(z)), 1) < 0));
        if (domain)
        {
            acc1 = acb_rel_accuracy_bits(r1);
            acc2 = acb_rel_accuracy_bits(r2);

            if (acc1 < acc2 - 2)
            {
                flint_printf("FAIL: accuracy (%d)\n\n", which);
                flint_printf("prec = %wd, acc1 = %wd, acc2 = %wd\n\n", prec, acc1, acc2);
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z);
        acb_clear(r1); acb_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
