/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "acb_modular.h"

TEST_FUNCTION_START(acb_modular_fundamental_domain_approx, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t x, y;
        psl2z_t g;
        arf_t one_minus_eps, tol;
        acb_t z, w, w2;
        arb_t t;
        slong prec;

        fmpq_init(x);
        fmpq_init(y);

        psl2z_init(g);
        acb_init(z);
        acb_init(w);
        acb_init(w2);
        arf_init(one_minus_eps);
        arf_init(tol);
        arb_init(t);

        /* pick an exact point in the upper half plane */
        fmpq_randtest(x, state, 1 + n_randint(state, 500));
        do {
            fmpq_randtest(y, state, 1 + n_randint(state, 500));
        } while (fmpz_sgn(fmpq_numref(y)) <= 0);

        /* pick a tolerance */
        arf_set_ui_2exp_si(tol, 1, -(slong) n_randint(state, 500));

        /* now increase the precision until convergence */
        for (prec = 32; ; prec *= 2)
        {
            if (prec > 16384)
            {
                flint_printf("FAIL (no convergence)\n");
                flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
                flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 50); flint_printf("\n\n");
                flint_printf("w = "); acb_printd(w, 50); flint_printf("\n\n");
                flint_printf("w2 = "); acb_printd(w2, 50); flint_printf("\n\n");
                flint_printf("g = "); psl2z_print(g); flint_printf("\n\n");
                flint_abort();
            }

            arb_set_fmpq(acb_realref(z), x, prec);
            arb_set_fmpq(acb_imagref(z), y, prec);

            arf_set_ui_2exp_si(one_minus_eps, 1, -prec / 4);
            arf_sub_ui(one_minus_eps, one_minus_eps, 1, prec, ARF_RND_DOWN);
            arf_neg(one_minus_eps, one_minus_eps);

            acb_modular_fundamental_domain_approx(w, g, z, one_minus_eps, prec);
            acb_modular_transform(w2, g, z, prec);

            if (!psl2z_is_correct(g) || !acb_overlaps(w, w2))
            {
                flint_printf("FAIL (incorrect transformation)\n");
                flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
                flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 50); flint_printf("\n\n");
                flint_printf("w = "); acb_printd(w, 50); flint_printf("\n\n");
                flint_printf("w2 = "); acb_printd(w2, 50); flint_printf("\n\n");
                flint_printf("g = "); psl2z_print(g); flint_printf("\n\n");
                flint_abort();
            }

            /* success */
            if (acb_modular_is_in_fundamental_domain(w, tol, prec))
                break;
        }

        /* check that g^(-1) * w contains x+yi */
        psl2z_inv(g, g);
        acb_modular_transform(w2, g, w, 2 + n_randint(state, 1000));
        if (!arb_contains_fmpq(acb_realref(w2), x) ||
            !arb_contains_fmpq(acb_imagref(w2), y))
        {
            flint_printf("FAIL (inverse containment)\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 50); flint_printf("\n\n");
            flint_printf("w = "); acb_printd(w, 50); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 50); flint_printf("\n\n");
            flint_printf("g = "); psl2z_print(g); flint_printf("\n\n");
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);

        psl2z_clear(g);
        acb_clear(z);
        acb_clear(w);
        acb_clear(w2);
        arf_clear(one_minus_eps);
        arf_clear(tol);
        arb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
