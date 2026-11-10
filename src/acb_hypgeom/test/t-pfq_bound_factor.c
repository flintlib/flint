/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <math.h>
#include "acb.h"
#include "acb_hypgeom.h"

/* Tail bounds with a lower parameter b near the negative real axis,
   truncating before the Pochhammer symbol (b)_k has passed the negative
   integers: compare with summing far enough past them. */
TEST_FUNCTION_START(acb_hypgeom_pfq_bound_factor, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_struct a[2], b[4];
        acb_t z, r1, r2;
        slong p, q, N, n1, n2, prec, i;
        mag_t C;

        for (i = 0; i < 2; i++)
            acb_init(a + i);
        for (i = 0; i < 4; i++)
            acb_init(b + i);
        acb_init(z);
        acb_init(r1);
        acb_init(r2);
        mag_init(C);

        prec = 2 + n_randint(state, 200);
        p = n_randint(state, 3);
        q = p + n_randint(state, 2);
        q = FLINT_MAX(q, 1);

        for (i = 0; i < p; i++)
            acb_randtest_precise(a + i, state, 1 + n_randint(state, 200), 3);
        for (i = 0; i < q; i++)
        {
            acb_randtest_precise(b + i, state, 1 + n_randint(state, 200), 3);
            arb_abs(acb_realref(b + i), acb_realref(b + i));
            arb_add_ui(acb_realref(b + i), acb_realref(b + i), 1, prec + 30);
        }

        /* b_0 = -N + delta with 2^-40 < |delta| < 1 (roughly) */
        N = n_randint(state, 2) ? n_randint(state, 5) : n_randint(state, 300);
        acb_randtest_precise(z, state, 1 + n_randint(state, 200), 1);
        acb_mul_2exp_si(z, z, -6);
        acb_set_d_d(b + 0, 0.5 + 0.5 * n_randint(state, 1000) / 1000.0,
            n_randint(state, 2) ? 0.0 : 0.5 * n_randint(state, 1000) / 1000.0);
        acb_mul_2exp_si(b + 0, b + 0, -(slong) n_randint(state, 40));
        if (n_randint(state, 2))
            acb_neg(b + 0, b + 0);
        acb_add(b + 0, b + 0, z, prec);
        acb_sub_ui(b + 0, b + 0, N, prec + 30);

        acb_randtest_precise(z, state, 1 + n_randint(state, 200), 1);
        acb_mul_2exp_si(z, z, -(slong) n_randint(state, 4));

        /* sometimes large enough for the terms to grow again near the
           pole, or at first */
        if (n_randint(state, 4) == 0)
            acb_mul_2exp_si(z, z, n_randint(state, 7));

        /* usually include the factorial (b = 1) */
        if (n_randint(state, 4) != 0)
        {
            acb_one(b + q);
            q++;
        }

        n1 = n_randint(state, 30);
        n2 = N + 50 + 8 * (slong) mag_get_d(arb_radref(acb_realref(z)))
            + 8 * (slong) fabs(arf_get_d(arb_midref(acb_realref(z)), ARF_RND_UP))
            + 8 * (slong) fabs(arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_UP))
            + n_randint(state, 100);

        acb_hypgeom_pfq_bound_factor(C, a, p, b, q, z, n1);

        if (mag_is_finite(C))
        {
            acb_hypgeom_pfq_direct(r1, a, p, b, q, z, n1, prec);
            acb_hypgeom_pfq_direct(r2, a, p, b, q, z, n2, prec);

            if (!acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("p = %wd, q = %wd, n1 = %wd, n2 = %wd\n\n", p, q, n1, n2);
                for (i = 0; i < p; i++)
                {
                    flint_printf("a = "); acb_printd(a + i, 20); flint_printf("\n");
                }
                for (i = 0; i < q; i++)
                {
                    flint_printf("b = "); acb_printd(b + i, 20); flint_printf("\n");
                }
                flint_printf("z = "); acb_printd(z, 20); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printd(r1, 20); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 20); flint_printf("\n\n");
                flint_abort();
            }
        }

        for (i = 0; i < 2; i++)
            acb_clear(a + i);
        for (i = 0; i < 4; i++)
            acb_clear(b + i);
        acb_clear(z);
        acb_clear(r1);
        acb_clear(r2);
        mag_clear(C);
    }

    TEST_FUNCTION_END(state);
}
