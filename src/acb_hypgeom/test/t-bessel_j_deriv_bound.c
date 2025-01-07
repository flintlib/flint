/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

static void
_acb_randtest_inner(acb_t z2, flint_rand_t state, const acb_t z1)
{
    acb_zero(z2);

    arf_set_mag(arb_midref(acb_realref(z2)), arb_radref(acb_realref(z1)));
    arf_set_mag(arb_midref(acb_imagref(z2)), arb_radref(acb_imagref(z1)));

    switch (n_randint(state, 5))
    {
        case 0:
            arf_add(arb_midref(acb_realref(z2)), arb_midref(acb_realref(z1)), arb_midref(acb_realref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_add(arb_midref(acb_imagref(z2)), arb_midref(acb_imagref(z1)), arb_midref(acb_imagref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
            break;
        case 1:
            arf_add(arb_midref(acb_realref(z2)), arb_midref(acb_realref(z1)), arb_midref(acb_realref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_sub(arb_midref(acb_imagref(z2)), arb_midref(acb_imagref(z1)), arb_midref(acb_imagref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
            break;
        case 2:
            arf_sub(arb_midref(acb_realref(z2)), arb_midref(acb_realref(z1)), arb_midref(acb_realref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_add(arb_midref(acb_imagref(z2)), arb_midref(acb_imagref(z1)), arb_midref(acb_imagref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
            break;
        case 3:
            arf_sub(arb_midref(acb_realref(z2)), arb_midref(acb_realref(z1)), arb_midref(acb_realref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_sub(arb_midref(acb_imagref(z2)), arb_midref(acb_imagref(z1)), arb_midref(acb_imagref(z2)), ARF_PREC_EXACT, ARF_RND_DOWN);
            break;
        default:
            arf_set(arb_midref(acb_realref(z2)), arb_midref(acb_realref(z1)));
            arf_set(arb_midref(acb_imagref(z2)), arb_midref(acb_imagref(z1)));
    }

    if (!acb_contains(z1, z2))
        flint_abort();
}

TEST_FUNCTION_START(acb_hypgeom_bessel_j_deriv_bound, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t f, nu1, nu2, z1, z2;
        slong prec;
        mag_t B, fB;
        ulong d;

        acb_init(z1);
        acb_init(z2);
        acb_init(nu1);
        acb_init(nu2);
        acb_init(f);
        mag_init(B);
        mag_init(fB);

        acb_randtest(z1, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        arb_mul_ui(acb_realref(z1), acb_realref(z1), n_randint(state, 300), 1 + n_randint(state, 200));
        arb_mul_ui(acb_imagref(z1), acb_imagref(z1), n_randint(state, 300), 1 + n_randint(state, 200));

        _acb_randtest_inner(z2, state, z1);

        /* Currently only integer nu is implemented, so don't
           waste time testing anything else. */
        if (n_randint(state, 2))
            acb_set_si(nu1, n_randint(state, 10));
        else
            acb_set_si(nu1, n_randint(state, 100));

        if (n_randint(state, 2))
            acb_neg(nu1, nu1);
        acb_set(nu2, nu1);

        d = n_randint(state, 2);

        acb_hypgeom_bessel_j_deriv_bound(B, nu1, z1, d);

        prec = MAG_BITS + 10;

        do {
            if (d == 0)
            {
                acb_hypgeom_bessel_j(f, nu2, z2, prec);
            }
            else
            {
                /* J'_nu = (J_{nu-1} - J_{nu+1})/2 */
                acb_t t;
                acb_init(t);
                acb_sub_ui(t, nu2, 1, prec);
                acb_hypgeom_bessel_j(f, t, z2, prec);
                acb_add_ui(t, nu2, 1, prec);
                acb_hypgeom_bessel_j(t, t, z2, prec);
                acb_sub(f, f, t, prec);
                acb_mul_2exp_si(f, f, -1);
                acb_clear(t);
            }

            if (acb_rel_accuracy_bits(f) >= MAG_BITS)
                break;

            prec *= 2;
        } while (1);

        acb_get_mag_lower(fB, f);

        if (mag_cmp(fB, B) > 0)
        {
            printf("FAIL\n");

            flint_printf("d = %wu\n\n", d);
            flint_printf("z1 = "); acb_printd(z1, 20); flint_printf("\n");
            flint_printf("z2 = "); acb_printd(z2, 20); flint_printf("\n\n");
            flint_printf("nu1 = "); acb_printd(nu1, 20); flint_printf("\n");
            flint_printf("nu2 = "); acb_printd(nu2, 20); flint_printf("\n\n");

            flint_printf("B = "); mag_printd(B, 10); printf("\n");
            flint_printf("f = "); acb_printd(f, 20); flint_printf("\n");
            flint_printf("fB = "); mag_printd(fB, 10); printf("\n\n");

            flint_abort();
        }

        acb_clear(z1);
        acb_clear(z2);
        acb_clear(nu1);
        acb_clear(nu2);
        acb_clear(f);
        mag_clear(B);
        mag_clear(fB);
    }

    TEST_FUNCTION_END(state);
}
