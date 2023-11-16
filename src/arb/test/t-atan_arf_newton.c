/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "arb.h"

#ifndef arb_atan_arf_via_mpfr
#define arb_atan_arf_via_mpfr arb_atan_arf_via_mpfr
void
arb_atan_arf_via_mpfr(arb_t z, const arf_t x, slong prec)
{
    mpfr_t t, u;
    int exact;

    mpfr_init2(t, 2 + arf_bits(x));
    mpfr_init2(u, prec);

    mpfr_set_emin(MPFR_EMIN_MIN);
    mpfr_set_emax(MPFR_EMAX_MAX);

    arf_get_mpfr(t, x, MPFR_RNDD);
    exact = (mpfr_atan(u, t, MPFR_RNDD) == 0);

    arf_set_mpfr(arb_midref(z), u);
    if (!exact)
        arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);

    mpfr_clear(t);
    mpfr_clear(u);
}
#endif

TEST_FUNCTION_START(arb_atan_arf_newton, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z;
        slong prec, prec2;

        arb_init(x);
        arb_init(y);
        arb_init(z);

        prec = 2 + n_randint(state, 8000);

        if (n_randint(state, 100) == 0)
            flint_set_num_threads(1 + n_randint(state, 4));
        else
            flint_set_num_threads(1);

        arb_randtest(x, state, 1 + n_randint(state, 8000), 3);
        mag_zero(arb_radref(x));

        if (n_randint(state, 2))
            arb_mul_2exp_si(x, x, n_randint(state, 1.5 * prec));
        else
            arb_mul_2exp_si(x, x, -n_randint(state, 1.5 * prec));

        if (!arf_is_special(arb_midref(x)))
            prec2 = prec + 100 + 2 * FLINT_MAX(0, -ARF_EXP(arb_midref(x)));
        else
            prec2 = prec + 100;

        arb_atan_arf_via_mpfr(y, arb_midref(x), prec2);
        arb_atan_arf_newton(z, arb_midref(x), prec);

        if (!arb_contains(z, y))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("prec = %wd\n\n", prec);
            flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 50); flint_printf("\n\n");
            flint_printf("z = "); arb_printd(z, 50); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_rel_accuracy_bits(z) < prec - 2)
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec = %wd,  acc = %wd\n\n", prec, arb_rel_accuracy_bits(z));
            flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 50); flint_printf("\n\n");
            flint_printf("z = "); arb_printd(z, 50); flint_printf("\n\n");
            flint_abort();
        }

        arb_atan_arf_newton(x, arb_midref(x), prec);

        if (!arb_overlaps(x, z))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
    }

    TEST_FUNCTION_END(state);
}
