/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_gamma_upper_fmpq, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t a;
        arb_t z, r1, r2;
        mag_t tol1, tol2, err1, err2;
        slong N1, N2, prec1, prec2;

        mag_init(tol1);
        mag_init(tol2);
        mag_init(err1);
        mag_init(err2);

        fmpq_init(a);
        arb_init(z);
        arb_init(r1);
        arb_init(r2);

        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);

        do {
            fmpq_randtest(a, state, 5);
        } while (fmpz_is_one(fmpq_denref(a)) && fmpz_sgn(fmpq_numref(a)) <= 0);

        arb_set_ui(z, 1 + n_randint(state, 100));
        arb_div_ui(z, z, 1 + n_randint(state, 100), 200);

        mag_set_ui_2exp_si(tol1, 1, -(slong) n_randint(state, 100));
        mag_set_ui_2exp_si(tol2, 1, -(slong) n_randint(state, 100));

        N1 = _arb_hypgeom_gamma_upper_fmpq_inf_choose_N(err1, a, z, tol1);
        N2 = _arb_hypgeom_gamma_upper_fmpq_inf_choose_N(err2, a, z, tol2);

        _arb_hypgeom_gamma_upper_fmpq_inf_bsplit(r1, a, z, N1, prec1);
        arb_add_error_mag(r1, err1);

        _arb_hypgeom_gamma_upper_fmpq_inf_bsplit(r2, a, z, N2, prec2);
        arb_add_error_mag(r2, err2);

        if (!arb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); fmpq_print(a); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_printf("r1 = "); arb_printn(r1, 100, 0); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printn(r2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        fmpq_clear(a);
        arb_clear(z);
        arb_clear(r1);
        arb_clear(r2);

        mag_clear(tol1);
        mag_clear(tol2);
        mag_clear(err1);
        mag_clear(err2);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t a;
        arb_t z, r1, r2;
        mag_t tol1, tol2, err1, err2;
        slong N1, N2, prec1, prec2;

        mag_init(tol1);
        mag_init(tol2);
        mag_init(err1);
        mag_init(err2);

        fmpq_init(a);
        arb_init(z);
        arb_init(r1);
        arb_init(r2);

        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);

        do {
            fmpq_randtest(a, state, 5);
        } while (fmpz_is_one(fmpq_denref(a)) && fmpz_sgn(fmpq_numref(a)) <= 0);

        arb_set_ui(z, 1 + n_randint(state, 100));
        arb_div_ui(z, z, 1 + n_randint(state, 100), 200);

        mag_set_ui_2exp_si(tol1, 1, -(slong) n_randint(state, 100));
        mag_set_ui_2exp_si(tol2, 1, -(slong) n_randint(state, 100));

        N1 = _arb_hypgeom_gamma_lower_fmpq_0_choose_N(err1, a, z, tol1);
        N2 = _arb_hypgeom_gamma_lower_fmpq_0_choose_N(err2, a, z, tol2);

        _arb_hypgeom_gamma_lower_fmpq_0_bsplit(r1, a, z, N1, prec1);
        arb_add_error_mag(r1, err1);

        _arb_hypgeom_gamma_lower_fmpq_0_bsplit(r2, a, z, N2, prec2);
        arb_add_error_mag(r2, err2);

        if (!arb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); fmpq_print(a); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_printf("r1 = "); arb_printn(r1, 100, 0); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printn(r2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        fmpq_clear(a);
        arb_clear(z);
        arb_clear(r1);
        arb_clear(r2);

        mag_clear(tol1);
        mag_clear(tol2);
        mag_clear(err1);
        mag_clear(err2);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong na;
        arb_t z, r1, r2;
        mag_t tol1, tol2, err1, err2;
        slong N1, N2, prec1, prec2;

        mag_init(tol1);
        mag_init(tol2);
        mag_init(err1);
        mag_init(err2);

        arb_init(z);
        arb_init(r1);
        arb_init(r2);

        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);

        na = n_randint(state, 50);

        arb_set_ui(z, 1 + n_randint(state, 100));
        arb_div_ui(z, z, 1 + n_randint(state, 100), 200);

        mag_set_ui_2exp_si(tol1, 1, -(slong) n_randint(state, 100));
        mag_set_ui_2exp_si(tol2, 1, -(slong) n_randint(state, 100));

        N1 = _arb_hypgeom_gamma_upper_singular_si_choose_N(err1, na, z, tol1);
        N2 = _arb_hypgeom_gamma_upper_singular_si_choose_N(err2, na, z, tol2);

        _arb_hypgeom_gamma_upper_singular_si_bsplit(r1, na, z, N1, prec1);
        arb_add_error_mag(r1, err1);

        _arb_hypgeom_gamma_upper_singular_si_bsplit(r2, na, z, N2, prec2);
        arb_add_error_mag(r2, err2);

        if (!arb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("na = %wd", na); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_printf("r1 = "); arb_printn(r1, 100, 0); flint_printf("\n\n");
            flint_printf("r2 = "); arb_printn(r2, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(z);
        arb_clear(r1);
        arb_clear(r2);

        mag_clear(tol1);
        mag_clear(tol2);
        mag_clear(err1);
        mag_clear(err2);
    }

    TEST_FUNCTION_END(state);
}
