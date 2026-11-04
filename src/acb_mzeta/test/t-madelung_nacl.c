/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mzeta.h"

TEST_FUNCTION_START(acb_mzeta_madelung_nacl, state)
{
    slong iter;
    arb_t ref;

    arb_init(ref);
    arb_set_str(ref, "1.747564594633182190636212035544397403485161436624741758152825350765040623532761179890758362694607889930832581538753710593282029944183828013036933002156599363282376607172297568659238037167203810410603 +/- 1e-197", 700);

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        slong prec = 2 + n_randint(state, (iter < 5) ? 6000 : 600);
        slong threads = 1 + n_randint(state, 3);
        arb_t M;
        acb_t z, nu;
        arb_mat_t A;
        arb_ptr x, y;

        arb_init(M);
        acb_init(z);
        acb_init(nu);
        arb_mat_init(A, 3, 3);
        x = _arb_vec_init(3);
        y = _arb_vec_init(3);

        flint_set_num_threads(threads);
        acb_mzeta_madelung_nacl(M, prec);
        flint_set_num_threads(1);

        if (!arb_overlaps(M, ref) || arb_rel_accuracy_bits(M) < prec - 4)
        {
            flint_printf("FAIL: madelung_nacl\n");
            flint_printf("prec = %wd, threads = %wd\n", prec, threads);
            arb_printd(M, 50); flint_printf("\n");
            flint_abort();
        }

        /* the same through the Epstein interface, x = (1, 0, 2) */
        arb_mat_one(A);
        acb_one(nu);
        arb_set_ui(x, 1);
        arb_set_ui(x + 2, 2);
        arb_set_d(y, 0.5);
        arb_set_d(y + 1, -0.5);
        arb_set_d(y + 2, 1.5);
        if (!acb_mzeta_epstein_shell(z, nu, A, x, y, 0, NULL, prec)
            || !arb_overlaps(acb_realref(z), M) || !arb_is_zero(acb_imagref(z)))
        {
            flint_printf("FAIL: epstein_shell\n");
            flint_printf("prec = %wd\n", prec);
            acb_printd(z, 50); flint_printf("\n");
            flint_abort();
        }

        arb_clear(M);
        acb_clear(z);
        acb_clear(nu);
        arb_mat_clear(A);
        _arb_vec_clear(x, 3);
        _arb_vec_clear(y, 3);
    }

    arb_clear(ref);
    TEST_FUNCTION_END(state);
}
