/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

TEST_FUNCTION_START(mp_real_atan_frac_bsplit, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        flint_set_num_threads(1 + n_randint(state, 4));
        /* both merge strategies: large q (the power table, from 64
           working limbs) and small q (the carried denominator) */
        slong n = 1 + n_randint(state, (iter % 8 == 0) ? 300 : 40);
        slong qn = 1 + n_randint(state, 3), pn;
        ulong p[3], q[3];
        int hyperbolic = n_randint(state, 2);
        fmpz_t P, Q;
        arb_t r, s;
        mp_real_t b;

        fmpz_init(P); fmpz_init(Q);
        arb_init(r); arb_init(s);
        mp_real_init(b);

        for (;;)
        {
            flint_mpn_rrandom(q, state, qn);
            if (n_randint(state, 4) == 0)
                q[0] = 2 + n_randint(state, 1000), qn = 1;
            if (q[qn - 1] == 0)
                q[qn - 1] = 1;

            if (n_randint(state, 2))
            {
                p[0] = 1;
                pn = 1;
            }
            else
            {
                pn = 1 + n_randint(state, qn);
                flint_mpn_rrandom(p, state, pn);
                if (pn == qn)
                    mpn_rshift(p, p, pn, 1 + n_randint(state, 10));
                while (pn > 0 && p[pn - 1] == 0)
                    pn--;
                if (pn == 0)
                    p[0] = 1, pn = 1;
            }

            fmpz_set_ui_array(P, p, pn);
            fmpz_set_ui_array(Q, q, qn);
            /* p/q <= 0.98 */
            fmpz_mul_ui(P, P, 100);
            fmpz_mul_ui(Q, Q, 98);
            if (fmpz_cmp(P, Q) <= 0)
                break;
        }
        fmpz_set_ui_array(P, p, pn);
        fmpz_set_ui_array(Q, q, qn);

        if (hyperbolic)
            _mp_real_atanh_frac_bsplit(b, p, pn, q, qn, n);
        else
            _mp_real_atan_frac_bsplit(b, p, pn, q, qn, n);
        mp_real_get_arb(r, b);
        arb_atan_frac_bsplit(s, P, Q, hyperbolic, FLINT_BITS * n + 64);

        if (!arb_overlaps(r, s) || arb_rel_accuracy_bits(r) < FLINT_BITS * n)
        {
            flint_printf("FAIL: n = %wd, hyperbolic = %d\np = ", n, hyperbolic);
            fmpz_print(P); flint_printf("\nq = "); fmpz_print(Q);
            flint_printf("\nr = "); arb_printd(r, 30);
            flint_printf("\ns = "); arb_printd(s, 30); flint_printf("\n");
            flint_abort();
        }

        fmpz_clear(P); fmpz_clear(Q);
        arb_clear(r); arb_clear(s);
        mp_real_clear(b);
    }

    flint_set_num_threads(1);
    TEST_FUNCTION_END(state);
}
