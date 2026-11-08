/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

/* arb_tan_pi_fmpq, and with it arb_sin_cos_pi_fmpq, arb_sin_pi_fmpq and
   arb_cos_pi_fmpq on small and word-size denominators (the mp_real
   paths, the square roots, the octant reduction): overlap with the
   functions of the exact rational at a higher precision, an
   indeterminate tangent exactly at the poles, and a relative accuracy
   within a few bits of prec */
TEST_FUNCTION_START(arb_tan_pi_fmpq, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        arb_t s, c, t, s2, c2, rs, rc, rt, u;
        fmpq_t x;
        fmpz_t p, q;
        slong prec, wp;

        arb_init(s);
        arb_init(c);
        arb_init(t);
        arb_init(s2);
        arb_init(c2);
        arb_init(rs);
        arb_init(rc);
        arb_init(rt);
        arb_init(u);
        fmpq_init(x);
        fmpz_init(p);
        fmpz_init(q);

        switch (n_randint(state, 5))
        {
            case 0:
                fmpz_set_ui(q, 1 + n_randint(state, 13));
                break;
            case 1:
                fmpz_set_ui(q, 1 + n_randint(state, 1000));
                break;
            case 2:
                /* 2^k q0, q0 = 1, 3, 5 */
                fmpz_set_ui(q, (n_randint(state, 3) == 0) ? 1 : (n_randint(state, 2) ? 3 : 5));
                fmpz_mul_2exp(q, q, n_randint(state, 20));
                break;
            case 3:
                fmpz_set_ui(q, n_randtest_not_zero(state));
                break;
            default:
                fmpz_randtest_not_zero(q, state, 100);
                fmpz_abs(q, q);
                break;
        }
        fmpz_randtest(p, state, n_randint(state, 2) ? 10 : 80);
        fmpz_set(fmpq_numref(x), p);
        fmpz_set(fmpq_denref(x), q);
        fmpq_canonicalise(x);

        /* (the square roots through arb_sqrt run up to 12288 bits) */
        prec = 2 + n_randint(state, (iter % 20 == 0) ? 12000 : (iter % 4 == 0) ? 12500 : 800);
        wp = 3 * prec + 100 + fmpz_bits(fmpq_denref(x));

        arb_tan_pi_fmpq(t, x, prec);
        arb_sin_cos_pi_fmpq(s, c, x, prec);
        arb_sin_pi_fmpq(s2, x, prec);
        arb_cos_pi_fmpq(c2, x, prec);

        /* the reference from the exact rational */
        arb_set_fmpq(u, x, wp);
        arb_sin_cos_pi(rs, rc, u, wp);
        arb_div(rt, rs, rc, wp);

        if (!arb_overlaps(s, rs) || !arb_overlaps(c, rc)
            || !arb_overlaps(s2, rs) || !arb_overlaps(c2, rc)
            || (arb_is_finite(t) && !arb_overlaps(t, rt))
            || (!arb_is_finite(t) && !arb_is_zero(c)))
        {
            flint_printf("FAIL: overlap, prec = %wd\n", prec);
            flint_printf("x = "); fmpq_print(x); flint_printf("\n");
            flint_printf("s = "); arb_printd(s, 30); flint_printf("  "); arb_printd(rs, 30); flint_printf("\n");
            flint_printf("c = "); arb_printd(c, 30); flint_printf("  "); arb_printd(rc, 30); flint_printf("\n");
            flint_printf("t = "); arb_printd(t, 30); flint_printf("  "); arb_printd(rt, 30); flint_printf("\n");
            flint_abort();
        }

        if ((!arb_is_zero(s) && arb_rel_accuracy_bits(s) < prec - 5)
            || (!arb_is_zero(c) && arb_rel_accuracy_bits(c) < prec - 5)
            || (arb_is_finite(t) && !arb_is_zero(t) && arb_rel_accuracy_bits(t) < prec - 6))
        {
            flint_printf("FAIL: accuracy, prec = %wd\n", prec);
            flint_printf("x = "); fmpq_print(x); flint_printf("\n");
            flint_printf("s = "); arb_printd(s, 30); flint_printf("\n");
            flint_printf("c = "); arb_printd(c, 30); flint_printf("\n");
            flint_printf("t = "); arb_printd(t, 30); flint_printf("\n");
            flint_abort();
        }

        arb_clear(s);
        arb_clear(c);
        arb_clear(t);
        arb_clear(s2);
        arb_clear(c2);
        arb_clear(rs);
        arb_clear(rc);
        arb_clear(rt);
        arb_clear(u);
        fmpq_clear(x);
        fmpz_clear(p);
        fmpz_clear(q);
    }

    TEST_FUNCTION_END(state);
}
