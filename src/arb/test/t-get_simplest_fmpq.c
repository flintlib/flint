/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "arb.h"
#include "fmpq.h"

TEST_FUNCTION_START(arb_get_simplest_fmpq, state)
{
    /* Test 1: exact zero ball -> 0/1 */
    {
        arb_t b;
        fmpq_t res, expected;
        arb_init(b);
        fmpq_init(res);
        fmpq_init(expected);

        arb_zero(b);

        if (!arb_get_simplest_fmpq(res, b) || !fmpq_equal(res, expected))
        {
            flint_printf("FAIL: exact zero\n");
            flint_printf("res = "); fmpq_print(res); flint_printf("\n");
            flint_abort();
        }

        arb_clear(b);
        fmpq_clear(res);
        fmpq_clear(expected);
    }

    /* Test 2: ball straddling zero with radius 0.5 -> 0/1 */
    {
        arb_t b;
        fmpq_t res, expected;
        arb_init(b);
        fmpq_init(res);
        fmpq_init(expected);

        /* b = 0.1 +/- 0.5 contains 0 */
        arb_set_d(b, 0.1);
        mag_set_d(arb_radref(b), 0.5);

        if (!arb_get_simplest_fmpq(res, b) || !fmpq_equal(res, expected))
        {
            flint_printf("FAIL: zero-crossing ball\n");
            flint_printf("res = "); fmpq_print(res); flint_printf("\n");
            flint_abort();
        }

        arb_clear(b);
        fmpq_clear(res);
        fmpq_clear(expected);
    }

    /* Test 3: exact integer 7 -> 7/1 */
    {
        arb_t b;
        fmpq_t res;
        arb_init(b);
        fmpq_init(res);

        arb_set_si(b, 7);

        if (!arb_get_simplest_fmpq(res, b)
            || fmpz_cmp_si(fmpq_numref(res), 7) != 0
            || fmpz_cmp_si(fmpq_denref(res), 1) != 0)
        {
            flint_printf("FAIL: exact integer 7\n");
            flint_printf("res = "); fmpq_print(res); flint_printf("\n");
            flint_abort();
        }

        arb_clear(b);
        fmpq_clear(res);
    }

    /* Test 4: ball [2.3, 3.7] contains integer 3 -> 3/1 */
    {
        arb_t b;
        fmpq_t res;
        arf_t lo_arf, hi_arf;
        arb_init(b);
        fmpq_init(res);
        arf_init(lo_arf);
        arf_init(hi_arf);

        arf_set_d(lo_arf, 2.3);
        arf_set_d(hi_arf, 3.7);
        arb_set_interval_arf(b, lo_arf, hi_arf, 53);

        if (!arb_get_simplest_fmpq(res, b)
            || fmpz_cmp_si(fmpq_numref(res), 3) != 0
            || fmpz_cmp_si(fmpq_denref(res), 1) != 0)
        {
            flint_printf("FAIL: ball [2.3, 3.7] -> 3\n");
            flint_printf("res = "); fmpq_print(res); flint_printf("\n");
            flint_abort();
        }

        arb_clear(b);
        fmpq_clear(res);
        arf_clear(lo_arf);
        arf_clear(hi_arf);
    }

    /* Test 5: invalid ball (NaN) -> returns 0 */
    {
        arb_t b;
        fmpq_t res;
        arb_init(b);
        fmpq_init(res);

        arb_indeterminate(b);

        if (arb_get_simplest_fmpq(res, b) != 0)
        {
            flint_printf("FAIL: NaN ball should return 0\n");
            flint_abort();
        }

        arb_clear(b);
        fmpq_clear(res);
    }

    /* Test 6: exact 3/2 -> 3/2 */
    {
        arb_t b;
        fmpq_t res;
        arb_init(b);
        fmpq_init(res);

        arb_set_str(b, "1.5 +/- 0", 53);

        if (!arb_get_simplest_fmpq(res, b)
            || fmpz_cmp_si(fmpq_numref(res), 3) != 0
            || fmpz_cmp_si(fmpq_denref(res), 2) != 0)
        {
            flint_printf("FAIL: exact 3/2\n");
            flint_printf("res = "); fmpq_print(res); flint_printf("\n");
            flint_abort();
        }

        arb_clear(b);
        fmpq_clear(res);
    }

    /* Test 7: tight ball around 355/113 at prec 53 -> 355/113 */
    {
        arb_t b, target;
        fmpq_t res;
        fmpz_t num, den;
        arb_init(b);
        arb_init(target);
        fmpq_init(res);
        fmpz_init_set_si(num, 355);
        fmpz_init_set_si(den, 113);

        arb_set_fmpz(target, num);
        arb_div_fmpz(b, target, den, 53);

        if (!arb_get_simplest_fmpq(res, b)
            || fmpz_cmp(fmpq_numref(res), num) != 0
            || fmpz_cmp(fmpq_denref(res), den) != 0)
        {
            flint_printf("FAIL: 355/113 at prec 53\n");
            flint_printf("res = "); fmpq_print(res); flint_printf("\n");
            flint_abort();
        }

        arb_clear(b);
        arb_clear(target);
        fmpq_clear(res);
        fmpz_clear(num);
        fmpz_clear(den);
    }

    /* Test 8: negative ball around -355/113 -> -355/113 */
    {
        arb_t b, target;
        fmpq_t res;
        fmpz_t num, den;
        arb_init(b);
        arb_init(target);
        fmpq_init(res);
        fmpz_init_set_si(num, -355);
        fmpz_init_set_si(den, 113);

        arb_set_fmpz(target, num);
        arb_div_fmpz(b, target, den, 53);

        if (!arb_get_simplest_fmpq(res, b)
            || fmpz_cmp(fmpq_numref(res), num) != 0
            || fmpz_cmp(fmpq_denref(res), den) != 0)
        {
            flint_printf("FAIL: -355/113 at prec 53\n");
            flint_printf("res = "); fmpq_print(res); flint_printf("\n");
            flint_abort();
        }

        arb_clear(b);
        arb_clear(target);
        fmpq_clear(res);
        fmpz_clear(num);
        fmpz_clear(den);
    }

    /* Test 9: exact 355/113 at prec 50000 -- stack-safety regression. */
    {
        arb_t b, target;
        fmpq_t res;
        fmpz_t num, den;
        arb_init(b);
        arb_init(target);
        fmpq_init(res);
        fmpz_init_set_si(num, 355);
        fmpz_init_set_si(den, 113);

        arb_set_fmpz(target, num);
        arb_div_fmpz(b, target, den, 50000);

        if (!arb_get_simplest_fmpq(res, b)
            || fmpz_cmp(fmpq_numref(res), num) != 0
            || fmpz_cmp(fmpq_denref(res), den) != 0)
        {
            flint_printf("FAIL: 355/113 at prec 50000\n");
            flint_printf("res = "); fmpq_print(res); flint_printf("\n");
            flint_abort();
        }

        arb_clear(b);
        arb_clear(target);
        fmpq_clear(res);
        fmpz_clear(num);
        fmpz_clear(den);
    }

    /* Test 10: ball around pi at prec 50000.
     * Asserts: function returns success and runs to completion without
     * crashing. We do not constrain the exact answer (it depends on
     * the precision-50000 enclosure of pi), only its denominator size. */
    {
        arb_t b;
        fmpq_t res;
        flint_bitcnt_t den_bits;
        arb_init(b);
        fmpq_init(res);

        arb_const_pi(b, 50000);

        if (!arb_get_simplest_fmpq(res, b))
        {
            flint_printf("FAIL: pi at prec 50000 returned 0\n");
            flint_abort();
        }

        /* Sanity: denominator should be < 2^25000 (CF length is roughly
         * half the precision for typical irrationals). */
        den_bits = fmpz_bits(fmpq_denref(res));
        if (den_bits > 25000)
        {
            flint_printf("FAIL: pi at prec 50000: denominator too large\n");
            flint_printf("den_bits = %wu\n", (ulong) den_bits);
            flint_abort();
        }

        arb_clear(b);
        fmpq_clear(res);
    }

    /* Test 11: ball with infinite radius (represents the whole real
     * line). Not finite, so arb_get_simplest_fmpq must return 0. */
    {
        arb_t b;
        fmpq_t res;
        arb_init(b);
        fmpq_init(res);

        arb_zero_pm_inf(b);

        if (arb_get_simplest_fmpq(res, b) != 0)
        {
            flint_printf("FAIL: infinite radius ball should return 0\n");
            flint_abort();
        }

        arb_clear(b);
        fmpq_clear(res);
    }

    /* Test 12: random fuzz. For random balls, the returned fmpq must
     * actually lie in the ball. */
    {
        slong iter;

        for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
        {
            arb_t b;
            fmpq_t res;
            int ok;

            arb_init(b);
            fmpq_init(res);

            arb_randtest(b, state, 200, 5);

            ok = arb_get_simplest_fmpq(res, b);

            if (ok && !arb_contains_fmpq(b, res))
            {
                flint_printf("FAIL: returned fmpq not contained in ball\n");
                flint_printf("b = "); arb_printd(b, 30); flint_printf("\n");
                flint_printf("res = "); fmpq_print(res); flint_printf("\n");
                flint_abort();
            }

            arb_clear(b);
            fmpq_clear(res);
        }
    }

    /* Test 13: theoretical-bound coverage.
     *
     * For p/q in lowest terms (q >= 1), any other rational r/s in
     * lowest terms satisfies |r/s - p/q| = |r*q - p*s|/(s*q) >= 1/(s*q)
     * because |r*q - p*s| >= 1 when the fractions differ. So a closed
     * ball around p/q of radius < 1/((q-1)*q) (for q >= 2) contains no
     * rational with denominator < q, and p/q is therefore the simplest
     * rational in the ball. The same bound with s >= 1 covers q = 1:
     * |r/s - p| >= 1/s, so any sub-unit radius around an integer p
     * leaves p/1 as the unique simplest rational (returned either by
     * the contains-zero fast path for p = 0 or by fmpq_simplest_between
     * otherwise).
     *
     * For the random p/q below, the canonical q is at most 1000, so
     * the theoretical bound demands radius < ~1e-6. arb_set_fmpq at
     * prec 100 yields radius below 2^-89 ~ 1.6e-27, comfortably inside
     * the bound. */
    {
        slong iter;

        for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
        {
            arb_t b;
            fmpq_t target, res;
            slong p_raw;
            ulong q_raw;

            arb_init(b);
            fmpq_init(target);
            fmpq_init(res);

            p_raw = (slong) n_randint(state, 2001) - WORD(1000);
            q_raw = 1 + n_randint(state, 1000);

            fmpq_set_si(target, p_raw, q_raw);
            fmpq_canonicalise(target);
            arb_set_fmpq(b, target, 100);

            if (!arb_get_simplest_fmpq(res, b)
                || !fmpq_equal(res, target))
            {
                flint_printf("FAIL: theoretical-bound test\n");
                flint_printf("target = ");
                fmpq_print(target); flint_printf("\n");
                flint_printf("b = ");
                arb_printd(b, 30); flint_printf("\n");
                flint_printf("res = ");
                fmpq_print(res); flint_printf("\n");
                flint_abort();
            }

            arb_clear(b);
            fmpq_clear(target);
            fmpq_clear(res);
        }
    }

    TEST_FUNCTION_END(state);
}
