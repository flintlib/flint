/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"
#include "arb_hypgeom.h"
#include "mp_real.h"

/* Gamma(1/2, y) for rational y > 0 by every method against
   sqrt(pi) erfc(sqrt(y)) in arb: small and large y (the asymptotic series from
   y of about FLINT_BITS n log 2, the continued fraction for large y),
   multi-limb numerators and denominators; the automatic choice must be
   accurate to about n limbs. */

/* the limbs of x > 0, with an extra zero limb on top if pad */
static void
_limbs(ulong * buf, slong * len, const fmpz_t x, int pad)
{
    *len = fmpz_size(x);
    fmpz_get_ui_array(buf, *len, x);
    if (pad)
        buf[(*len)++] = 0;
}

TEST_FUNCTION_START(mp_real_gamma_upper_half, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        slong n, acc;
        fmpq_t y;
        mp_real_t x;
        arb_t r, s, yy, a;
        double yd;
        ulong pd[8], qd[8];
        slong pn, qn;

        fmpq_init(y);
        mp_real_init(x);
        arb_init(r);
        arb_init(s);
        arb_init(yy);
        arb_init(a);

        n = 1 + n_randint(state, (iter % 10 == 0) ? 100 : 20);

        /* y = p/q of various sizes */
        switch (n_randint(state, 4))
        {
            case 0:
                fmpz_randtest_not_zero(fmpq_numref(y), state, 8);
                fmpz_randtest_not_zero(fmpq_denref(y), state, 8);
                break;
            case 1:
                fmpz_randtest_not_zero(fmpq_numref(y), state, 20);
                fmpz_one(fmpq_denref(y));
                break;
            case 2:
                /* large: around and beyond the asymptotic regime */
                fmpz_set_ui(fmpq_numref(y), 1 + n_randint(state, 3 * FLINT_BITS * n));
                fmpz_set_ui(fmpq_denref(y), 1 + n_randint(state, 8));
                break;
            default:
                fmpz_randtest_not_zero(fmpq_numref(y), state, 150);
                fmpz_randtest_not_zero(fmpq_denref(y), state, 150);
        }
        fmpz_abs(fmpq_numref(y), fmpq_numref(y));
        fmpz_abs(fmpq_denref(y), fmpq_denref(y));
        fmpq_canonicalise(y);
        if (n_randint(state, 4) == 0)
        {
            /* not in lowest terms */
            fmpz_mul_ui(fmpq_numref(y), fmpq_numref(y), 3);
            fmpz_mul_ui(fmpq_denref(y), fmpq_denref(y), 3);
        }
        /* the supported range is y < 2^(FLINT_BITS - 6) */
        if (fmpz_bits(fmpq_numref(y)) > fmpz_bits(fmpq_denref(y)) + FLINT_BITS - 8)
            fmpq_inv(y, y);
        yd = fmpq_get_d(y);
        _limbs(pd, &pn, fmpq_numref(y), n_randint(state, 4) == 0);
        _limbs(qd, &qn, fmpq_denref(y), n_randint(state, 4) == 0);

        arb_set_fmpq(yy, y, FLINT_BITS * n + 400);
        /* sqrt(pi) erfc(sqrt(y)) (arb_hypgeom_gamma_upper can be slow
           here) */
        arb_sqrt(s, yy, FLINT_BITS * n + 300);
        arb_hypgeom_erfc(r, s, FLINT_BITS * n + 200);
        arb_const_sqrt_pi(s, FLINT_BITS * n + 200);
        arb_mul(r, r, s, FLINT_BITS * n + 200);

        mp_real_gamma_upper_half_frac(x, pd, pn, qd, qn, n);
        mp_real_get_arb(a, x);
        acc = arb_rel_accuracy_bits(a);
        if (!arb_overlaps(a, r) || acc < FLINT_BITS * (n - 1) - 16)
        {
            flint_printf("FAIL: automatic, n = %wd, acc = %wd\n", n, acc);
            flint_printf("y = "); fmpq_print(y); flint_printf("\n");
            flint_printf("a = "); arb_printn(a, 50, 0); flint_printf("\n");
            flint_printf("r = "); arb_printn(r, 50, 0); flint_printf("\n");
            flint_abort();
        }

        if (yd < 1e5)
        {
            _mp_real_gamma_upper_half_frac_lower(x, pd, pn, qd, qn, n);
            mp_real_get_arb(a, x);
            if (!arb_overlaps(a, r))
            {
                flint_printf("FAIL: lower, n = %wd\n", n);
                flint_printf("y = "); fmpq_print(y); flint_printf("\n");
                flint_abort();
            }
        }

        if (_mp_real_gamma_upper_half_frac_asymp(x, pd, pn, qd, qn, n_randint(state, 2) ? 0 : 1 + n_randint(state, 100), n))
        {
            mp_real_get_arb(a, x);
            if (!arb_overlaps(a, r))
            {
                flint_printf("FAIL: asymp, n = %wd\n", n);
                flint_printf("y = "); fmpq_print(y); flint_printf("\n");
                flint_abort();
            }
        }

        if (yd >= 0.25 && yd < 1e12)
        {
            slong N = 1 + n_randint(state, 1000);
            if (n_randint(state, 2))
                N += (slong) (pow(FLINT_BITS * n / 2.4, 2.0) / yd);
            N = FLINT_MIN(N, WORD(100000));
            if (_mp_real_gamma_upper_half_frac_cf(x, pd, pn, qd, qn, N, n))
            {
                mp_real_get_arb(a, x);
                if (!arb_overlaps(a, r))
                {
                    flint_printf("FAIL: cf, n = %wd, N = %wd\n", n, N);
                    flint_printf("y = "); fmpq_print(y); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        /* now and then at a few thousand limbs with y = p/2^64 of about a
           third of the precision, as in the dual terms of Khinchin's
           constant: the fraction, long enough, must reach the accuracy
           (its convergents differ in size by about 2^bits(p), which the
           transformed products of the evaluation must preserve) */
        if (iter % 150 == 1)
        {
            slong n2 = 1500 + n_randint(state, 1000), N;
            double y2;
            ulong p2[2], q2[2] = { 0, 1 };
            arb_t r2;

            arb_init(r2);
            p2[0] = n_randlimb(state);
            p2[1] = (ulong) (FLINT_BITS * n2 / 3) + n_randint(state, 100);
            y2 = (double) p2[1];
            N = (slong) (2.0 * pow(FLINT_BITS * n2 / 2.4, 2.0) / y2) + 1000;

            fmpz_set_ui_array(fmpq_numref(y), p2, 2);
            fmpz_set_ui_array(fmpq_denref(y), q2, 2);
            arb_set_fmpq(yy, y, FLINT_BITS * n2 + 400);
            arb_sqrt(s, yy, FLINT_BITS * n2 + 300);
            arb_hypgeom_erfc(r2, s, FLINT_BITS * n2 + 200);
            arb_const_sqrt_pi(s, FLINT_BITS * n2 + 200);
            arb_mul(r2, r2, s, FLINT_BITS * n2 + 200);

            if (!_mp_real_gamma_upper_half_frac_cf(x, p2, 2, q2, 2, N, n2))
            {
                flint_printf("FAIL: cf not accurate, n = %wd, N = %wd, y = %g\n", n2, N, y2);
                flint_abort();
            }
            mp_real_get_arb(a, x);
            if (!arb_overlaps(a, r2))
            {
                flint_printf("FAIL: cf (large), n = %wd, N = %wd, y = %g\n", n2, N, y2);
                flint_abort();
            }
            arb_clear(r2);
        }

        fmpq_clear(y);
        mp_real_clear(x);
        arb_clear(r);
        arb_clear(s);
        arb_clear(yy);
        arb_clear(a);
    }

    TEST_FUNCTION_END(state);
}
