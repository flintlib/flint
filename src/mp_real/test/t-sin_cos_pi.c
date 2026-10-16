/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpq.h"
#include "arb.h"
#include "mp_real.h"

/* _mp_real_sin_cos_pi_ui_div_ui (every algorithm and order, the closed
   forms for q = 3, 4, 5, 6, 8, 10, 12) and
   mp_real_sin_cos_pi_ui_div_ui against arb_sin_cos_pi_fmpq: the
   fixed-point outputs lie within their bound, which respects the
   documented maximum (a few ulps for the Chebyshev iteration); the
   balls contain the values (overlap the reference), have a relative accuracy of about n limbs,
   and are exact for 0, +-1, +-1/2; either output may be NULL.  Also
   mp_real_tan_pi_ui_div_ui. */

/* a denominator: small, moderate, a word, near a power of two */
static ulong
scp_rand_q(flint_rand_t state)
{
    ulong q;

    switch (n_randint(state, 6))
    {
        case 0:
            q = 1 + n_randint(state, 16);
            break;
        case 1:
            q = 1 + n_randint(state, 1000);
            break;
        case 2:
            q = 1 + n_randint(state, 1000000);
            break;
        case 3:
            q = n_randtest_not_zero(state);
            break;
        case 4:
            q = UWORD(1) << n_randint(state, FLINT_BITS - 1);
            q += n_randint(state, 3);
            break;
        default:
            q = n_randprime(state, 2 + n_randint(state, FLINT_BITS - 3), 0);
            break;
    }

    return FLINT_MAX(q, 1);
}

/* the fixed-point value (y, n + 1) +- err ulps overlaps ref */
static int
scp_fixed_ok(nn_srcptr y, ulong err, slong n, const arb_t ref)
{
    arb_t t;
    fmpz_t z;
    mag_t e;
    int ok;

    arb_init(t);
    fmpz_init(z);
    mag_init(e);

    fmpz_set_ui_array(z, y, n + 1);
    arb_set_fmpz(t, z);
    arb_mul_2exp_si(t, t, -FLINT_BITS * n);
    mag_set_ui(e, err);
    mag_mul_2exp_si(e, e, -FLINT_BITS * n);
    mag_add(arb_radref(t), arb_radref(t), e);
    ok = arb_overlaps(t, ref);

    arb_clear(t);
    fmpz_clear(z);
    mag_clear(e);
    return ok;
}

TEST_FUNCTION_START(mp_real_sin_cos_pi, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        ulong p, q, err, maxerr;
        slong n;
        int alg, r, which;
        nn_ptr ys, yc;
        arb_t sref, cref;
        fmpq_t pq;

        q = scp_rand_q(state);
        p = n_randint(state, q / 2 + 1);
        if (n_randint(state, 4) == 0)
            p = n_randint(state, 2) ? FLINT_MIN(1, q / 2) : q / 2;

        if (iter % 50 == 0)
            n = 1 + n_randint(state, 1500);
        else if (iter % 5 == 0)
            n = 1 + n_randint(state, 300);
        else
            n = 1 + n_randint(state, 40);

        alg = n_randint(state, 4);
        r = n_randint(state, 2) ? 0 : 2 + n_randint(state, 15);
        which = n_randint(state, 3);    /* 0 both, 1 sin only, 2 cos only */

        arb_init(sref);
        arb_init(cref);
        fmpq_init(pq);

        fmpq_set_ui(pq, p, q);
        arb_sin_cos_pi_fmpq(sref, cref, pq, FLINT_BITS * n + 100);

        ys = flint_malloc((n + 1) * sizeof(ulong));
        yc = flint_malloc((n + 1) * sizeof(ulong));

        _mp_real_sin_cos_pi_ui_div_ui_tune(which != 2 ? ys : NULL,
            which != 1 ? yc : NULL, &err, p, q, n, alg, r);

        /* the kernel's bound, or a few ulps for the iteration and the
           closed forms */
        maxerr = 6 * 768 + 130;
        {
            ulong qr = q / n_gcd(p, q);
            if ((alg == 2 && qr <= (UWORD(1) << (FLINT_BITS - 9)))
                || qr <= 4 || qr == 6
                || (alg == 3 && (qr == 5 || qr == 8 || qr == 10 || qr == 12)))
                maxerr = 4;
        }

        if (err > maxerr
            || (which != 2 && !scp_fixed_ok(ys, err, n, sref))
            || (which != 1 && !scp_fixed_ok(yc, err, n, cref)))
        {
            flint_printf("FAIL: fixed, p = %wu, q = %wu, n = %wd, alg = %d, r = %d, which = %d, err = %wu\n",
                p, q, n, alg, r, which, err);
            flint_printf("sin = "); arb_printd(sref, 30); flint_printf("\n");
            flint_printf("cos = "); arb_printd(cref, 30); flint_printf("\n");
            flint_abort();
        }

        flint_free(ys);
        flint_free(yc);

        /* the balls, for any p */
        {
            mp_real_t s, c;
            arb_t ts, tc;
            slong m = 1 + n_randint(state, (iter % 10 == 0) ? 200 : 20);

            mp_real_init(s);
            mp_real_init(c);
            arb_init(ts);
            arb_init(tc);

            p = n_randint(state, 2) ? n_randtest(state) : n_randint(state, 4 * q + 1);
            fmpq_set_ui(pq, p, q);
            arb_sin_cos_pi_fmpq(sref, cref, pq, FLINT_BITS * (m + 4) + 100);

            mp_real_sin_cos_pi_ui_div_ui(which != 2 ? s : NULL,
                which != 1 ? c : NULL, p, q, m);
            mp_real_get_arb(ts, s);
            mp_real_get_arb(tc, c);

            if ((which != 2 && !arb_overlaps(ts, sref))
                || (which != 1 && !arb_overlaps(tc, cref)))
            {
                flint_printf("FAIL: ball containment, p = %wu, q = %wu, n = %wd\n", p, q, m);
                flint_printf("sin = "); arb_printd(ts, 30); flint_printf("\n");
                flint_printf("      "); arb_printd(sref, 30); flint_printf("\n");
                flint_printf("cos = "); arb_printd(tc, 30); flint_printf("\n");
                flint_printf("      "); arb_printd(cref, 30); flint_printf("\n");
                flint_abort();
            }

            if ((which != 2 && (arb_is_exact(sref) ? !arb_is_exact(ts)
                    : arb_rel_accuracy_bits(ts) < FLINT_BITS * m - 2))
                || (which != 1 && (arb_is_exact(cref) ? !arb_is_exact(tc)
                    : arb_rel_accuracy_bits(tc) < FLINT_BITS * m - 2)))
            {
                flint_printf("FAIL: ball accuracy, p = %wu, q = %wu, n = %wd\n", p, q, m);
                flint_printf("sin = "); arb_printd(ts, 30); flint_printf("\n");
                flint_printf("cos = "); arb_printd(tc, 30); flint_printf("\n");
                flint_abort();
            }

            mp_real_clear(s);
            mp_real_clear(c);
            arb_clear(ts);
            arb_clear(tc);
        }

        arb_clear(sref);
        arb_clear(cref);
        fmpq_clear(pq);
    }

    /* mp_real_tan_pi_ui_div_ui: status 0 exactly at the poles, the value
       otherwise, exact for 0, +-1, with a relative accuracy of about n
       limbs */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        ulong p, q;
        slong n;
        int ok;
        mp_real_t y;
        arb_t s, c, t, a;
        fmpq_t pq;

        q = scp_rand_q(state);
        p = n_randint(state, 2) ? n_randtest(state) : n_randint(state, 4 * q + 1);
        n = 1 + n_randint(state, (iter % 20 == 0) ? 1500 : 40);

        mp_real_init(y);
        arb_init(s);
        arb_init(c);
        arb_init(t);
        arb_init(a);
        fmpq_init(pq);

        ok = mp_real_tan_pi_ui_div_ui(y, p, q, n);
        fmpq_set_ui(pq, p, q);
        arb_sin_cos_pi_fmpq(s, c, pq, FLINT_BITS * n + 200);
        mp_real_get_arb(a, y);

        if (arb_is_zero(c))
        {
            if (ok || !mp_real_is_zero(y))
            {
                flint_printf("FAIL: tan pole, p = %wu, q = %wu\n", p, q);
                flint_abort();
            }
        }
        else
        {
            arb_div(t, s, c, FLINT_BITS * n + 200);
            if (!ok || !arb_overlaps(a, t) || (arb_is_exact(t) ? !arb_equal(a, t)
                    : arb_rel_accuracy_bits(a) < FLINT_BITS * n - 2))
            {
                flint_printf("FAIL: tan, p = %wu, q = %wu, n = %wd, ok = %d\n", p, q, n, ok);
                flint_printf("y = "); arb_printd(a, 30); flint_printf("\n");
                flint_printf("    "); arb_printd(t, 30); flint_printf("\n");
                flint_abort();
            }
        }

        mp_real_clear(y);
        arb_clear(s);
        arb_clear(c);
        arb_clear(t);
        arb_clear(a);
        fmpq_clear(pq);
    }

    TEST_FUNCTION_END(state);
}
