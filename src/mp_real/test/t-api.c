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
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"
#include "ball_helpers.h"

/* the small API entry points: setters, swap, exact-zero predicate,
   mpn import, fixed-point export with its rigorous bound, and the
   explicit error-widening helpers */
static void
test_api(flint_rand_t state, slong iters)
{
    slong iter;

    for (iter = 0; iter < iters; iter++)
    {
        mp_real_t x, y;
        arb_t ax, ay;
        arf_t t;

        mp_real_init(x); mp_real_init(y);
        arb_init(ax); arb_init(ay);
        arf_init(t);

        /* set_ui / set_si / is_zero_exact */
        {
            ulong c = n_randtest(state);
            slong d = (slong) n_randtest(state);

            mp_real_set_ui(x, c);
            mp_real_get_arb(ax, x);
            if (!arb_equal_si(ax, (slong) 0) && c == 0)
                { flint_printf("set_ui zero\n"); flint_abort(); }
            arb_set_ui(ay, c);
            if (!arb_equal(ax, ay))
                { flint_printf("set_ui: c = %wu\n", c); flint_abort(); }
            if (!mp_real_is_zero(x) != !(c == 0))
                { flint_printf("is_zero_exact ui\n"); flint_abort(); }

            mp_real_set_si(x, d);
            mp_real_get_arb(ax, x);
            arb_set_si(ay, d);
            if (!arb_equal(ax, ay))
                { flint_printf("set_si: d = %wd\n", d); flint_abort(); }
            _mp_real_add_error_ulps(x, 2.0);
            if (mp_real_is_zero(x))
                { flint_printf("is_zero_exact with radius\n"); flint_abort(); }
        }

        /* swap */
        {
            mp_real_randtest(x, state, 6, 1);
            mp_real_randtest(y, state, 6, 1);
            mp_real_get_arb(ax, x);
            mp_real_get_arb(ay, y);
            mp_real_swap(x, y);
            {
                arb_t bx, by;
                arb_init(bx); arb_init(by);
                mp_real_get_arb(bx, x);
                mp_real_get_arb(by, y);
                if (!arb_equal(bx, ay) || !arb_equal(by, ax))
                    { flint_printf("swap\n"); flint_abort(); }
                arb_clear(bx); arb_clear(by);
            }
        }

        /* set_mpn_2exp: exact import */
        {
            ulong p[4];
            slong len = 1 + (slong) n_randint(state, 4);
            slong ebits = (slong) n_randint(state, 400) - 200;
            slong i;

            for (i = 0; i < len; i++)
                p[i] = n_randtest(state);
            _mp_real_set_mpn_2exp(x, p, len, ebits);
            mp_real_get_arb(ax, x);
            {
                fmpz_t f;
                fmpz_init(f);
                fmpz_set_ui_array(f, p, len);
                arb_set_fmpz(ay, f);
                arb_mul_2exp_si(ay, ay, ebits);
                if (!arb_equal(ax, ay))
                    { flint_printf("set_mpn_2exp: len = %wd, "
                        "ebits = %wd\n", len, ebits); flint_abort(); }
                fmpz_clear(f);
            }
        }

        /* add_error / add_error_ulps: the original point value must
           remain contained after widening */
        {
            mp_real_randtest(x, state, 6, 1);
            mp_real_random_point(t, x, state);
            _mp_real_add_error_ulps(x, (double) (1 + n_randint(state, 100)));
            check_contains(x, t, "add_error_ulps", iter);
            _mp_real_add_error_ulps_at(x, (double) (1 + n_randint(state, 100)),
                x->exp - x->size - (slong) n_randint(state, 3));
            check_contains(x, t, "add_error", iter);
        }

        /* get_fixed: build a ball inside [0, 1), export at wn limbs,
           and check the true point against the returned ulp bound */
        {
            ulong p[3], f[6];
            slong len = 1 + (slong) n_randint(state, 3);
            slong wn = 1 + (slong) n_randint(state, 5);
            slong i;
            ulong bound;

            for (i = 0; i < len; i++)
                p[i] = n_randtest(state);
            p[len - 1] |= UWORD(1) << (FLINT_BITS - 1);
            /* value in [1/2, 1) * 2^-shift */
            _mp_real_set_mpn_2exp(x, p, len,
                -FLINT_BITS * len - (slong) n_randint(state, 40));
            if (n_randint(state, 2))
                _mp_real_add_error_ulps(x,
                    (double) (1 + n_randint(state, 50)));
            mp_real_random_point(t, x, state);
            if (arf_sgn(t) < 0)
                arf_zero(t);

            _mp_real_get_fixed(f, &bound, x, wn);

            {
                arf_t u, w;
                arf_init(u); arf_init(w);
                {
                    fmpz_t g;
                    fmpz_init(g);
                    fmpz_set_ui_array(g, f, wn);
                    arf_set_fmpz(u, g);
                    arf_mul_2exp_si(u, u, -FLINT_BITS * wn);
                    fmpz_clear(g);
                }
                arf_sub(u, t, u, ARF_PREC_EXACT, ARF_RND_NEAR);
                arf_abs(u, u);
                arf_set_ui(w, bound);
                arf_mul_2exp_si(w, w, -FLINT_BITS * wn);
                if (bound != UWORD_MAX && arf_cmp(u, w) > 0)
                    { flint_printf("get_fixed bound: wn = %wd, "
                        "bound = %wu\n", wn, bound); flint_abort(); }
                arf_clear(u); arf_clear(w);
            }
        }

        /* get_fixed_floor: on success the output must equal the
           EXACT floor of every point of the ball -- checked against
           the floor of a random true point via arf -- and a radius
           spanning a grid line must be rejected */
        {
            ulong p[3], f[6];
            slong len = 1 + (slong) n_randint(state, 3);
            slong wn = 1 + (slong) n_randint(state, 5);
            slong i;

            for (i = 0; i < len; i++)
                p[i] = n_randtest(state);
            p[len - 1] |= UWORD(1) << (FLINT_BITS - 1);
            _mp_real_set_mpn_2exp(x, p, len,
                -FLINT_BITS * len - (slong) n_randint(state, 3));
            if (n_randint(state, 2))
                _mp_real_add_error_ulps(x,
                    (double) (1 + n_randint(state, 1000)));

            if (_mp_real_get_fixed_floor(f, wn, x))
            {
                fmpz_t g, h;
                fmpz_init(g); fmpz_init(h);
                mp_real_random_point(t, x, state);
                if (arf_sgn(t) < 0)
                    arf_zero(t);
                arf_mul_2exp_si(t, t, FLINT_BITS * wn);
                arf_get_fmpz(g, t, ARF_RND_FLOOR);
                fmpz_set_ui_array(h, f, wn);
                if (!fmpz_equal(g, h))
                    { flint_printf("FAIL: get_fixed_floor value "
                        "(iter %wd)\n", iter); flint_abort(); }
                fmpz_clear(g); fmpz_clear(h);
            }

            /* a ball straddling a grid line must be refused: center
               the value ON a grid multiple with a nonzero radius */
            _mp_real_set_mpn_2exp(x, p, 1, -FLINT_BITS);
            _mp_real_add_error_ulps(x, 1.0);
            if (p[0] != 0 && _mp_real_get_fixed_floor(f, 1, x))
                { flint_printf("FAIL: get_fixed_floor accepted a "
                    "grid-straddling ball (iter %wd)\n", iter);
                  flint_abort(); }
        }

        mp_real_clear(x); mp_real_clear(y);
        arb_clear(ax); arb_clear(ay);
        arf_clear(t);
    }
}

TEST_FUNCTION_START(mp_real_api, state)
{
    test_api(state, 1000 + 1000 * flint_test_multiplier());

    TEST_FUNCTION_END(state);
}
