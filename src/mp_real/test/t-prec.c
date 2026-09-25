/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

/* mp_real_prec_bits and mp_real_prec_limbs: at the precision they
   return, the basic operations give a relative radius of at most
   4 eps on exact operands and 16 eps on operands of relative radius at
   most eps (eps = 2^-b, B^-p), at every size (measured worst cases up
   to 20000 limbs: about 2 and 12). */

/* (relative radius of x) / eps, eps = 2^-b; x nonzero */
static double
prec_rel(const mp_real_t x, slong b)
{
    arb_t a;
    arf_t r, m;
    double d;

    arb_init(a);
    arf_init(r);
    arf_init(m);
    mp_real_get_arb(a, x);
    arf_set_mag(r, arb_radref(a));
    arf_abs(m, arb_midref(a));
    arf_div(r, r, m, 60, ARF_RND_UP);
    arf_mul_2exp_si(r, r, b);
    d = arf_get_d(r, ARF_RND_UP);
    arb_clear(a);
    arf_clear(r);
    arf_clear(m);
    return d;
}

/* a random nonzero operand; if inexact, of relative radius <= 2^-b */
static void
prec_rand(mp_real_t x, flint_rand_t state, slong b, int inexact)
{
    for (;;)
    {
        slong L = (inexact ? (b + FLINT_BITS - 1) / FLINT_BITS : 0) + 1
            + n_randint(state, (b + FLINT_BITS - 1) / FLINT_BITS + 4);
        nn_ptr p = flint_malloc(L * sizeof(ulong));

        flint_mpn_rrandom(p, state, L);
        if (p[L - 1] == 0)
            p[L - 1] = 1;
        _mp_real_set_mpn_2exp(x, p, L, (slong) n_randint(state, 400) - 200 - FLINT_BITS * L);
        flint_free(p);
        if (n_randint(state, 2))
            mp_real_neg(x, x);
        if (!inexact)
            return;
        mp_real_add_rel_error_2exp_si(x, -b - 1);
        if (prec_rel(x, b) <= 1.0)
            return;
    }
}

TEST_FUNCTION_START(mp_real_prec, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        int inexact = n_randint(state, 2), op, use_limbs = n_randint(state, 2);
        slong b, n;
        ulong c = n_randtest_not_zero(state);
        mp_real_t x, y, ax, ay, z;

        /* sqrt and division need radii below 2^-30 */
        if (iter % 200 == 0)
        {
            /* the windowed middle products (from 64 kept limbs) and the
               Newton division and square roots (from about 900 and
               1900 limbs) */
            slong p = 60 + n_randint(state, 2500);
            b = FLINT_BITS * p - n_randint(state, FLINT_BITS);
            n = use_limbs ? mp_real_prec_limbs(p) : mp_real_prec_bits(b);
            if (use_limbs)
                b = FLINT_BITS * p;
        }
        else if (use_limbs)
        {
            slong p = 1 + n_randint(state, (iter % 10 == 0) ? 60 : 10);
            b = FLINT_BITS * p;
            n = mp_real_prec_limbs(p);
        }
        else
        {
            b = (inexact ? 40 : 1) + n_randint(state, (iter % 10 == 0) ? 4000 : 600);
            n = mp_real_prec_bits(b);
        }

        mp_real_init(x);
        mp_real_init(y);
        mp_real_init(ax);
        mp_real_init(ay);
        mp_real_init(z);

        prec_rand(x, state, b, inexact);
        prec_rand(y, state, b, inexact);
        mp_real_set(ax, x);
        if (ax->negative)
            mp_real_neg(ax, ax);
        mp_real_set(ay, y);
        if (ay->negative)
            mp_real_neg(ay, ay);

        for (op = 0; op < 10; op++)
        {
            double d;

            switch (op)
            {
                case 0: mp_real_add(z, ax, ay, n); break;
                case 1:
                    /* |x'| >= 2 |y|: no cancellation */
                    mp_real_mul_2exp_si(z, ax, 2 + mp_real_abs_bound_lt_2exp_si(ay)
                        - mp_real_abs_bound_lt_2exp_si(ax));
                    mp_real_sub(z, z, ay, n);
                    break;
                case 2: mp_real_mul(z, x, y, n); break;
                case 3: mp_real_mul(z, x, x, n); break;
                case 4: mp_real_div(z, x, y, n); break;
                case 5: mp_real_sqrt(z, ax, n); break;
                case 6: mp_real_rsqrt(z, ax, n); break;
                case 7: mp_real_mul_ui(z, x, c, n); break;
                case 8: mp_real_div_ui(z, x, c, n); break;
                default: mp_real_addmul_ui(z, ax, ay, c, n); break;
            }

            d = prec_rel(z, b);
            if (d > (inexact ? 16.0 : 4.0))
                TEST_FUNCTION_FAIL("op %d, %s operands, b = %wd, n = %wd: "
                    "relative radius %g eps\n", op, inexact ? "inexact" : "exact",
                    b, n, d);
        }

        mp_real_clear(x);
        mp_real_clear(y);
        mp_real_clear(ax);
        mp_real_clear(ay);
        mp_real_clear(z);
    }

    TEST_FUNCTION_END(state);
}
