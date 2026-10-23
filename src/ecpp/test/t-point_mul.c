/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "ecpp.h"

/* affine coordinates of a Jacobian point */
static void
_affine(fmpz_t x, fmpz_t y, const ecpp_point_t P, const fmpz_mod_ctx_t ctx)
{
    fmpz_t zi, t;
    fmpz_init(zi);
    fmpz_init(t);
    fmpz_mod_inv(zi, P->Z, ctx);
    fmpz_mod_mul(t, zi, zi, ctx);
    fmpz_mod_mul(x, P->X, t, ctx);
    fmpz_mod_mul(t, t, zi, ctx);
    fmpz_mod_mul(y, P->Y, t, ctx);
    fmpz_clear(zi);
    fmpz_clear(t);
}

TEST_FUNCTION_START(ecpp_point_mul, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_t n, a, b, x, y, k, l, kl, acc, g, x1, y1, x2, y2, t;
        fmpz_mod_ctx_t ctx;
        ecpp_point_t P, R1, R2, R3;
        int ok;

        fmpz_init(n); fmpz_init(a); fmpz_init(b); fmpz_init(x); fmpz_init(y);
        fmpz_init(k); fmpz_init(l); fmpz_init(kl); fmpz_init(acc); fmpz_init(g);
        fmpz_init(x1); fmpz_init(y1); fmpz_init(x2); fmpz_init(y2); fmpz_init(t);

        fmpz_randprime(n, state, 20 + n_randint(state, 200), 0);
        fmpz_mod_ctx_init(ctx, n);
        ecpp_point_init(P); ecpp_point_init(R1); ecpp_point_init(R2); ecpp_point_init(R3);

        /* random curve through a random point: b = y^2 - x^3 - a x */
        fmpz_randm(a, state, n);
        fmpz_randm(x, state, n);
        fmpz_randm(y, state, n);
        fmpz_mod_mul(b, y, y, ctx);
        fmpz_mod_mul(t, x, x, ctx);
        fmpz_mod_add(t, t, a, ctx);
        fmpz_mod_mul(t, t, x, ctx);
        fmpz_mod_sub(b, b, t, ctx);

        ecpp_point_set_affine(P, x, y);
        fmpz_randtest_unsigned(k, state, 100);
        fmpz_randtest_unsigned(l, state, 100);
        fmpz_add(kl, k, l);

        /* (k + l) P == k P + l P, checked in affine coordinates */
        fmpz_one(acc);
        ecpp_point_mul(R1, P, k, a, acc, ctx);
        ecpp_point_mul(R2, P, l, a, acc, ctx);
        ecpp_point_mul(R3, P, kl, a, acc, ctx);
        fmpz_gcd(g, acc, n);
        if (!fmpz_is_one(g))
        {
            flint_printf("FAIL: acc not a unit for prime n\n");
            fflush(stdout);
            flint_abort();
        }

        ok = 1;

        /* compare (k + l) P with k P + l P added in affine coordinates */
        {
            if (!ecpp_point_is_zero(R1) && !ecpp_point_is_zero(R2) && !ecpp_point_is_zero(R3))
            {
                /* add R1 and R2 in affine coordinates */
                fmpz_t lam, x3, y3;
                fmpz_init(lam); fmpz_init(x3); fmpz_init(y3);
                _affine(x1, y1, R1, ctx);
                _affine(x2, y2, R2, ctx);
                if (fmpz_equal(x1, x2))
                {
                    if (!fmpz_equal(y1, y2))
                        ok = ecpp_point_is_zero(R3);   /* R1 = -R2 */
                    else
                    {
                        /* doubling: lam = (3 x^2 + a) / (2 y) */
                        fmpz_mod_mul(lam, x1, x1, ctx);
                        fmpz_mod_mul_ui(lam, lam, 3, ctx);
                        fmpz_mod_add(lam, lam, a, ctx);
                        fmpz_mod_add(t, y1, y1, ctx);
                        fmpz_mod_inv(t, t, ctx);
                        fmpz_mod_mul(lam, lam, t, ctx);
                    }
                }
                else
                {
                    fmpz_mod_sub(lam, y2, y1, ctx);
                    fmpz_mod_sub(t, x2, x1, ctx);
                    fmpz_mod_inv(t, t, ctx);
                    fmpz_mod_mul(lam, lam, t, ctx);
                }
                if (ok && !ecpp_point_is_zero(R3) && !(fmpz_equal(x1, x2) && !fmpz_equal(y1, y2)))
                {
                    fmpz_mod_mul(x3, lam, lam, ctx);
                    fmpz_mod_sub(x3, x3, x1, ctx);
                    fmpz_mod_sub(x3, x3, x2, ctx);
                    fmpz_mod_sub(y3, x1, x3, ctx);
                    fmpz_mod_mul(y3, y3, lam, ctx);
                    fmpz_mod_sub(y3, y3, y1, ctx);
                    _affine(x1, y1, R3, ctx);
                    ok = fmpz_equal(x1, x3) && fmpz_equal(y1, y3);
                }
                fmpz_clear(lam); fmpz_clear(x3); fmpz_clear(y3);
            }
        }

        if (!ok)
        {
            flint_printf("FAIL: (k + l) P != k P + l P\n");
            flint_printf("n = "); fmpz_print(n); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        ecpp_point_clear(P); ecpp_point_clear(R1); ecpp_point_clear(R2); ecpp_point_clear(R3);
        fmpz_mod_ctx_clear(ctx);
        fmpz_clear(n); fmpz_clear(a); fmpz_clear(b); fmpz_clear(x); fmpz_clear(y);
        fmpz_clear(k); fmpz_clear(l); fmpz_clear(kl); fmpz_clear(acc); fmpz_clear(g);
        fmpz_clear(x1); fmpz_clear(y1); fmpz_clear(x2); fmpz_clear(y2); fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
