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
#include "fmpz.h"

TEST_FUNCTION_START(flint_mpz_div, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        mpz_t a, b, q, r, q2, r2, t, x, y;
        mpz_ptr qq, rr;
        mpz_srcptr aa, bb;
        slong abits, bbits;
        int op, alias;

        mpz_init(a); mpz_init(b); mpz_init(q); mpz_init(r); mpz_init(q2); mpz_init(r2); mpz_init(t);
        mpz_init(x); mpz_init(y);

        /* operands of a hundred thousand bits cost a thousand times more
           than the small ones, so they are only used in a few iterations */
        if (n_randint(state, 50) == 0)
        {
            switch (n_randint(state, 3))
            {
                case 0: bbits = 1 + n_randint(state, 100000); abits = n_randint(state, 200000); break;
                case 1: bbits = 1 + n_randint(state, 5000); abits = n_randint(state, 150000); break;
                default: bbits = 1 + n_randint(state, 100000); abits = bbits + n_randint(state, 3000); break;
            }
        }
        else
        {
            switch (n_randint(state, 3))
            {
                case 0: bbits = 1 + n_randint(state, 200); abits = n_randint(state, 400); break;
                case 1: bbits = 1 + n_randint(state, 20000); abits = bbits + n_randint(state, 20000); break;
                default: bbits = 1 + n_randint(state, 8 * FLINT_BITS); abits = n_randint(state, 20 * FLINT_BITS); break;
            }
        }

        {
            fmpz_t fa, fb;
            fmpz_init(fa); fmpz_init(fb);
            fmpz_randbits(fa, state, abits);
            fmpz_randbits(fb, state, bbits);
            fmpz_get_mpz(a, fa);
            fmpz_get_mpz(b, fb);
            fmpz_clear(fa); fmpz_clear(fb);
        }
        if (mpz_sgn(b) == 0)
            mpz_set_ui(b, 1);
        if (n_randint(state, 2)) mpz_neg(a, a);
        if (n_randint(state, 2)) mpz_neg(b, b);
        if (n_randint(state, 5) == 0)
        {
            mpz_mul(a, a, b);
            if (n_randint(state, 2))
                mpz_add_ui(a, a, n_randint(state, 3));
        }
        if (n_randint(state, 20) == 0)
            mpz_set(a, b);

        /* garbage in the outputs */
        {
            fmpz_t g;
            fmpz_init(g);
            fmpz_randbits(g, state, n_randint(state, 300));
            fmpz_get_mpz(q, g);
            fmpz_randbits(g, state, n_randint(state, 300));
            fmpz_get_mpz(r, g);
            fmpz_clear(g);
        }

        op = n_randint(state, 11);

        /* alias: 0 none, 1 q = a, 2 q = b, 3 r = a, 4 r = b, 5 a = b */
        alias = n_randint(state, 6);
        if (op >= 3 && op <= 5 && alias >= 3 && alias <= 4)
            alias -= 2;
        if (op >= 6 && op <= 9 && (alias == 1 || alias == 2))
            alias += 2;

        if (alias == 5)
        {
            /* a is also the divisor */
            if (mpz_sgn(a) == 0)
                mpz_set_ui(a, 1);
            mpz_set(b, a);
        }
        if (mpz_sgn(b) == 0)
            mpz_set_ui(b, 1);

        qq = q; rr = r; aa = a; bb = b;
        if (alias == 1) { mpz_set(q, a); aa = q; }
        if (alias == 2) { mpz_set(q, b); bb = q; }
        if (alias == 3) { mpz_set(r, a); aa = r; }
        if (alias == 4) { mpz_set(r, b); bb = r; }
        if (alias == 5) { bb = aa; }

        if (op <= 2)
        {
            if (op == 0) { mpz_tdiv_qr(q2, r2, a, b); flint_mpz_tdiv_qr(qq, rr, aa, bb); }
            if (op == 1) { mpz_fdiv_qr(q2, r2, a, b); flint_mpz_fdiv_qr(qq, rr, aa, bb); }
            if (op == 2) { mpz_cdiv_qr(q2, r2, a, b); flint_mpz_cdiv_qr(qq, rr, aa, bb); }
        }
        else if (op <= 5)
        {
            if (op == 3) { mpz_tdiv_q(q2, a, b); flint_mpz_tdiv_q(qq, aa, bb); }
            if (op == 4) { mpz_fdiv_q(q2, a, b); flint_mpz_fdiv_q(qq, aa, bb); }
            if (op == 5) { mpz_cdiv_q(q2, a, b); flint_mpz_cdiv_q(qq, aa, bb); }
            mpz_set_ui(r, 0); mpz_set_ui(r2, 0);
        }
        else if (op <= 9)
        {
            if (op == 6) { mpz_tdiv_r(r2, a, b); flint_mpz_tdiv_r(rr, aa, bb); }
            if (op == 7) { mpz_fdiv_r(r2, a, b); flint_mpz_fdiv_r(rr, aa, bb); }
            if (op == 8) { mpz_cdiv_r(r2, a, b); flint_mpz_cdiv_r(rr, aa, bb); }
            if (op == 9) { mpz_mod(r2, a, b); flint_mpz_mod(rr, aa, bb); }
            mpz_set_ui(q, 0); mpz_set_ui(q2, 0);
        }
        else
        {
            /* exact division: q = a b / b */
            mpz_mul(x, a, b);
            mpz_set(y, b);
            mpz_divexact(q2, x, y);
            if (alias == 1 || alias == 3) { mpz_set(q, x); flint_mpz_divexact(q, q, y); }
            else if (alias == 2 || alias == 4) { mpz_set(q, y); flint_mpz_divexact(q, x, q); }
            else if (alias == 5) { flint_mpz_divexact(q, y, y); mpz_set_ui(q2, 1); }
            else flint_mpz_divexact(q, x, y);
            mpz_set_ui(r, 0); mpz_set_ui(r2, 0);
        }

        if (mpz_cmp(q, q2) != 0 || mpz_cmp(r, r2) != 0)
            TEST_FUNCTION_FAIL("op = %d, alias = %d, abits = %wd, bbits = %wd\n", op, alias, abits, bbits);

        /* square root: always for small operands, rarely for large ones
           (four square roots per iteration are the dominant cost there) */
        mpz_abs(a, a);
        if (mpz_size(a) > 100 && n_randint(state, 20) != 0)
            goto cleanup;

        mpz_sqrtrem(q2, r2, a);
        mpz_sqrt(t, a);
        alias = n_randint(state, 3);
        if (alias == 1) { mpz_set(q, a); flint_mpz_sqrtrem(q, r, q); }
        else if (alias == 2) { mpz_set(r, a); flint_mpz_sqrtrem(q, r, r); }
        else flint_mpz_sqrtrem(q, r, a);
        if (mpz_cmp(q, q2) != 0 || mpz_cmp(r, r2) != 0)
            TEST_FUNCTION_FAIL("sqrtrem: abits = %wd, alias = %d\n", abits, alias);
        if (n_randint(state, 2)) { mpz_set(q, a); flint_mpz_sqrt(q, q); }
        else flint_mpz_sqrt(q, a);
        if (mpz_cmp(q, t) != 0)
            TEST_FUNCTION_FAIL("sqrt: abits = %wd\n", abits);

cleanup:
        mpz_clear(a); mpz_clear(b); mpz_clear(q); mpz_clear(r); mpz_clear(q2); mpz_clear(r2); mpz_clear(t);
        mpz_clear(x); mpz_clear(y);
    }

    TEST_FUNCTION_END(state);
}
