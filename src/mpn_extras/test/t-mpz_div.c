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

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        mpz_t a, b, q, r, q2, r2, t;
        slong abits, bbits;
        int op, alias;

        mpz_init(a); mpz_init(b); mpz_init(q); mpz_init(r); mpz_init(q2); mpz_init(r2); mpz_init(t);

        switch (n_randint(state, 4))
        {
            case 0: bbits = 1 + n_randint(state, 200); abits = n_randint(state, 400); break;
            case 1: bbits = 1 + n_randint(state, 100000); abits = n_randint(state, 200000); break;
            case 2: bbits = 1 + n_randint(state, 5000); abits = n_randint(state, 150000); break;
            default: bbits = 1 + n_randint(state, 100000); abits = bbits + n_randint(state, 3000); break;
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

        op = n_randint(state, 10);
        alias = n_randint(state, 3);

        if (op <= 2)
        {
            if (op == 0) { mpz_tdiv_qr(q2, r2, a, b); }
            if (op == 1) { mpz_fdiv_qr(q2, r2, a, b); }
            if (op == 2) { mpz_cdiv_qr(q2, r2, a, b); }
            if (alias == 1) mpz_set(q, a); else if (alias == 2) mpz_set(r, b);
            {
                mpz_srcptr aa = (alias == 1) ? q : a, bb = (alias == 2) ? r : b;
                if (op == 0) flint_mpz_tdiv_qr(q, r, aa, bb);
                if (op == 1) flint_mpz_fdiv_qr(q, r, aa, bb);
                if (op == 2) flint_mpz_cdiv_qr(q, r, aa, bb);
            }
        }
        else if (op <= 5)
        {
            if (op == 3) { mpz_tdiv_q(q2, a, b); flint_mpz_tdiv_q(q, a, b); }
            if (op == 4) { mpz_fdiv_q(q2, a, b); flint_mpz_fdiv_q(q, a, b); }
            if (op == 5) { mpz_cdiv_q(q2, a, b); flint_mpz_cdiv_q(q, a, b); }
            mpz_set_ui(r, 0); mpz_set_ui(r2, 0);
        }
        else if (op <= 8)
        {
            if (op == 6) { mpz_tdiv_r(r2, a, b); flint_mpz_tdiv_r(r, a, b); }
            if (op == 7) { mpz_fdiv_r(r2, a, b); flint_mpz_fdiv_r(r, a, b); }
            if (op == 8) { mpz_mod(r2, a, b); flint_mpz_mod(r, a, b); }
            mpz_set_ui(q, 0); mpz_set_ui(q2, 0);
        }
        else
        {
            mpz_mul(a, a, b);
            mpz_divexact(q2, a, b);
            if (alias == 1) { mpz_set(q, a); flint_mpz_divexact(q, q, b); }
            else flint_mpz_divexact(q, a, b);
            mpz_set_ui(r, 0); mpz_set_ui(r2, 0);
        }

        if (mpz_cmp(q, q2) != 0 || mpz_cmp(r, r2) != 0)
            TEST_FUNCTION_FAIL("op = %d, alias = %d, abits = %wd, bbits = %wd\n", op, alias, abits, bbits);

        /* square root */
        mpz_abs(a, a);
        mpz_sqrtrem(q2, r2, a);
        if (alias == 1) { mpz_set(q, a); flint_mpz_sqrtrem(q, r, q); }
        else flint_mpz_sqrtrem(q, r, a);
        mpz_sqrt(t, a);
        if (mpz_cmp(q, q2) != 0 || mpz_cmp(r, r2) != 0)
            TEST_FUNCTION_FAIL("sqrtrem: abits = %wd\n", abits);
        flint_mpz_sqrt(q, a);
        if (mpz_cmp(q, t) != 0)
            TEST_FUNCTION_FAIL("sqrt: abits = %wd\n", abits);

        mpz_clear(a); mpz_clear(b); mpz_clear(q); mpz_clear(r); mpz_clear(q2); mpz_clear(r2); mpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
