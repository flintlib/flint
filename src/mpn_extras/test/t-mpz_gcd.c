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

/* uses _flint_mpn_gcd_test_pair and _flint_mpn_gcd_test_len (t-gcd.c) */

/* flint_mpz_gcd, flint_mpz_gcdext and flint_mpz_invert agree with GMP,
   for operands of any signs (including zero) and with aliasing */
TEST_FUNCTION_START(flint_mpz_gcd, state)
{
    slong iter;
    mpz_t a, b, g, s, t, rg, rs, rt, x, y;

    mpz_init(a);
    mpz_init(b);
    mpz_init(g);
    mpz_init(s);
    mpz_init(t);
    mpz_init(rg);
    mpz_init(rs);
    mpz_init(rt);
    mpz_init(x);
    mpz_init(y);

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        slong na, nb;
        int alias, r1, r2;

        na = _flint_mpn_gcd_test_len(state);
        nb = n_randint(state, 2) ? na : 1 + n_randint(state, na);
        _flint_mpn_gcd_test_pair(a, b, state, FLINT_MIN(na, 1000), FLINT_MIN(nb, 1000));

        if (n_randint(state, 2))
            mpz_swap(a, b);
        if (n_randint(state, 2))
            mpz_neg(a, a);
        if (n_randint(state, 2))
            mpz_neg(b, b);
        if (n_randint(state, 20) == 0)
            mpz_set_ui(a, 0);
        if (n_randint(state, 20) == 0)
            mpz_set_ui(b, 0);
        if (n_randint(state, 20) == 0)
            mpz_set_si(b, n_randint(state, 2) ? 1 : -1);

        /* gcd */
        mpz_gcd(rg, a, b);
        alias = n_randint(state, 3);
        if (alias == 0)
        {
            flint_mpz_gcd(g, a, b);
        }
        else if (alias == 1)
        {
            mpz_set(g, a);
            flint_mpz_gcd(g, g, b);
        }
        else
        {
            mpz_set(g, b);
            flint_mpz_gcd(g, a, g);
        }
        if (mpz_cmp(g, rg) != 0)
            TEST_FUNCTION_FAIL("gcd (alias %d)\na = %{mpz}\nb = %{mpz}\n", alias, a, b);

        /* gcdext */
        mpz_gcdext(rg, rs, rt, a, b);
        alias = n_randint(state, 5);
        if (alias == 0)
        {
            flint_mpz_gcdext(g, s, t, a, b);
        }
        else if (alias == 1)
        {
            /* outputs on top of the inputs */
            mpz_set(s, a);
            mpz_set(t, b);
            flint_mpz_gcdext(g, s, t, s, t);
        }
        else if (alias == 2)
        {
            mpz_set(t, a);
            mpz_set(s, b);
            flint_mpz_gcdext(g, s, t, t, s);
        }
        else if (alias == 3)
        {
            mpz_set(g, a);
            mpz_set(t, b);
            flint_mpz_gcdext(g, s, t, g, t);
        }
        else
        {
            /* only some outputs wanted */
            flint_mpz_gcdext(g, s, NULL, a, b);
            mpz_set(t, rt);
        }
        if (mpz_cmp(g, rg) != 0 || mpz_cmp(s, rs) != 0 || mpz_cmp(t, rt) != 0)
            TEST_FUNCTION_FAIL("gcdext (alias %d)\na = %{mpz}\nb = %{mpz}\n"
                "g = %{mpz}\ns = %{mpz}\nt = %{mpz}\nexpected\ng = %{mpz}\ns = %{mpz}\nt = %{mpz}\n",
                alias, a, b, g, s, t, rg, rs, rt);

        flint_mpz_gcdext(NULL, NULL, t, a, b);
        if (mpz_cmp(t, rt) != 0)
            TEST_FUNCTION_FAIL("gcdext (t only)\na = %{mpz}\nb = %{mpz}\n", a, b);

        /* invert */
        if (mpz_sgn(b) != 0)
        {
            r1 = mpz_invert(x, a, b);
            if (n_randint(state, 2))
            {
                r2 = flint_mpz_invert(y, a, b);
            }
            else
            {
                mpz_set(y, a);
                r2 = flint_mpz_invert(y, y, b);
            }
            if ((r1 != 0) != (r2 != 0) || (r1 && mpz_cmp(x, y) != 0))
                TEST_FUNCTION_FAIL("invert\na = %{mpz}\nb = %{mpz}\nx = %{mpz}\ny = %{mpz}\n",
                    a, b, x, y);
        }
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(g);
    mpz_clear(s);
    mpz_clear(t);
    mpz_clear(rg);
    mpz_clear(rs);
    mpz_clear(rt);
    mpz_clear(x);
    mpz_clear(y);

    TEST_FUNCTION_END(state);
}
