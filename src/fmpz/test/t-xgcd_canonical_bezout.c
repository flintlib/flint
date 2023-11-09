/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_xgcd_canonical_bezout, state)
{
    int ix, result;
    fmpz_t maxval;
    fmpz_t nd, na, nb, nf, ng;

    fmpz_init(maxval);

    /* For uniformly random distributions,
     * about half the numbers should be represented as slongs */
    fmpz_set_d_2exp(maxval, 1.0, FLINT_BITS - 1);

    fmpz_init(nd);
    fmpz_init(na);
    fmpz_init(nb);
    fmpz_init(nf);
    fmpz_init(ng);

    /* Check that xgcd(0, 0) = (0, 0, 0) */
    if (1)
    {
        fmpz_zero(nf);
        fmpz_zero(ng);

        fmpz_xgcd_canonical_bezout(nd, na, nb, nf, ng);

        result = (fmpz_is_zero(nd)
               && fmpz_is_zero(na)
               && fmpz_is_zero(nb)
               && fmpz_is_zero(nf)
               && fmpz_is_zero(ng));

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("xgcd(0, 0) is not equal to (0, 0, 0):\n");
            flint_printf("d = "), fmpz_print(nd), flint_printf("\n");
            flint_printf("a = "), fmpz_print(na), flint_printf("\n");
            flint_printf("b = "), fmpz_print(nb), flint_printf("\n");
            flint_printf("f = "), fmpz_print(nf), flint_printf("\n");
            flint_printf("g = "), fmpz_print(ng), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    /* Check that xgcd(±g, g) = (|g|, 0, sgn(g)) */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fmpz_randm(ng, state, maxval);
        if (n_randint(state, 2))
            fmpz_neg(ng, ng);
        fmpz_set(nf, ng);
        if (n_randint(state, 2))
            fmpz_neg(nf, nf);

        fmpz_xgcd_canonical_bezout(nd, na, nb, nf, ng);

        result = (fmpz_cmp_si(nd, 0) >= 0
               && fmpz_cmpabs(nd, ng) == 0
               && fmpz_is_zero(na)
               && fmpz_equal_si(nb, fmpz_sgn(ng)));

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("xgcd(+/-g, g) is not equal to (|g|, 0, sgn(g)):\n");
            flint_printf("d = "), fmpz_print(nd), flint_printf("\n");
            flint_printf("a = "), fmpz_print(na), flint_printf("\n");
            flint_printf("b = "), fmpz_print(nb), flint_printf("\n");
            flint_printf("f = "), fmpz_print(nf), flint_printf("\n");
            flint_printf("g = "), fmpz_print(ng), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    /* Check that xgcd(f, 0) = (|f|, sgn(f), 0) */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fmpz_randm(nf, state, maxval);
        if (n_randint(state, 2))
            fmpz_neg(nf, nf);
        fmpz_zero(ng);

        fmpz_xgcd_canonical_bezout(nd, na, nb, nf, ng);

        result = (fmpz_cmp_si(nd, 0) >= 0
               && fmpz_cmpabs(nd, nf) == 0
               && fmpz_equal_si(na, fmpz_sgn(nf))
               && fmpz_is_zero(nb));

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("xgcd(f, 0) is not equal to (|f|, sgn(f), 0):\n");
            flint_printf("d = "), fmpz_print(nd), flint_printf("\n");
            flint_printf("a = "), fmpz_print(na), flint_printf("\n");
            flint_printf("b = "), fmpz_print(nb), flint_printf("\n");
            flint_printf("f = "), fmpz_print(nf), flint_printf("\n");
            flint_printf("g = "), fmpz_print(ng), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    /* Check that xgcd(0, g) = (|g|, 0, sgn(g)) */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fmpz_zero(nf);
        fmpz_randm(ng, state, maxval);
        if (n_randint(state, 2))
            fmpz_neg(ng, ng);

        fmpz_xgcd_canonical_bezout(nd, na, nb, nf, ng);

        result = (fmpz_cmp_si(nd, 0) >= 0
               && fmpz_cmpabs(nd, ng) == 0
               && fmpz_is_zero(na)
               && fmpz_equal_si(nb, fmpz_sgn(ng)));

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("xgcd(0, g) is not equal to (|g|, 0, sgn(g)):\n");
            flint_printf("d = "), fmpz_print(nd), flint_printf("\n");
            flint_printf("a = "), fmpz_print(na), flint_printf("\n");
            flint_printf("b = "), fmpz_print(nb), flint_printf("\n");
            flint_printf("f = "), fmpz_print(nf), flint_printf("\n");
            flint_printf("g = "), fmpz_print(ng), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    /* Check that xgcd(f, ±1) = (1, 0, ±1) */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fmpz_randm(nf, state, maxval);
        if (n_randint(state, 2))
            fmpz_neg(nf, nf);
        fmpz_one(ng);
        if (n_randint(state, 2))
            fmpz_neg(ng, ng);

        fmpz_xgcd_canonical_bezout(nd, na, nb, nf, ng);

        result = (fmpz_is_one(nd)
               && fmpz_is_zero(na)
               && fmpz_equal_si(nb, fmpz_sgn(ng)));

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("xgcd(f, +/-1) is not equal to (1, 0, +/-1):\n");
            flint_printf("d = "), fmpz_print(nd), flint_printf("\n");
            flint_printf("a = "), fmpz_print(na), flint_printf("\n");
            flint_printf("b = "), fmpz_print(nb), flint_printf("\n");
            flint_printf("f = "), fmpz_print(nf), flint_printf("\n");
            flint_printf("g = "), fmpz_print(ng), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    /* Check that xgcd(±1, g) = (1, ±1, 0) when g is not 0 or ±1 */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fmpz_one(nf);
        if (n_randint(state, 2))
            fmpz_neg(nf, nf);
        fmpz_randm(ng, state, maxval);
        if (n_randint(state, 2))
            fmpz_neg(ng, ng);
        if (fmpz_is_zero(ng) || fmpz_is_pm1(ng))
            fmpz_set_si(ng, 3);

        fmpz_xgcd_canonical_bezout(nd, na, nb, nf, ng);

        result = (fmpz_is_one(nd)
               && fmpz_equal_si(na, fmpz_sgn(nf))
               && fmpz_is_zero(nb));

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("xgcd(+/-1, g) is not equal to (1, +/-1, 0):\n");
            flint_printf("d = "), fmpz_print(nd), flint_printf("\n");
            flint_printf("a = "), fmpz_print(na), flint_printf("\n");
            flint_printf("b = "), fmpz_print(nb), flint_printf("\n");
            flint_printf("f = "), fmpz_print(nf), flint_printf("\n");
            flint_printf("g = "), fmpz_print(ng), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }
    }

    /* Check that xgcd(±2d, g) = (d, X, sgn(g)) */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fmpz_t dsave;
        flint_bitcnt_t div2;

        fmpz_randm(nd, state, maxval);
        fmpz_init_set(dsave, nd);
        fmpz_mul_si(nf, nd, 2);
        if (n_randint(state, 2))
            fmpz_neg(nf, nf);
        fmpz_randm(ng, state, maxval);
        div2 = fmpz_val2(ng);
        fmpz_tdiv_q_2exp(ng, ng, div2); /* we still want gcd(2d, g) = d */
        fmpz_mul(ng, ng, nd);
        if (n_randint(state, 2))
            fmpz_neg(ng, ng);

        fmpz_xgcd_canonical_bezout(nd, na, nb, nf, ng);

        result = (fmpz_equal(nd, dsave)
               && fmpz_equal_si(nb, fmpz_sgn(ng)));

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("xgcd(+/-2 d, g) is not equal to (d, X, sgn(g)):\n");
            flint_printf("d = "), fmpz_print(nd), flint_printf("\n");
            flint_printf("a = "), fmpz_print(na), flint_printf("\n");
            flint_printf("b = "), fmpz_print(nb), flint_printf("\n");
            flint_printf("f = "), fmpz_print(nf), flint_printf("\n");
            flint_printf("g = "), fmpz_print(ng), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(dsave);
    }

    /* Check that xgcd(f, ±2d) = (d, sgn(f), X) */
    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        fmpz_t dsave;
        flint_bitcnt_t div2;

        fmpz_randm(nd, state, maxval);
        fmpz_init_set(dsave, nd);
        fmpz_mul_si(ng, nd, 2);
        if (n_randint(state, 2))
            fmpz_neg(ng, ng);
        fmpz_randm(nf, state, maxval);
        div2 = fmpz_val2(nf);
        fmpz_tdiv_q_2exp(nf, nf, div2); /* we still want gcd(2d, f) = d */
        fmpz_mul(nf, nf, nd);
        if (n_randint(state, 2))
            fmpz_neg(nf, nf);

        fmpz_xgcd_canonical_bezout(nd, na, nb, nf, ng);

        result = (fmpz_equal(nd, dsave)
               && fmpz_equal_si(na, fmpz_sgn(nf)));

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("xgcd(f, +/-2 d) is not equal to (d, sgn(f), X):\n");
            flint_printf("d = "), fmpz_print(nd), flint_printf("\n");
            flint_printf("a = "), fmpz_print(na), flint_printf("\n");
            flint_printf("b = "), fmpz_print(nb), flint_printf("\n");
            flint_printf("f = "), fmpz_print(nf), flint_printf("\n");
            flint_printf("g = "), fmpz_print(ng), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(dsave);
    }

    /* For the other cases, check that
     * a f + b g = d,
     * |a| < |g / (2 d)|,
     * |b| < |f / (2 d)|. */
    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++)
    {
        fmpz_t tmp;
        int aliasing;
        fmpz_init(tmp);

        fmpz_randm(nf, state, maxval);
        fmpz_randm(ng, state, maxval);

        fmpz_xgcd_canonical_bezout(nd, na, nb, nf, ng);

        aliasing = n_randint(state, 5);

        if (aliasing == 0)
        {
            fmpz_xgcd_canonical_bezout(nd, na, nb, nf, ng);
        }
        else if (aliasing == 1)
        {
            /* Test aliasing of d and f, a and g */
            fmpz_set(nd, nf);
            fmpz_set(na, ng);
            fmpz_xgcd_canonical_bezout(nd, na, nb, nd, na);
        }
        else if (aliasing == 2)
        {
            /* Test aliasing of a and f, d and g */
            fmpz_set(na, nf);
            fmpz_set(nd, ng);
            fmpz_xgcd_canonical_bezout(nd, na, nb, na, nd);
        }
        else if (aliasing == 3)
        {
            /* Test aliasing of d and f, b and g */
            fmpz_set(nd, nf);
            fmpz_set(nb, ng);
            fmpz_xgcd_canonical_bezout(nd, na, nb, nd, nb);
        }
        else
        {
            /* Test aliasing of b and f, d and g */
            fmpz_set(nb, nf);
            fmpz_set(nd, ng);
            fmpz_xgcd_canonical_bezout(nd, na, nb, nb, nd);
        }

        fmpz_mul(tmp, na, nf);
        fmpz_addmul(tmp, nb, ng);

        fmpz_mul_si(na, na, 2);
        fmpz_mul(na, na, nd);
        fmpz_mul_si(nb, nb, 2);
        fmpz_mul(nb, nb, nd);

        result = (fmpz_equal(tmp, nd)
               && fmpz_cmpabs(na, ng) < 0
               && fmpz_cmpabs(nb, nf) < 0)
               && _fmpz_is_canonical(nd)
               && _fmpz_is_canonical(na)
               && _fmpz_is_canonical(nb);

        if (!result)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("d = "), fmpz_print(nd), flint_printf("\n");
            flint_printf("a = "), fmpz_print(na), flint_printf("\n");
            flint_printf("b = "), fmpz_print(nb), flint_printf("\n");
            flint_printf("f = "), fmpz_print(nf), flint_printf("\n");
            flint_printf("g = "), fmpz_print(ng), flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(tmp);
    }

    fmpz_clear(maxval);
    fmpz_clear(nd);
    fmpz_clear(na);
    fmpz_clear(nb);
    fmpz_clear(nf);
    fmpz_clear(ng);

    TEST_FUNCTION_END(state);
}
