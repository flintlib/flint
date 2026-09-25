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

/* the two-argument neg and mul_2exp_si (out of place and aliased agree,
   the input is untouched), rsqrt_ui at c = 1, and the exponent safety
   predicates */

static int
misc_equal(const mp_real_t x, const mp_real_t y)
{
    return x->size == y->size && x->exp == y->exp && x->err == y->err
        && x->negative == y->negative
        && (x->size == 0 || mpn_cmp(x->d, y->d, x->size) == 0);
}

TEST_FUNCTION_START(mp_real_misc, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        mp_real_t x, x0, y, z;
        arb_t a, b;
        slong L = 1 + n_randint(state, 6), e;
        nn_ptr p = flint_malloc(L * sizeof(ulong));

        mp_real_init(x);
        mp_real_init(x0);
        mp_real_init(y);
        mp_real_init(z);
        arb_init(a);
        arb_init(b);

        flint_mpn_rrandom(p, state, L);
        _mp_real_set_mpn_2exp(x, p, L, (slong) n_randint(state, 300) - 150);
        if (n_randint(state, 2))
            mp_real_add_error_2exp_si(x, (slong) n_randint(state, 300) - 400);
        if (n_randint(state, 2))
            mp_real_neg(x, x);
        mp_real_set(x0, x);

        /* neg */
        mp_real_neg(y, x);
        if (!misc_equal(x, x0))
            TEST_FUNCTION_FAIL("neg modified its input\n");
        mp_real_set(z, x);
        mp_real_neg(z, z);
        if (!misc_equal(y, z))
            TEST_FUNCTION_FAIL("neg: aliased and out of place differ\n");
        mp_real_get_arb(a, x);
        mp_real_get_arb(b, y);
        arb_neg(b, b);
        if (!arb_equal(a, b))
            TEST_FUNCTION_FAIL("neg: wrong value\n");

        /* mul_2exp_si */
        e = (slong) n_randint(state, 1000) - 500;
        mp_real_mul_2exp_si(y, x, e);
        if (!misc_equal(x, x0))
            TEST_FUNCTION_FAIL("mul_2exp_si modified its input\n");
        mp_real_set(z, x);
        mp_real_mul_2exp_si(z, z, e);
        if (!misc_equal(y, z))
            TEST_FUNCTION_FAIL("mul_2exp_si: aliased and out of place differ\n");
        mp_real_get_arb(a, x);
        mp_real_get_arb(b, y);
        arb_mul_2exp_si(a, a, e);
        if (!arb_contains(b, a))
            TEST_FUNCTION_FAIL("mul_2exp_si: wrong value, e = %wd\n", e);

        /* safety */
        if (!mp_real_is_safe(x) || !mp_real_exp_is_safe(MP_REAL_EXP_MAX)
            || mp_real_exp_is_safe(MP_REAL_EXP_MAX + 1)
            || mp_real_exp_is_safe(-MP_REAL_EXP_MAX - 1))
            TEST_FUNCTION_FAIL("safety predicates\n");

        flint_free(p);
        mp_real_clear(x);
        mp_real_clear(x0);
        mp_real_clear(y);
        mp_real_clear(z);
        arb_clear(a);
        arb_clear(b);
    }

    {
        mp_real_t x;
        mp_real_init(x);
        mp_real_rsqrt_ui(x, 1, 5);
        if (x->size != 1 || x->d[0] != 1 || x->exp != 1 || x->err != 0 || x->negative)
            TEST_FUNCTION_FAIL("rsqrt_ui(1)\n");
        mp_real_clear(x);
    }

    TEST_FUNCTION_END(state);
}
