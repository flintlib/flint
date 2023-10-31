/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpf_vec.h"

#define MPF_VEC_SMM_ASSOC_BITS (65)

TEST_FUNCTION_START(mpf_vec_scalar_mul_mpf, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpf *a, *b;
        mpf_t n;
        slong len = n_randint(state, 100);
        mpf_init(n);

        _flint_rand_init_gmp(state);
        mpf_urandomb(n, state->gmp_state, 100);
        if (n_randint(state, 2))
            mpf_neg(n, n);

        a = _mpf_vec_init(len, 200);
        b = _mpf_vec_init(len, 200);
        _mpf_vec_randtest(a, state, len, 200);

        _mpf_vec_scalar_mul_mpf(b, a, len, n);
        _mpf_vec_scalar_mul_mpf(a, a, len, n);

        result = (_mpf_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("aliasing failed\n");
            fflush(stdout);
            flint_abort();
        }

        _mpf_vec_clear(a, len);
        _mpf_vec_clear(b, len);
        mpf_clear(n);
    }

    /* Check that n (a + b) == na + nb */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpf *a, *b, *lhs, *rhs;
        mpf_t n;
        slong len = n_randint(state, 100);
        mpf_init(n);

        _flint_rand_init_gmp(state);
        mpf_urandomb(n, state->gmp_state, 100);
        if (n_randint(state, 2))
            mpf_neg(n, n);

        a = _mpf_vec_init(len, 200);
        b = _mpf_vec_init(len, 200);
        lhs = _mpf_vec_init(len, 300);
        rhs = _mpf_vec_init(len, 300);
        _mpf_vec_randtest(a, state, len, 200);
        _mpf_vec_randtest(b, state, len, 200);

        _mpf_vec_scalar_mul_mpf(lhs, a, len, n);
        _mpf_vec_scalar_mul_mpf(rhs, b, len, n);
        _mpf_vec_add(rhs, lhs, rhs, len);
        _mpf_vec_add(lhs, a, b, len);
        _mpf_vec_scalar_mul_mpf(lhs, lhs, len, n);

        result = (_mpf_vec_equal(lhs, rhs, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n ( a + b ) test failed\n");
            fflush(stdout);
            flint_abort();
        }

        _mpf_vec_clear(a, len);
        _mpf_vec_clear(b, len);
        _mpf_vec_clear(lhs, len);
        _mpf_vec_clear(rhs, len);
        mpf_clear(n);
    }

    /* Check that n2 * (n1 a) == (n1 * n2) a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        mpf *a, *b;
        mpf_t n1, n2, n;
        slong len = n_randint(state, 100);
        mpf_init(n1);
        mpf_init(n2);
        mpf_init(n);

        _flint_rand_init_gmp(state);
        mpf_urandomb(n1, state->gmp_state, 100);
        mpf_urandomb(n2, state->gmp_state, 100);
        if (n_randint(state, 2))
            mpf_neg(n1, n1);
        if (n_randint(state, 2))
            mpf_neg(n2, n2);

        a = _mpf_vec_init(len, 200);
        b = _mpf_vec_init(len, 200);
        _mpf_vec_randtest(a, state, len, 200);

        _mpf_vec_scalar_mul_mpf(b, a, len, n1);
        _mpf_vec_scalar_mul_mpf(b, b, len, n2);
        mpf_mul(n, n1, n2);
        _mpf_vec_scalar_mul_mpf(a, a, len, n);

        result = (_mpf_vec_approx_equal(a, b, len, MPF_VEC_SMM_ASSOC_BITS));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n2 * (n1 a) test failed\n");
            mpf_out_str(stdout, 10, 0, n1);
            flint_printf("\n");
            mpf_out_str(stdout, 10, 0, n2);
            flint_printf("\n");
            mpf_out_str(stdout, 10, 0, n);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        _mpf_vec_clear(a, len);
        _mpf_vec_clear(b, len);
        mpf_clear(n1);
        mpf_clear(n2);
        mpf_clear(n);
    }

    TEST_FUNCTION_END(state);
}
