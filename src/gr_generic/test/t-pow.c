/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_generic.h"

/* todo: have a proper interface to test a given powering function */
int gr_test_pow_ui_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags);
int gr_test_pow_ui_base_scalar_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags);
int gr_test_pow_ui_base_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags);
int gr_test_pow_ui_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags);
int gr_test_pow_fmpz_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags);

TEST_FUNCTION_START(gr_generic_pow, state)
{
    gr_ctx_t ZZn;
    ulong n;
    fmpz_t m;
    int status = GR_SUCCESS;

    for (n = 0; n < 1000 * flint_test_multiplier(); n++)
    {
        fmpz_init(m);
        fmpz_randtest_not_zero(m, state, 100);
        fmpz_abs(m, m);

        gr_ctx_init_fmpz_mod(ZZn, m);

        ZZn->methods[GR_METHOD_POW_UI] = (gr_funcptr) gr_generic_pow_ui_binexp;
        ZZn->methods[GR_METHOD_POW_FMPZ] = (gr_funcptr) gr_generic_pow_fmpz_binexp;

        status |= gr_test_pow_ui_exponent_addition(ZZn, state, 0);
        status |= gr_test_pow_ui_base_scalar_multiplication(ZZn, state, 0);
        status |= gr_test_pow_ui_base_multiplication(ZZn, state, 0);
        status |= gr_test_pow_ui_aliasing(ZZn, state, 0);
        status |= gr_test_pow_fmpz_exponent_addition(ZZn, state, 0);

        ZZn->methods[GR_METHOD_POW_UI] = (gr_funcptr) gr_generic_pow_ui_sliding;
        ZZn->methods[GR_METHOD_POW_FMPZ] = (gr_funcptr) gr_generic_pow_fmpz_sliding;

        status |= gr_test_pow_ui_exponent_addition(ZZn, state, 0);
        status |= gr_test_pow_ui_base_scalar_multiplication(ZZn, state, 0);
        status |= gr_test_pow_ui_base_multiplication(ZZn, state, 0);
        status |= gr_test_pow_ui_aliasing(ZZn, state, 0);
        status |= gr_test_pow_fmpz_exponent_addition(ZZn, state, 0);

        if (status & GR_TEST_FAIL)
        {
            fflush(stdout);
            flint_abort();
        }

        gr_ctx_clear(ZZn);
        fmpz_clear(m);
    }

    TEST_FUNCTION_END(state);
}
