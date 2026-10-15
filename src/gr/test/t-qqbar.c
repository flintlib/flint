/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"

TEST_FUNCTION_START(gr_qqbar, state)
{
    gr_ctx_t QQbar_real, QQbar;
    int flags = GR_TEST_ALWAYS_ABLE;

    gr_ctx_init_real_qqbar(QQbar_real);
    gr_test_ring(QQbar_real, 100, flags);
    gr_ctx_clear(QQbar_real);

    gr_ctx_init_complex_qqbar(QQbar);
    gr_test_ring(QQbar, 100, flags);

    /* trigonometric functions of pi x for x with a bignum numerator */
    {
        gr_ptr x, y, z;
        int status = GR_SUCCESS;

        GR_TMP_INIT3(x, y, z, QQbar);

        /* exp(pi i 10^40/7) = exp(4 pi i/7), exp(-pi i 10^40/7) = exp(-4 pi i/7) */
        status |= gr_set_str(x, "10^40/7", QQbar);
        status |= gr_exp_pi_i(y, x, QQbar);
        status |= gr_set_str(z, "4/7", QQbar);
        status |= gr_exp_pi_i(z, z, QQbar);
        if (status != GR_SUCCESS || gr_equal(y, z, QQbar) != T_TRUE)
            flint_abort();
        status |= gr_neg(x, x, QQbar);
        status |= gr_exp_pi_i(y, x, QQbar);
        status |= gr_set_str(z, "-4/7", QQbar);
        status |= gr_exp_pi_i(z, z, QQbar);
        if (status != GR_SUCCESS || gr_equal(y, z, QQbar) != T_TRUE)
            flint_abort();

        /* cos(pi 10^30/3) = -1/2, sin(pi 10^30) = 0 */
        status |= gr_set_str(x, "10^30/3", QQbar);
        status |= gr_cos_pi(y, x, QQbar);
        status |= gr_set_str(z, "-1/2", QQbar);
        if (status != GR_SUCCESS || gr_equal(y, z, QQbar) != T_TRUE)
            flint_abort();
        status |= gr_set_str(x, "10^30", QQbar);
        status |= gr_sin_pi(y, x, QQbar);
        if (status != GR_SUCCESS || gr_is_zero(y, QQbar) != T_TRUE)
            flint_abort();

        GR_TMP_CLEAR3(x, y, z, QQbar);
    }

    gr_ctx_clear(QQbar);

    TEST_FUNCTION_END(state);
}
