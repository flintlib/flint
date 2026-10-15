/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpoly.h"
#include "fmpz_mpoly_q.h"
#include "gr.h"

TEST_FUNCTION_START(gr_fmpz_mpoly_q_factored, state)
{
    gr_ctx_t ZZxy, QQxy;
    slong iter;
    int flags = 0;

    for (iter = 0; iter < 4; iter++)
    {
        gr_ctx_init_fmpz_mpoly_q_factored(ZZxy, n_randint(state, 3), mpoly_ordering_randtest(state));
        ZZxy->size_limit = 100;
        gr_test_ring(ZZxy, 100, flags);
        gr_ctx_clear(ZZxy);
    }

    /* conversions to and from fmpz_mpoly_q */
    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        gr_ptr a, b, c;
        slong nvars = n_randint(state, 4);
        ordering_t ord = mpoly_ordering_randtest(state);
        int status = GR_SUCCESS;

        gr_ctx_init_fmpz_mpoly_q_factored(ZZxy, nvars, ord);
        gr_ctx_init_fmpz_mpoly_q(QQxy, nvars, ord);

        a = gr_heap_init(ZZxy);
        b = gr_heap_init(QQxy);
        c = gr_heap_init(ZZxy);

        status |= gr_randtest(a, state, ZZxy);
        status |= gr_set_other(b, a, ZZxy, QQxy);
        status |= gr_set_other(c, b, QQxy, ZZxy);

        if (status != GR_SUCCESS || gr_equal(a, c, ZZxy) != T_TRUE)
        {
            flint_printf("FAIL: conversion\n");
            flint_printf("a = "); gr_println(a, ZZxy);
            flint_printf("b = "); gr_println(b, QQxy);
            flint_printf("c = "); gr_println(c, ZZxy);
            flint_abort();
        }

        gr_heap_clear(a, ZZxy);
        gr_heap_clear(b, QQxy);
        gr_heap_clear(c, ZZxy);

        gr_ctx_clear(ZZxy);
        gr_ctx_clear(QQxy);
    }

    TEST_FUNCTION_END(state);
}
