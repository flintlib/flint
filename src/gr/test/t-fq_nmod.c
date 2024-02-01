/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr.h"

TEST_FUNCTION_START(gr_fq_nmod, state)
{
    gr_ctx_t Fq;
    ulong p;
    slong d;
    slong iter;
    int flags = GR_TEST_ALWAYS_ABLE;

    for (iter = 0; iter < 30; iter++)
    {
        p = n_randprime(state, 2 + n_randint(state, FLINT_BITS - 1), 0);
        d = 1 + n_randint(state, 5);
        gr_ctx_init_fq_nmod(Fq, p, d, "a");
        gr_test_ring(Fq, 100, flags);
        gr_ctx_clear(Fq);
    }

    TEST_FUNCTION_END(state);
}
