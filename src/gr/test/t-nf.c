/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "gr.h"

TEST_FUNCTION_START(gr_nf, state)
{
    gr_ctx_t QQa;
    fmpq_poly_t f;
    fmpz_poly_t g;
    slong iter;
    int flags = GR_TEST_ALWAYS_ABLE;

    fmpq_poly_init(f);
    fmpz_poly_init(g);

    for (iter = 0; iter < 30; iter++)
    {
        do
        {
            fmpz_poly_randtest_irreducible(g, state, 2 + n_randint(state, 5), 1 + n_randint(state, 10));
        } while (g->length < 2);

        fmpq_poly_set_fmpz_poly(f, g);
        fmpq_poly_scalar_div_ui(f, f, 1 + n_randtest(state) % 256);

        gr_ctx_init_nf(QQa, f);
        gr_test_ring(QQa, 100, flags);
        gr_ctx_clear(QQa);
    }

    fmpq_poly_clear(f);
    fmpz_poly_clear(g);

    TEST_FUNCTION_END(state);
}
