/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_default_mat.h"

TEST_FUNCTION_START(fq_default_mat_set_fmpz_mat, state)
{
    int i, result;

    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_mat_t fq_mat1, fq_mat2;
        fmpz_mat_t mat;
        fmpz_t p;
        slong rows, cols;

        rows = n_randint(state, 20);
        cols = n_randint(state, 20);

        fmpz_init(p);

        fmpz_set_ui(p, 3);

        fq_default_ctx_init(ctx, p, 3, "x");

        fmpz_mat_init(mat, rows, cols);
        fmpz_mat_one(mat);

        fq_default_mat_init(fq_mat1, rows, cols, ctx);
        fq_default_mat_init(fq_mat2, rows, cols, ctx);

        fq_default_mat_set_fmpz_mat(fq_mat1, mat, ctx);

        fq_default_mat_one(fq_mat2, ctx);

        result = fq_default_mat_equal(fq_mat1, fq_mat2, ctx);
        if (!result)
        {
            flint_printf("FAIL: matrices not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(mat);

        fq_default_mat_clear(fq_mat1, ctx);
        fq_default_mat_clear(fq_mat2, ctx);

        fq_default_ctx_clear(ctx);

        fmpz_clear(p);
    }

    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_mat_t fq_mat1, fq_mat2;
        fmpz_mat_t mat;
        fmpz_t p;
        slong rows, cols;

        rows = n_randint(state, 20);
        cols = n_randint(state, 20);

        fmpz_init(p);

        fmpz_set_ui(p, 3);

        fq_default_ctx_init(ctx, p, 16, "x");

        fmpz_mat_init(mat, rows, cols);
        fmpz_mat_one(mat);

        fq_default_mat_init(fq_mat1, rows, cols, ctx);
        fq_default_mat_init(fq_mat2, rows, cols, ctx);

        fq_default_mat_set_fmpz_mat(fq_mat1, mat, ctx);

        fq_default_mat_one(fq_mat2, ctx);

        result = fq_default_mat_equal(fq_mat1, fq_mat2, ctx);
        if (!result)
        {
            flint_printf("FAIL: matrices not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(mat);

        fq_default_mat_clear(fq_mat1, ctx);
        fq_default_mat_clear(fq_mat2, ctx);

        fq_default_ctx_clear(ctx);

        fmpz_clear(p);
    }

    for (i = 0; i < 1 * flint_test_multiplier(); i++)
    {
        fq_default_ctx_t ctx;
        fq_default_mat_t fq_mat1, fq_mat2;
        fmpz_mat_t mat;
        fmpz_t p;
        slong rows, cols;

        rows = n_randint(state, 20);
        cols = n_randint(state, 20);

        fmpz_init(p);

        fmpz_set_str(p, "73786976294838206473", 10);

        fq_default_ctx_init(ctx, p, 1, "x");

        fmpz_mat_init(mat, rows, cols);
        fmpz_mat_one(mat);

        fq_default_mat_init(fq_mat1, rows, cols, ctx);
        fq_default_mat_init(fq_mat2, rows, cols, ctx);

        fq_default_mat_set_fmpz_mat(fq_mat1, mat, ctx);

        fq_default_mat_one(fq_mat2, ctx);

        result = fq_default_mat_equal(fq_mat1, fq_mat2, ctx);
        if (!result)
        {
            flint_printf("FAIL: matrices not equal\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mat_clear(mat);

        fq_default_mat_clear(fq_mat1, ctx);
        fq_default_mat_clear(fq_mat2, ctx);

        fq_default_ctx_clear(ctx);

        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
