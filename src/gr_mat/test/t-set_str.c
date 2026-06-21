/*
    Copyright (C) 2025 FLINT contributors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr.h"
#include "gr_mat.h"

static char *
mat_to_str(const gr_mat_t m, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_str(out);
    GR_MUST_SUCCEED(gr_mat_write(out, m, ctx));
    return out->s;
}

TEST_FUNCTION_START(gr_mat_set_str, state)
{
    slong iter;

    /* Fixed cases over ZZ. */
    {
        gr_ctx_t ZZ;
        gr_mat_t m, expected;

        gr_ctx_init_fmpz(ZZ);

        /* 2x2 value check, including the newline form produced by the printer */
        {
            const char * forms[] = {
                "[[1, 2], [3, 4]]",
                "[[1, 2],\n[3, 4]]",
                "[ [1, 2] , [3, 4] ]",
                NULL,
            };
            slong i;

            gr_mat_init(expected, 2, 2, ZZ);
            GR_MUST_SUCCEED(gr_set_si(gr_mat_entry_ptr(expected, 0, 0, ZZ), 1, ZZ));
            GR_MUST_SUCCEED(gr_set_si(gr_mat_entry_ptr(expected, 0, 1, ZZ), 2, ZZ));
            GR_MUST_SUCCEED(gr_set_si(gr_mat_entry_ptr(expected, 1, 0, ZZ), 3, ZZ));
            GR_MUST_SUCCEED(gr_set_si(gr_mat_entry_ptr(expected, 1, 1, ZZ), 4, ZZ));

            gr_mat_init(m, 0, 0, ZZ);
            for (i = 0; forms[i] != NULL; i++)
            {
                if (gr_mat_set_str(m, forms[i], 1, ZZ) != GR_SUCCESS)
                    TEST_FUNCTION_FAIL("expected success: %s\n", forms[i]);
                if (gr_mat_equal(m, expected, ZZ) != T_TRUE)
                    TEST_FUNCTION_FAIL("wrong value for %s\n", forms[i]);
            }
            gr_mat_clear(m, ZZ);
            gr_mat_clear(expected, ZZ);
        }

        /* Shapes. */
        {
            gr_mat_init(m, 0, 0, ZZ);

            if (gr_mat_set_str(m, "[[1, 2, 3]]", 1, ZZ) != GR_SUCCESS ||
                gr_mat_nrows(m, ZZ) != 1 || gr_mat_ncols(m, ZZ) != 3)
                TEST_FUNCTION_FAIL("expected a 1x3 matrix\n");

            if (gr_mat_set_str(m, "[[1], [2], [3]]", 1, ZZ) != GR_SUCCESS ||
                gr_mat_nrows(m, ZZ) != 3 || gr_mat_ncols(m, ZZ) != 1)
                TEST_FUNCTION_FAIL("expected a 3x1 matrix\n");

            if (gr_mat_set_str(m, "[]", 1, ZZ) != GR_SUCCESS ||
                gr_mat_nrows(m, ZZ) != 0 || gr_mat_ncols(m, ZZ) != 0)
                TEST_FUNCTION_FAIL("expected a 0x0 matrix\n");

            /* n x 0 is unambiguous: n empty rows. */
            if (gr_mat_set_str(m, "[[], [], []]", 1, ZZ) != GR_SUCCESS ||
                gr_mat_nrows(m, ZZ) != 3 || gr_mat_ncols(m, ZZ) != 0)
                TEST_FUNCTION_FAIL("expected a 3x0 matrix\n");

            gr_mat_clear(m, ZZ);
        }

        /* Rejected input. */
        {
            const char * bad[] = {
                "[[1, 2], [3]]",     /* non-rectangular */
                "[[1, 2], 3]",       /* row not bracketed */
                "[1, 2, 3]",         /* a vector, not a matrix */
                "[[1, 2]",           /* unbalanced */
                "[[1, 2]] junk",     /* trailing junk */
                NULL,
            };
            slong i;

            gr_mat_init(m, 0, 0, ZZ);
            for (i = 0; bad[i] != NULL; i++)
                if (gr_mat_set_str(m, bad[i], 1, ZZ) == GR_SUCCESS)
                    TEST_FUNCTION_FAIL("expected failure: %s\n", bad[i]);
            gr_mat_clear(m, ZZ);
        }

        /* resize = 0: the parsed shape must match the pre-initialized shape,
           with "[]" compatible with any 0 x c. */
        {
            gr_mat_init(m, 2, 2, ZZ);
            if (gr_mat_set_str(m, "[[1, 2], [3, 4]]", 0, ZZ) != GR_SUCCESS)
                TEST_FUNCTION_FAIL("resize=0 matching shape should succeed\n");
            if (gr_mat_nrows(m, ZZ) != 2 || gr_mat_ncols(m, ZZ) != 2)
                TEST_FUNCTION_FAIL("resize=0 must not change the shape\n");
            if (gr_mat_set_str(m, "[[1, 2, 3]]", 0, ZZ) != GR_DOMAIN)
                TEST_FUNCTION_FAIL("resize=0 wrong shape should return GR_DOMAIN\n");
            if (gr_mat_set_str(m, "[]", 0, ZZ) != GR_DOMAIN)
                TEST_FUNCTION_FAIL("resize=0 [] into a 2x2 should return GR_DOMAIN\n");
            gr_mat_clear(m, ZZ);

            /* "[]" is compatible with any 0 x c when resize = 0. */
            gr_mat_init(m, 0, 5, ZZ);
            if (gr_mat_set_str(m, "[]", 0, ZZ) != GR_SUCCESS)
                TEST_FUNCTION_FAIL("resize=0 [] into a 0x5 should succeed\n");
            if (gr_mat_nrows(m, ZZ) != 0 || gr_mat_ncols(m, ZZ) != 5)
                TEST_FUNCTION_FAIL("resize=0 [] must leave the 0x5 shape intact\n");
            gr_mat_clear(m, ZZ);

            /* n x 0 is exact, not the degenerate case. */
            gr_mat_init(m, 3, 0, ZZ);
            if (gr_mat_set_str(m, "[[], [], []]", 0, ZZ) != GR_SUCCESS)
                TEST_FUNCTION_FAIL("resize=0 3x0 matching should succeed\n");
            if (gr_mat_set_str(m, "[[], []]", 0, ZZ) != GR_DOMAIN)
                TEST_FUNCTION_FAIL("resize=0 2x0 into 3x0 should return GR_DOMAIN\n");
            gr_mat_clear(m, ZZ);
        }

        gr_ctx_clear(ZZ);
    }
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ZZ;
        gr_mat_t a, b;
        char * str;
        slong r, c;

        gr_ctx_init_fmpz(ZZ);
        r = n_randint(state, 5);
        c = n_randint(state, 5);
        if (r == 0 || c == 0)
            r = c = 0;   /* a 0xc or rx0 matrix prints as "[]"; normalize */

        gr_mat_init(a, r, c, ZZ);
        gr_mat_init(b, 0, 0, ZZ);
        GR_MUST_SUCCEED(gr_mat_randtest(a, state, ZZ));

        str = mat_to_str(a, ZZ);

        if (gr_mat_set_str(b, str, 1, ZZ) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("round-trip parse failed: %s\n", str);

        if (gr_mat_equal(a, b, ZZ) != T_TRUE)
            TEST_FUNCTION_FAIL("round-trip mismatch: %s\n", str);

        flint_free(str);
        gr_mat_clear(a, ZZ);
        gr_mat_clear(b, ZZ);
        gr_ctx_clear(ZZ);
    }

    TEST_FUNCTION_END(state);
}
