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
#include "gr_vec.h"
#include "gr_series.h"

static char *
vec_to_str(const gr_vec_t v, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_str(out);
    GR_MUST_SUCCEED(gr_vec_write(out, v, ctx));
    return out->s;
}

TEST_FUNCTION_START(gr_vec_set_str, state)
{
    slong iter;

    /* Fixed cases over ZZ, plus malformed input. */
    {
        gr_ctx_t ZZ;
        gr_vec_t v;
        slong i;

        const char * good[] = {
            "[1, 2, 3]",
            "[]",
            "[ 7 ]",
            "[ 1 ,2,  3 ]",
            "[-5]",
            NULL,
        };
        const slong good_len[] = { 3, 0, 1, 3, 1 };

        const char * bad[] = {
            "1, 2, 3",      /* no brackets */
            "[1, 2",        /* unterminated */
            "[1, , 2]",     /* empty entry */
            "[1, 2]]",      /* trailing junk */
            "[1, 2] x",     /* trailing junk */
            "[[1, 2]",      /* unbalanced */
            NULL,
        };

        gr_ctx_init_fmpz(ZZ);
        gr_vec_init(v, 0, ZZ);

        for (i = 0; good[i] != NULL; i++)
        {
            if (gr_vec_set_str(v, good[i], 1, ZZ) != GR_SUCCESS)
                TEST_FUNCTION_FAIL("expected success: %s\n", good[i]);
            if (gr_vec_length(v, ZZ) != good_len[i])
                TEST_FUNCTION_FAIL("wrong length for %s: got %wd\n", good[i], gr_vec_length(v, ZZ));
        }

        for (i = 0; bad[i] != NULL; i++)
            if (gr_vec_set_str(v, bad[i], 1, ZZ) == GR_SUCCESS)
                TEST_FUNCTION_FAIL("expected failure: %s\n", bad[i]);

        gr_vec_clear(v, ZZ);
        gr_ctx_clear(ZZ);
    }

    /* resize = 0: the parsed length must match the pre-initialized length. */
    {
        gr_ctx_t ZZ;
        gr_vec_t v;

        gr_ctx_init_fmpz(ZZ);

        gr_vec_init(v, 3, ZZ);
        if (gr_vec_set_str(v, "[10, 20, 30]", 0, ZZ) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("resize=0 with matching length should succeed\n");
        if (gr_vec_length(v, ZZ) != 3)
            TEST_FUNCTION_FAIL("resize=0 must not change the length\n");
        if (gr_vec_set_str(v, "[1, 2]", 0, ZZ) != GR_DOMAIN)
            TEST_FUNCTION_FAIL("resize=0 with wrong length should return GR_DOMAIN\n");
        gr_vec_clear(v, ZZ);

        gr_vec_init(v, 0, ZZ);
        if (gr_vec_set_str(v, "[]", 0, ZZ) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("resize=0 empty into length-0 should succeed\n");
        gr_vec_clear(v, ZZ);

        gr_ctx_clear(ZZ);
    }

    /* Direct exercise of the shape/evaluation helpers used by gr_mat_set_str. */
    {
        gr_ctx_t ZZ;
        gr_vec_t v;
        slong count;

        gr_ctx_init_fmpz(ZZ);

        if (gr_vec_str_count_entries(&count, "[10, 20, 30]", ZZ) != GR_SUCCESS || count != 3)
            TEST_FUNCTION_FAIL("count_entries [10, 20, 30] should be 3\n");
        if (gr_vec_str_count_entries(&count, "[]", ZZ) != GR_SUCCESS || count != 0)
            TEST_FUNCTION_FAIL("count_entries [] should be 0\n");
        /* trailing characters after the matching ']' are ignored */
        if (gr_vec_str_count_entries(&count, "[1, 2], [3, 4]]", ZZ) != GR_SUCCESS || count != 2)
            TEST_FUNCTION_FAIL("count_entries should ignore trailing characters\n");

        gr_vec_init(v, 2, ZZ);
        if (_gr_vec_set_str(v->entries, "[5, 6]", 2, ZZ) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("_gr_vec_set_str with len 2 should succeed\n");
        if (_gr_vec_set_str(v->entries, "[5, 6, 7]", 2, ZZ) != GR_DOMAIN)
            TEST_FUNCTION_FAIL("_gr_vec_set_str too many entries should be GR_DOMAIN\n");
        if (_gr_vec_set_str(v->entries, "[5]", 2, ZZ) != GR_DOMAIN)
            TEST_FUNCTION_FAIL("_gr_vec_set_str too few entries should be GR_DOMAIN\n");
        gr_vec_clear(v, ZZ);

        gr_ctx_clear(ZZ);
    }

    /* A top-level comma must be detected even when entries contain parentheses:
       O(x^5) etc. exercises paren-depth tracking. */
    {
        gr_ctx_t QQ, Rx;
        gr_vec_t v;

        gr_ctx_init_fmpq(QQ);
        gr_ctx_init_gr_series(Rx, QQ, 20);
        GR_MUST_SUCCEED(gr_ctx_set_gen_name(Rx, "x"));

        gr_vec_init(v, 0, Rx);
        if (gr_vec_set_str(v, "[x + O(x^5), 1 + O(x^3)]", 1, Rx) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("failed to parse vector of series\n");
        if (gr_vec_length(v, Rx) != 2)
            TEST_FUNCTION_FAIL("expected 2 series entries, got %wd\n", gr_vec_length(v, Rx));
        gr_vec_clear(v, Rx);

        gr_ctx_clear(Rx);
        gr_ctx_clear(QQ);
    }

    /* Random round-trip over ZZ: write, then parse back. */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ZZ;
        gr_vec_t a, b;
        char * str;
        slong n;

        gr_ctx_init_fmpz(ZZ);
        n = n_randint(state, 8);

        gr_vec_init(a, n, ZZ);
        gr_vec_init(b, 0, ZZ);
        GR_MUST_SUCCEED(_gr_vec_randtest(a->entries, state, n, ZZ));

        str = vec_to_str(a, ZZ);

        if (gr_vec_set_str(b, str, 1, ZZ) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("round-trip parse failed: %s\n", str);

        if (b->length != n || _gr_vec_equal(a->entries, b->entries, n, ZZ) != T_TRUE)
            TEST_FUNCTION_FAIL("round-trip mismatch: %s\n", str);

        flint_free(str);
        gr_vec_clear(a, ZZ);
        gr_vec_clear(b, ZZ);
        gr_ctx_clear(ZZ);
    }

    TEST_FUNCTION_END(state);
}
