/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_equal_is_fmpq, state)
{
    int result;

    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        fmpq_t q, r;
        fmpz_t a, b;
        const char * vars[] = {"x","y","z"};

        fmpq_mpoly_ctx_init(ctx, 3, ORD_DEGLEX);
        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_init(q);
        fmpq_init(r);
        fmpz_init(a);
        fmpz_init(b);

        result = 1;

        fmpq_mpoly_set_str_pretty(f, "0", vars, ctx);

        fmpq_set_si(q, WORD(0), WORD(1));
        result = result && fmpq_mpoly_equal_fmpq(f, q, ctx);

        fmpz_set_si(a, WORD(0));
        result = result && fmpq_mpoly_equal_fmpz(f, a, ctx);

        result = result && fmpq_mpoly_equal_ui(f, UWORD(0), ctx);

        result = result && fmpq_mpoly_equal_si(f, WORD(0), ctx);

        result = result && fmpq_mpoly_is_fmpq(f, ctx);
                           fmpq_mpoly_get_fmpq(r, f, ctx);
        result = result && fmpq_is_zero(r);

        result = result && !fmpq_mpoly_is_gen(f, WORD(0), ctx);

        result = result && !fmpq_mpoly_is_gen(f, -WORD(1), ctx);

        result = result && (fmpq_mpoly_length(f, ctx) == WORD(0));

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("test 1\n");
            fflush(stdout);
            flint_abort();
        }

        result = 1;

        fmpq_mpoly_set_str_pretty(f, "1", vars, ctx);

        fmpq_set_si(q, WORD(1), WORD(1));
        result = result && fmpq_mpoly_equal_fmpq(f, q, ctx);

        fmpz_set_si(a, WORD(1));
        result = result && fmpq_mpoly_equal_fmpz(f, a, ctx);

        result = result && fmpq_mpoly_equal_ui(f, UWORD(1), ctx);

        result = result && fmpq_mpoly_equal_si(f, WORD(1), ctx);

        result = result && fmpq_mpoly_is_fmpq(f, ctx);
                           fmpq_mpoly_get_fmpq(r, f, ctx);
        result = result && fmpq_is_one(r);

        fmpq_mpoly_get_coeff_fmpq_monomial(r, f, f, ctx);
        result = result && fmpq_is_one(r);

        result = result && !fmpq_mpoly_is_gen(f, WORD(0), ctx);

        result = result && !fmpq_mpoly_is_gen(f, -WORD(1), ctx);

        result = result && (fmpq_mpoly_length(f, ctx) == WORD(1));

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("test 2\n");
            fflush(stdout);
            flint_abort();
        }

        result = 1;

        fmpq_mpoly_set_str_pretty(f, "2/3", vars, ctx);

        fmpq_set_si(q, WORD(2), WORD(3));
        result = result && fmpq_mpoly_equal_fmpq(f, q, ctx);

        fmpz_set_si(a, WORD(2));
        result = result && !fmpq_mpoly_equal_fmpz(f, a, ctx);

        result = result && !fmpq_mpoly_equal_ui(f, UWORD(2), ctx);

        result = result && !fmpq_mpoly_equal_si(f, WORD(2), ctx);

        result = result && fmpq_mpoly_is_fmpq(f, ctx);
                           fmpq_mpoly_get_fmpq(r, f, ctx);
        result = result && fmpq_equal(r, q);

        result = result && !fmpq_mpoly_is_gen(f, WORD(0), ctx);

        result = result && !fmpq_mpoly_is_gen(f, -WORD(1), ctx);

        result = result && (fmpq_mpoly_length(f, ctx) == WORD(1));

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("test 3\n");
            fflush(stdout);
            flint_abort();
        }

        result = 1;

        fmpq_mpoly_set_str_pretty(f, "x", vars, ctx);

        fmpq_set_si(q, WORD(1), WORD(1));
        result = result && !fmpq_mpoly_equal_fmpq(f, q, ctx);

        fmpz_set_si(a, WORD(1));
        result = result && !fmpq_mpoly_equal_fmpz(f, a, ctx);

        result = result && !fmpq_mpoly_equal_ui(f, UWORD(1), ctx);

        result = result && !fmpq_mpoly_equal_si(f, WORD(1), ctx);

        result = result && !fmpq_mpoly_is_fmpq(f, ctx);

        result = result && fmpq_mpoly_is_gen(f, WORD(0), ctx);

        result = result && fmpq_mpoly_is_gen(f, -WORD(1), ctx);

        result = result && (fmpq_mpoly_length(f, ctx) == WORD(1));

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("test 4\n");
            fflush(stdout);
            flint_abort();
        }

        result = 1;

        fmpq_mpoly_set_str_pretty(f, "2/3*x*z", vars, ctx);
        fmpq_mpoly_set_str_pretty(g, "2*x*z", vars, ctx);
        fmpq_mpoly_set_str_pretty(h, "2*x*z*y^2", vars, ctx);

        fmpq_set_si(q, WORD(2), WORD(3));
        result = result && !fmpq_mpoly_equal_fmpq(f, q, ctx);

        fmpq_mpoly_get_denominator(a, f, ctx);
        result = result && fmpz_equal_si(a, WORD(3));
        fmpq_mpoly_get_denominator(a, g, ctx);
        result = result && fmpz_equal_si(a, WORD(1));

        fmpz_set_si(a, WORD(1));
        result = result && !fmpq_mpoly_equal_fmpz(f, a, ctx);

        result = result && !fmpq_mpoly_equal_ui(f, UWORD(1), ctx);

        result = result && !fmpq_mpoly_equal_si(f, WORD(1), ctx);

        result = result && !fmpq_mpoly_is_gen(f, WORD(0), ctx);

        result = result && !fmpq_mpoly_is_gen(f, -WORD(1), ctx);

        result = result && (fmpq_mpoly_length(f, ctx) == WORD(1));

        fmpq_mpoly_get_coeff_fmpq_monomial(r, f, g, ctx);
        result = result && fmpq_equal(r, q);

        fmpq_mpoly_get_coeff_fmpq_monomial(r, f, h, ctx);
        result = result && fmpq_is_zero(r);

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("test 5\n");
            fflush(stdout);
            flint_abort();
        }

        result = 1;

        fmpq_mpoly_set_str_pretty(f, "2/3*x*z+3/5*x^999999999999999999999*y*z", vars, ctx);
        fmpq_mpoly_set_str_pretty(g, "-x*z", vars, ctx);
        fmpq_mpoly_set_str_pretty(h, "x^999999999999999999999*y*z", vars, ctx);

        fmpq_set_si(q, WORD(2), WORD(3));
        result = result && !fmpq_mpoly_equal_fmpq(f, q, ctx);

        fmpq_mpoly_get_denominator(a, f, ctx);
        result = result && fmpz_equal_si(a, WORD(15));
        fmpq_mpoly_get_denominator(a, g, ctx);
        result = result && fmpz_equal_si(a, WORD(1));

        fmpz_set_si(a, WORD(1));
        result = result && !fmpq_mpoly_equal_fmpz(f, a, ctx);

        result = result && !fmpq_mpoly_equal_ui(f, UWORD(1), ctx);

        result = result && !fmpq_mpoly_equal_si(f, WORD(1), ctx);

        result = result && !fmpq_mpoly_is_gen(f, WORD(0), ctx);

        result = result && !fmpq_mpoly_is_gen(f, -WORD(1), ctx);

        result = result && (fmpq_mpoly_length(f, ctx) != WORD(1));

        fmpq_mpoly_get_coeff_fmpq_monomial(r, f, g, ctx);
        fmpq_set_si(q, WORD(2), WORD(3));
        result = result && fmpq_equal(r, q);

        fmpq_mpoly_get_coeff_fmpq_monomial(r, f, h, ctx);
        fmpq_set_si(q, WORD(3), WORD(5));
        result = result && fmpq_equal(r, q);

        result = result && fmpq_mpoly_divides(g, h, g, ctx);
        fmpq_mpoly_get_coeff_fmpq_monomial(r, f, g, ctx);
        result = result && fmpq_is_zero(r);

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("test 6\n");
            fflush(stdout);
            flint_abort();
        }

        result = 1;

        fmpq_mpoly_set_str_pretty(f, "0", vars, ctx);

        fmpq_mpoly_set_str_pretty(g, "2*x^2", vars, ctx);
        fmpq_set_si(q, WORD(2), WORD(3));
        fmpq_mpoly_set_coeff_fmpq_monomial(f, q, g, ctx);

        fmpq_mpoly_set_str_pretty(g, "2*y^3", vars, ctx);
        fmpq_set_si(q, WORD(3), WORD(5));
        fmpq_mpoly_set_coeff_fmpq_monomial(f, q, g, ctx);

        fmpq_mpoly_set_str_pretty(g, "2*z", vars, ctx);
        fmpq_set_si(q, WORD(7), WORD(1));
        fmpq_mpoly_set_coeff_fmpq_monomial(f, q, g, ctx);

        fmpq_mpoly_set_str_pretty(g, "2*y^3", vars, ctx);
        fmpq_set_si(q, -WORD(3), WORD(5));
        fmpq_mpoly_set_coeff_fmpq_monomial(f, q, g, ctx);

        fmpq_mpoly_set_str_pretty(g, "-2", vars, ctx);
        fmpq_set_si(q, -WORD(6), WORD(11));
        fmpq_mpoly_set_coeff_fmpq_monomial(f, q, g, ctx);

        fmpq_mpoly_set_str_pretty(h, "2/3*x^2 - 3/5*y^3 + 7*z - 6/11", vars, ctx);

        result = fmpq_mpoly_equal(f, h, ctx);
        if (!result)
        {
            printf("FAIL\n");
            flint_printf("test 7\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(b);
        fmpz_clear(a);
        fmpq_clear(r);
        fmpq_clear(q);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
