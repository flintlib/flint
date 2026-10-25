/*
    Copyright (C) 2026 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <string.h>
#include "fmpq.h"
#include "arf.h"
#include "gr.h"
#include "decimal.h"

/* number of significant digits in a string written by arf_get_str_rnd */
static slong
_num_digits(const char * s)
{
    slong n = 0;
    int leading = 1;

    for ( ; *s != '\0' && *s != 'e'; s++)
    {
        if (*s >= '0' && *s <= '9')
        {
            if (*s != '0')
                leading = 0;
            if (!leading)
                n++;
        }
    }

    return n;
}

TEST_FUNCTION_START(arf_get_str_rnd, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t y, z;
        arf_t x;
        fmpq_t q;
        char * s;
        char * s2;
        slong d;
        arf_rnd_t rnd;

        arf_init(x);
        fmpq_init(q);

        arf_randtest_special(x, state, 2 + n_randint(state, 200), 2 + n_randint(state, 12));
        d = 1 + n_randint(state, 60);
        rnd = n_randint(state, 5);

        s = arf_get_str_rnd(x, d, rnd);

        if (!arf_is_finite(x))
        {
            s2 = arf_get_str(x, d);

            if (strcmp(s, s2))
            {
                flint_printf("FAIL (special)\n");
                arf_print(x); flint_printf("\n");
                flint_printf("%s\n%s\n", s, s2);
                flint_abort();
            }

            flint_free(s2);
        }
        else
        {
            gr_ctx_init_decfloat(ctx, d, 0);
            decfloat_init(y, ctx);
            decfloat_init(z, ctx);

            /* the printed value must be the correctly rounded value */
            arf_get_fmpq(q, x);
            GR_MUST_SUCCEED(decfloat_set_round_fmpq_reference(y, q, d, rnd, ctx));
            GR_MUST_SUCCEED(decfloat_set_round_str(z, s, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));

            if (decfloat_equal(y, z, ctx) != T_TRUE ||
                (!arf_is_zero(x) && _num_digits(s) != d))
            {
                flint_printf("FAIL (value)\n");
                arf_print(x); flint_printf("\n");
                flint_printf("d = %wd, rnd = %d\n", d, (int) rnd);
                flint_printf("%s\n", s);
                s2 = decfloat_get_str(y, ctx);
                flint_printf("%s\n", s2);
                flint_abort();
            }

            /* same format as arf_get_str whenever that prints the same
               value (it is only faithful, so it may differ by 1 ulp) */
            s2 = arf_get_str(x, d);
            GR_MUST_SUCCEED(decfloat_set_round_str(z, s2, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));

            if (decfloat_equal(y, z, ctx) == T_TRUE && strcmp(s, s2))
            {
                flint_printf("FAIL (format)\n");
                arf_print(x); flint_printf("\n");
                flint_printf("d = %wd, rnd = %d\n", d, (int) rnd);
                flint_printf("%s\n%s\n", s, s2);
                flint_abort();
            }

            flint_free(s2);

            decfloat_clear(y, ctx);
            decfloat_clear(z, ctx);
            gr_ctx_clear(ctx);
        }

        flint_free(s);
        arf_clear(x);
        fmpq_clear(q);
    }

    /* huge exponents, where an exact rational reference is infeasible */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        decfloat_t y, z;
        arf_t x;
        fmpz_t e;
        char * s;
        char * s2;
        slong d;
        arf_rnd_t rnd;

        arf_init(x);
        fmpz_init(e);

        arf_randtest_not_zero(x, state, 2 + n_randint(state, 200), 8);
        fmpz_randbits(e, state, 63 + n_randint(state, 140));
        arf_mul_2exp_fmpz(x, x, e);
        d = 1 + n_randint(state, 60);
        rnd = n_randint(state, 5);

        s = arf_get_str_rnd(x, d, rnd);

        gr_ctx_init_decfloat(ctx, d, 0);
        decfloat_init(y, ctx);
        decfloat_init(z, ctx);

        /* reference: round a faithful approximation with 30 guard digits;
           this can only fail if x is within a relative distance 10^(-d-29)
           of a rounding boundary, which does not happen in practice */
        s2 = arf_get_str(x, d + 30);
        GR_MUST_SUCCEED(decfloat_set_round_str(y, s2, d, rnd, ctx));
        flint_free(s2);
        GR_MUST_SUCCEED(decfloat_set_round_str(z, s, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));

        if (decfloat_equal(y, z, ctx) != T_TRUE || _num_digits(s) != d)
        {
            flint_printf("FAIL (huge value)\n");
            arf_print(x); flint_printf("\n");
            flint_printf("d = %wd, rnd = %d\n", d, (int) rnd);
            flint_printf("%s\n", s);
            s2 = decfloat_get_str(y, ctx);
            flint_printf("%s\n", s2);
            flint_abort();
        }

        s2 = arf_get_str(x, d);
        GR_MUST_SUCCEED(decfloat_set_round_str(z, s2, DECIMAL_PREC_EXACT, DECIMAL_RND_DOWN, ctx));

        if (decfloat_equal(y, z, ctx) == T_TRUE && strcmp(s, s2))
        {
            flint_printf("FAIL (huge format)\n");
            arf_print(x); flint_printf("\n");
            flint_printf("d = %wd, rnd = %d\n", d, (int) rnd);
            flint_printf("%s\n%s\n", s, s2);
            flint_abort();
        }

        flint_free(s2);
        flint_free(s);

        decfloat_clear(y, ctx);
        decfloat_clear(z, ctx);
        gr_ctx_clear(ctx);
        arf_clear(x);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
