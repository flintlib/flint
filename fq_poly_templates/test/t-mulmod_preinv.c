/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("mulmod_preinv....");
    fflush(stdout);

    /* Aliasing res and a */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, res, t, f, finv;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (finv, ctx);
        TEMPLATE(T, poly_init) (res, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 50), ctx);
        do
        {
            TEMPLATE(T, poly_randtest) (f, state, n_randint(state, 50), ctx);
        } while (TEMPLATE(T, poly_is_zero) (f, ctx));
        if (a->length >= f->length)
            TEMPLATE(T, poly_rem) (a, a, f, ctx);
        if (b->length >= f->length)
            TEMPLATE(T, poly_rem) (b, b, f, ctx);

        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        TEMPLATE(T, poly_mulmod_preinv) (res, a, b, f, finv, ctx);
        TEMPLATE(T, poly_mulmod_preinv) (a, a, b, f, finv, ctx);

        result = (TEMPLATE(T, poly_equal) (res, a, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("b:\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n");
            TEMPLATE(T, poly_print) (res, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (finv, ctx);
        TEMPLATE(T, poly_clear) (res, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Aliasing res and b */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, res, t, f, finv;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (finv, ctx);
        TEMPLATE(T, poly_init) (res, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 50), ctx);
        do
        {
            TEMPLATE(T, poly_randtest) (f, state, n_randint(state, 50), ctx);
        } while (TEMPLATE(T, poly_is_zero) (f, ctx));
        if (a->length >= f->length)
            TEMPLATE(T, poly_rem) (a, a, f, ctx);
        if (b->length >= f->length)
            TEMPLATE(T, poly_rem) (b, b, f, ctx);

        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        TEMPLATE(T, poly_mulmod_preinv) (res, a, b, f, finv, ctx);
        TEMPLATE(T, poly_mulmod_preinv) (b, a, b, f, finv, ctx);

        result = (TEMPLATE(T, poly_equal) (res, b, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("b:\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n");
            TEMPLATE(T, poly_print) (res, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (finv, ctx);
        TEMPLATE(T, poly_clear) (res, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Aliasing res and f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, res, t, f, finv;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (finv, ctx);
        TEMPLATE(T, poly_init) (res, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 50), ctx);
        do
        {
            TEMPLATE(T, poly_randtest) (f, state, n_randint(state, 50), ctx);
        } while (TEMPLATE(T, poly_is_zero) (f, ctx));
        if (a->length >= f->length)
            TEMPLATE(T, poly_rem) (a, a, f, ctx);
        if (b->length >= f->length)
            TEMPLATE(T, poly_rem) (b, b, f, ctx);

        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        TEMPLATE(T, poly_mulmod_preinv) (res, a, b, f, finv, ctx);
        TEMPLATE(T, poly_mulmod_preinv) (f, a, b, f, finv, ctx);

        result = (TEMPLATE(T, poly_equal) (res, f, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("b:\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n");
            TEMPLATE(T, poly_print) (res, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (finv, ctx);
        TEMPLATE(T, poly_clear) (res, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* Aliasing res and finv */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, res, t, f, finv;
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (finv, ctx);
        TEMPLATE(T, poly_init) (res, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 50), ctx);
        do
        {
            TEMPLATE(T, poly_randtest) (f, state, n_randint(state, 50), ctx);
        } while (TEMPLATE(T, poly_is_zero) (f, ctx));
        if (a->length >= f->length)
            TEMPLATE(T, poly_rem) (a, a, f, ctx);
        if (b->length >= f->length)
            TEMPLATE(T, poly_rem) (b, b, f, ctx);

        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        TEMPLATE(T, poly_mulmod_preinv) (res, a, b, f, finv, ctx);
        TEMPLATE(T, poly_mulmod_preinv) (finv, a, b, f, finv, ctx);

        result = (TEMPLATE(T, poly_equal) (res, finv, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("b:\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("finv:\n");
            TEMPLATE(T, poly_print) (finv, ctx), flint_printf("\n\n");
            flint_printf("res1:\n");
            TEMPLATE(T, poly_print) (res, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (finv, ctx);
        TEMPLATE(T, poly_clear) (res, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    /* No aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, poly_t) a, b, res1, res2, t, f, finv;

        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, ctx_randtest) (ctx, state);

        TEMPLATE(T, poly_init) (a, ctx);
        TEMPLATE(T, poly_init) (b, ctx);
        TEMPLATE(T, poly_init) (f, ctx);
        TEMPLATE(T, poly_init) (finv, ctx);
        TEMPLATE(T, poly_init) (res1, ctx);
        TEMPLATE(T, poly_init) (res2, ctx);
        TEMPLATE(T, poly_init) (t, ctx);

        TEMPLATE(T, poly_randtest) (a, state, n_randint(state, 50), ctx);
        TEMPLATE(T, poly_randtest) (b, state, n_randint(state, 50), ctx);
        do
        {
            TEMPLATE(T, poly_randtest) (f, state, n_randint(state, 50), ctx);
        } while (TEMPLATE(T, poly_is_zero) (f, ctx));
        if (a->length >= f->length)
            TEMPLATE(T, poly_rem) (a, a, f, ctx);
        if (b->length >= f->length)
            TEMPLATE(T, poly_rem) (b, b, f, ctx);

        TEMPLATE(T, poly_reverse) (finv, f, f->length, ctx);
        TEMPLATE(T, poly_inv_series_newton) (finv, finv, f->length, ctx);

        TEMPLATE(T, poly_mulmod_preinv) (res1, a, b, f, finv, ctx);

        TEMPLATE(T, poly_mul) (res2, a, b, ctx);
        TEMPLATE(T, poly_divrem) (t, res2, res2, f, ctx);

        result = (TEMPLATE(T, poly_equal) (res1, res2, ctx));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a:\n");
            TEMPLATE(T, poly_print) (a, ctx), flint_printf("\n\n");
            flint_printf("b:\n");
            TEMPLATE(T, poly_print) (b, ctx), flint_printf("\n\n");
            flint_printf("f:\n");
            TEMPLATE(T, poly_print) (f, ctx), flint_printf("\n\n");
            flint_printf("res1:\n");
            TEMPLATE(T, poly_print) (res1, ctx), flint_printf("\n\n");
            flint_printf("res2:\n");
            TEMPLATE(T, poly_print) (res2, ctx), flint_printf("\n\n");
            abort();
        }

        TEMPLATE(T, poly_clear) (a, ctx);
        TEMPLATE(T, poly_clear) (b, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_clear) (finv, ctx);
        TEMPLATE(T, poly_clear) (res1, ctx);
        TEMPLATE(T, poly_clear) (res2, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
        TEMPLATE(T, ctx_clear) (ctx);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}


#endif
