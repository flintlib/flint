/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, frobenius, state)
{
    int i, result;

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, t) a, b, c;
        slong e;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, set)(b, a, ctx);
        e = n_randint(state, 10) % TEMPLATE(T, ctx_degree)(ctx);

        TEMPLATE(T, frobenius)(c, b, e, ctx);
        TEMPLATE(T, frobenius)(b, b, e, ctx);

        result = (TEMPLATE(T, equal)(b, c, ctx));
        if (!result)
        {
            flint_printf("FAIL (alias):\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            flint_printf("e = %wd\n", e);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check sigma^e(x) == x^{p^e}  */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, b, c;
        slong e;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(c, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        e = n_randint(state, 10) % TEMPLATE(T, ctx_degree)(ctx);

        TEMPLATE(T, frobenius)(b, a, e, ctx);
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, TEMPLATE(T, ctx_prime)(ctx), e);
            TEMPLATE(T, pow)(c, a, t, ctx);
            fmpz_clear(t);
        }

        result = (TEMPLATE(T, equal)(b,c,ctx));
        if (!result)
        {
            flint_printf("FAIL (sigma^e(x) = x^{p^e}):\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("c = "), TEMPLATE(T, print_pretty)(c, ctx), flint_printf("\n");
            flint_printf("e = %wd\n", e);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(c, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check sigma^e(x + y) = sigma^e(x) + sigma^e(y) */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, t) a, b, s, s1, s2, lhs, rhs;
        slong e;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(s, ctx);
        TEMPLATE(T, init)(s1, ctx);
        TEMPLATE(T, init)(s2, ctx);
        TEMPLATE(T, init)(lhs, ctx);
        TEMPLATE(T, init)(rhs, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, randtest)(b, state, ctx);
        e = n_randint(state, 10) % TEMPLATE(T, ctx_degree)(ctx);

        TEMPLATE(T, add)(s, a, b, ctx);
        TEMPLATE(T, frobenius)(lhs, s, e, ctx);
        TEMPLATE(T, frobenius)(s1, a, e, ctx);
        TEMPLATE(T, frobenius)(s2, b, e, ctx);
        TEMPLATE(T, add)(rhs, s1, s2, ctx);

        result = (TEMPLATE(T, equal)(lhs, rhs, ctx));
        if (!result)
        {
            flint_printf("FAIL (sigma(a+b) = sigma(a) + sigma(b)):\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("s = "), TEMPLATE(T, print_pretty)(s, ctx), flint_printf("\n");
            flint_printf("s1 = "), TEMPLATE(T, print_pretty)(s1, ctx), flint_printf("\n");
            flint_printf("s2 = "), TEMPLATE(T, print_pretty)(s2, ctx), flint_printf("\n");
            flint_printf("lhs = "), TEMPLATE(T, print_pretty)(lhs, ctx), flint_printf("\n");
            flint_printf("rhs = "), TEMPLATE(T, print_pretty)(rhs, ctx), flint_printf("\n");
            flint_printf("e = %wd\n", e);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(s, ctx);
        TEMPLATE(T, clear)(s1, ctx);
        TEMPLATE(T, clear)(s2, ctx);
        TEMPLATE(T, clear)(lhs, ctx);
        TEMPLATE(T, clear)(rhs, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    /* Check sigma^e(x * y) = sigma^e(x) * sigma^e(y) on Zq */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        TEMPLATE(T, ctx_t) ctx;

        TEMPLATE(T, t) a, b, s, s1, s2, lhs, rhs;
        slong e;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(b, ctx);
        TEMPLATE(T, init)(s, ctx);
        TEMPLATE(T, init)(s1, ctx);
        TEMPLATE(T, init)(s2, ctx);
        TEMPLATE(T, init)(lhs, ctx);
        TEMPLATE(T, init)(rhs, ctx);

        TEMPLATE(T, randtest)(a, state, ctx);
        TEMPLATE(T, randtest)(b, state, ctx);
        e = n_randint(state, 10) % TEMPLATE(T, ctx_degree)(ctx);

        TEMPLATE(T, mul)(s, a, b, ctx);
        TEMPLATE(T, frobenius)(lhs, s, e, ctx);
        TEMPLATE(T, frobenius)(s1, a, e, ctx);
        TEMPLATE(T, frobenius)(s2, b, e, ctx);
        TEMPLATE(T, mul)(rhs, s1, s2, ctx);

        result = (TEMPLATE(T, equal)(lhs, rhs, ctx));
        if (!result)
        {
            flint_printf("FAIL (sigma(a*b) = sigma(a) * sigma(b)):\n\n");
            flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
            flint_printf("b = "), TEMPLATE(T, print_pretty)(b, ctx), flint_printf("\n");
            flint_printf("s = "), TEMPLATE(T, print_pretty)(s, ctx), flint_printf("\n");
            flint_printf("s1 = "), TEMPLATE(T, print_pretty)(s1, ctx), flint_printf("\n");
            flint_printf("s2 = "), TEMPLATE(T, print_pretty)(s2, ctx), flint_printf("\n");
            flint_printf("lhs = "), TEMPLATE(T, print_pretty)(lhs, ctx), flint_printf("\n");
            flint_printf("rhs = "), TEMPLATE(T, print_pretty)(rhs, ctx), flint_printf("\n");
            flint_printf("e = %wd\n", e);
            fflush(stdout);
            flint_abort();
        }

        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(b, ctx);
        TEMPLATE(T, clear)(s, ctx);
        TEMPLATE(T, clear)(s1, ctx);
        TEMPLATE(T, clear)(s2, ctx);
        TEMPLATE(T, clear)(lhs, ctx);
        TEMPLATE(T, clear)(rhs, ctx);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
