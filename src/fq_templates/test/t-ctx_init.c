/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2023 Albin Ahlbäck

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
#include "ulong_extras.h"

TEST_TEMPLATE_FUNCTION_START(T, ctx_init, state)
{
    int i, k, result;

    for (i = 0; i < 3 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d;
        TEMPLATE(T, ctx_t) ctx;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state,
                                   FLINT_MIN(FLINT_BITS - 1, 50)), 1));
        d = n_randint(state, 20) + 1;

        TEMPLATE(T, ctx_init)(ctx, p, d, "a");
        TEMPLATE(T, ctx_clear)(ctx);
    }

    for (i = 0; i < 3 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        slong d;
        TEMPLATE(T, ctx_t) ctx_conway, ctx_mod;

        TEMPLATE(T, t) a, b, lhs, rhs;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        TEMPLATE(T, ctx_init_conway)(ctx_conway, p, d, "a");

#ifdef FQ_H
        TEMPLATE(T, ctx_init_modulus)(ctx_mod, ctx_conway->modulus, ctx_conway->ctxp, "a");
#else
        TEMPLATE(T, ctx_init_modulus)(ctx_mod, ctx_conway->modulus, "a");
#endif

        TEMPLATE(T, init)(a, ctx_conway);
        TEMPLATE(T, init)(b, ctx_mod);
        TEMPLATE(T, init)(lhs, ctx_conway);
        TEMPLATE(T, init)(rhs, ctx_mod);

        for (k = 0; k < 30; k++)
        {
            TEMPLATE(T, randtest)(a, state, ctx_conway);
            TEMPLATE(T, set)(b, a, ctx_mod);

            TEMPLATE(T, mul)(lhs, a, a, ctx_conway);
            TEMPLATE(T, mul)(rhs, b, b, ctx_mod);

            result = (TEMPLATE(T, equal)(lhs, rhs, ctx_mod));
            if (!result)
            {
                flint_printf("FAIL:\n\n");
                flint_printf("a   = "), TEMPLATE(T, print_pretty)(a, ctx_conway), flint_printf("\n");
                flint_printf("b   = "), TEMPLATE(T, print_pretty)(b, ctx_mod), flint_printf("\n");
                flint_printf("lhs = "), TEMPLATE(T, print_pretty)(lhs, ctx_conway), flint_printf("\n");
                flint_printf("rhs = "), TEMPLATE(T, print_pretty)(rhs, ctx_mod), flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }

        }

        TEMPLATE(T, clear)(a, ctx_conway);
        TEMPLATE(T, clear)(b, ctx_mod);
        TEMPLATE(T, clear)(lhs, ctx_conway);
        TEMPLATE(T, clear)(rhs, ctx_mod);

        TEMPLATE(T, ctx_clear)(ctx_conway);
        TEMPLATE(T, ctx_clear)(ctx_mod);
    }

    TEST_FUNCTION_END(state);
}
#endif
