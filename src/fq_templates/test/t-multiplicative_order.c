/*
    Copyright (C) 2018 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "test_helpers.h"
#include "templates.h"
#include "fmpz.h"

TEST_TEMPLATE_FUNCTION_START(T, multiplicative_order, state)
{
    slong ix;
    int result;

    for (ix = 0; ix < 100 * flint_test_multiplier(); ix++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, t) a, X;
        fmpz_t ord;
        int type;

        type = n_randint(state, 2);

        /* Test big primes and small degrees */
        TEMPLATE(T, ctx_init_randtest)(ctx, state, 3);

        fmpz_init(ord);
        TEMPLATE(T, init)(a, ctx);
        TEMPLATE(T, init)(X, ctx);

        if (type == 0)
        {
            /* Check that multiplicative order is a multiple of the real one */
            TEMPLATE(T, randtest)(a, state, ctx);
            result = TEMPLATE(T, multiplicative_order)(ord, a, ctx);
            TEMPLATE(T, pow)(X, a, ord, ctx);

            if (result && !TEMPLATE(T, is_one)(X, ctx))
            {
                flint_printf("FAIL:\n\n");
                flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
                flint_printf("ord = "), fmpz_print(ord), flint_printf("\n");
                TEMPLATE(T, ctx_print)(ctx);
                flint_abort();
            }
        }
        else
        {
            /* Check multiplicative order is coherent with powering by p - 1 */
            fmpz_t size;
#if defined(FQ_NMOD_H) || defined(FQ_ZECH_H)
            ulong pm1;
#else
            fmpz_t pm1;
#endif

            fmpz_init(size);
#if !(defined(FQ_NMOD_H) || defined(FQ_ZECH_H))
            fmpz_init(pm1);
#endif
            TEMPLATE(T, gen)(X, ctx);
            if (TEMPLATE(T, is_primitive)(X, ctx))
            {
#if defined(FQ_NMOD_H) || defined(FQ_ZECH_H)
                pm1 = TEMPLATE(T, ctx_prime)(ctx) - 1;
                TEMPLATE(T, pow_ui)(a, X, pm1, ctx);
                result = TEMPLATE(T, multiplicative_order)(ord, a, ctx);
                fmpz_mul_ui(ord, ord, pm1);
#else
                fmpz_sub_ui(pm1, TEMPLATE(T, ctx_prime)(ctx), 1);
                TEMPLATE(T, pow)(a, X, pm1, ctx);
                result = TEMPLATE(T, multiplicative_order)(ord, a, ctx);
                fmpz_mul(ord, ord, pm1);
#endif

                TEMPLATE(T, ctx_order)(size, ctx);
                fmpz_sub(size, size, ord);

                if (result && !fmpz_is_one(size))
                {
                    flint_printf("FAIL:\n\n");
                    flint_printf("a = "), TEMPLATE(T, print_pretty)(a, ctx), flint_printf("\n");
                    flint_printf("ord = "), fmpz_print(ord), flint_printf("\n");
                    TEMPLATE(T, ctx_print)(ctx);
                    flint_abort();
                }
            }
#if !(defined(FQ_NMOD_H) || defined(FQ_ZECH_H))
            fmpz_clear(pm1);
#endif
            fmpz_clear(size);
        }

        fmpz_clear(ord);
        TEMPLATE(T, clear)(a, ctx);
        TEMPLATE(T, clear)(X, ctx);
        TEMPLATE(T, ctx_clear)(ctx);
    }

    TEST_FUNCTION_END(state);
}
#endif
