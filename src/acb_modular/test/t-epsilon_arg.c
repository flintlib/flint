/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "acb.h"
#include "acb_modular.h"
#include "arith.h"

static void
acb_modular_epsilon_arg_naive(fmpq_t arg, const psl2z_t g)
{
#define a (&g->a)
#define b (&g->b)
#define c (&g->c)
#define d (&g->d)

    if (fmpz_is_zero(c))
    {
        fmpz_set(fmpq_numref(arg), b);
        fmpz_set_ui(fmpq_denref(arg), 12);
        fmpq_canonicalise(arg);
    }
    else
    {
        fmpq_t t;
        fmpq_init(t);

        fmpz_add(fmpq_numref(arg), a, d);
        fmpz_submul_ui(fmpq_numref(arg), c, 3);
        fmpz_mul_ui(fmpq_denref(arg), c, 12);
        fmpq_canonicalise(arg);

        fmpq_dedekind_sum(t, d, c);
        fmpq_sub(arg, arg, t);

        fmpq_clear(t);
    }
#undef a
#undef b
#undef c
#undef d
}

TEST_FUNCTION_START(acb_modular_epsilon_arg, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        psl2z_t g;
        fmpq_t x, y;
        acb_t a, b;

        psl2z_init(g);
        fmpq_init(x);
        fmpq_init(y);
        acb_init(a);
        acb_init(b);

        psl2z_randtest(g, state, n_randint(state, 200));
        fmpq_randtest(x, state, 100);
        fmpq_randtest(y, state, 100);

        fmpq_set_si(x, acb_modular_epsilon_arg(g), 12);
        acb_modular_epsilon_arg_naive(y, g);

        arb_sin_cos_pi_fmpq(acb_imagref(a), acb_realref(a), x, 200);
        arb_sin_cos_pi_fmpq(acb_imagref(b), acb_realref(b), x, 200);

        if (!acb_overlaps(a, b))
        {
            flint_printf("FAIL\n");
            flint_printf("g = "); psl2z_print(g); flint_printf("\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n");
            flint_abort();
        }

        psl2z_clear(g);
        fmpq_clear(x);
        fmpq_clear(y);
        acb_clear(a);
        acb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
