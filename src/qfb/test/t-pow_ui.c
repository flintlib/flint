/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "qfb.h"

TEST_FUNCTION_START(qfb_pow_ui, state)
{
    int result;
    slong i, j;

    for (i = 1; i < 1000; i++)
    {
        fmpz_t D, L;
        qfb_t r, s, t;
        ulong exp;

        fmpz_init(D);
        fmpz_init(L);
        qfb_init(r);
        qfb_init(s);
        qfb_init(t);

        do
        {
           fmpz_randtest_unsigned(r->a, state, 100);
           if (fmpz_is_zero(r->a))
              fmpz_set_ui(r->a, 1);

           fmpz_randtest(r->b, state, 100);
           fmpz_randtest(r->c, state, 100);

           qfb_discriminant(D, r);
        } while (fmpz_sgn(D) >= 0 || !qfb_is_primitive(r));

        fmpz_abs(L, D);
        fmpz_root(L, L, 4);

        qfb_reduce(r, r, D);

        exp = n_randint(state, 1000);

        qfb_pow_ui(s, r, D, exp);
        qfb_reduce(s, s, D);

        qfb_principal_form(t, D);
        for (j = 0; j < exp; j++)
        {
           qfb_nucomp(t, t, r, D, L);
           qfb_reduce(t, t, D);
        }

        result = (qfb_equal(s, t));
        if (!result)
        {
           printf("FAIL:\n");
           flint_printf("exp = %wu\n", exp);
           qfb_print(r); printf("\n");
           qfb_print(s); printf("\n");
           qfb_print(t); printf("\n");
           flint_abort();
        }

        fmpz_clear(D);
        fmpz_clear(L);
        qfb_clear(r);
        qfb_clear(s);
        qfb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
