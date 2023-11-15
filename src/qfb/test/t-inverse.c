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

TEST_FUNCTION_START(qfb_inverse, state)
{
    int result;
    slong i;

    for (i = 1; i < 100000; i++)
    {
        fmpz_t D, L;
        qfb_t r, s, t;

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
        } while (fmpz_sgn(D) >= 0);

        qfb_reduce(r, r, D);
        qfb_inverse(s, r);
        fmpz_abs(L, D);
        fmpz_root(L, L, 4);
        qfb_nucomp(t, r, s, D, L);
        qfb_reduce(t, t, D);

        result = (qfb_is_principal_form(t, D));
        if (!result)
        {
           printf("FAIL:\n");
           qfb_print(r); printf("\n");
           qfb_print(s); printf("\n");
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
