/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq.h"
#include "fmpq_vec.h"

TEST_FUNCTION_START(fmpq_vec_get_fmpz_vec_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        fmpq * a;
        fmpz * num, * num2;
        fmpz_t den, den2, D, t;
        slong len, i;
        int shape;

        len = n_randint(state, 3) == 0 ? n_randint(state, 200) : n_randint(state, 20);
        shape = n_randint(state, 4);

        a = _fmpq_vec_init(len);
        num = _fmpz_vec_init(len);
        num2 = _fmpz_vec_init(len);
        fmpz_init(den);
        fmpz_init(den2);
        fmpz_init(D);
        fmpz_init(t);

        fmpz_randtest_not_zero(D, state, 1 + n_randint(state, 300));
        fmpz_abs(D, D);

        for (i = 0; i < len; i++)
        {
            fmpz_randtest(fmpq_numref(a + i), state, 1 + n_randint(state, 200));

            switch (shape)
            {
                case 0:
                    fmpz_randtest_not_zero(fmpq_denref(a + i), state, 1 + n_randint(state, 200));
                    fmpz_abs(fmpq_denref(a + i), fmpq_denref(a + i));
                    break;
                case 1:
                    if (n_randint(state, 4))
                        fmpz_set(fmpq_denref(a + i), D);
                    else
                        fmpz_one(fmpq_denref(a + i));
                    break;
                case 2:
                    fmpz_set_ui(fmpq_denref(a + i), 1 + n_randint(state, 30));
                    break;
                default:
                    fmpz_set_ui(t, 1 + n_randint(state, 1000));
                    fmpz_gcd(fmpq_denref(a + i), D, t);
                    if (n_randint(state, 2))
                        fmpz_mul_ui(fmpq_denref(a + i), fmpq_denref(a + i), 1 + n_randint(state, 3));
                    break;
            }

            fmpq_canonicalise(a + i);
        }

        _fmpq_vec_get_fmpz_vec_fmpz(num, den, a, len);

        /* reference */
        fmpz_one(den2);
        for (i = 0; i < len; i++)
            fmpz_lcm(den2, den2, fmpq_denref(a + i));
        for (i = 0; i < len; i++)
        {
            fmpz_divexact(t, den2, fmpq_denref(a + i));
            fmpz_mul(num2 + i, fmpq_numref(a + i), t);
        }

        if (!fmpz_equal(den, den2) || !_fmpz_vec_equal(num, num2, len))
        {
            flint_printf("FAIL: len = %wd, shape = %d\n", len, shape);
            fmpz_print(den), flint_printf("\n");
            fmpz_print(den2), flint_printf("\n");
            flint_abort();
        }

        _fmpq_vec_clear(a, len);
        _fmpz_vec_clear(num, len);
        _fmpz_vec_clear(num2, len);
        fmpz_clear(den);
        fmpz_clear(den2);
        fmpz_clear(D);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
