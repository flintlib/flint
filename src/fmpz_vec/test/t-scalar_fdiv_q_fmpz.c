/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_scalar_fdiv_q_fmpz, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *c;
        fmpz_t n;
        mpz_t d, e, f, m;
        slong i;
        slong len = n_randint(state, 100);

        fmpz_init(n);
        fmpz_randtest_not_zero(n, state, 100);
        if (n_randint(state, 2))
            fmpz_neg(n, n);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        c = _fmpz_vec_init(len);

        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_set(b, a, len);

        _fmpz_vec_scalar_fdiv_q_fmpz(c, a, len, n);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(m);

        for (i = 0; i < len; i++)
        {
            fmpz_get_mpz(m, n);
            fmpz_get_mpz(d, b + i);

            mpz_fdiv_q(e, d, m);

            fmpz_get_mpz(f, c + i);

            result = (mpz_cmp(f, e) == 0);
            if (!result)
            {
                flint_printf("FAIL:\n");
                gmp_printf("d = %Zd, m = %Zd, e = %Zd, f = %Zd\n", d, m, e, f);
                fflush(stdout);
                flint_abort();
            }
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        _fmpz_vec_clear(c, len);

        fmpz_clear(n);
        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(m);
    }

    /* Test aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        fmpz_t n;
        mpz_t d, e, f, m;
        slong i;
        slong len = n_randint(state, 100);

        fmpz_init(n);
        fmpz_randtest_not_zero(n, state, 100);
        if (n_randint(state, 2))
            fmpz_neg(n, n);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);

        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_set(b, a, len);

        _fmpz_vec_scalar_fdiv_q_fmpz(a, a, len, n);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(m);

        for (i = 0; i < len; i++)
        {
            fmpz_get_mpz(m, n);
            fmpz_get_mpz(d, b + i);

            mpz_fdiv_q(e, d, m);

            fmpz_get_mpz(f, a + i);

            result = (mpz_cmp(f, e) == 0);
            if (!result)
            {
                flint_printf("FAIL:\n");
                gmp_printf("d = %Zd, m = %Zd, e = %Zd, f = %Zd\n", d, m, e, f);
                fflush(stdout);
                flint_abort();
            }
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);

        fmpz_clear(n);
        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(m);
    }

    TEST_FUNCTION_END(state);
}
