/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gmpcompat.h"
#include "long_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_addmul_si)
{
    int i, result;

    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, b;
        mpz_t d, e, f, xx;
        slong x;

        fmpz_init(a);
        fmpz_init(b);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);

        fmpz_randtest(a, state, 200);
        fmpz_randtest(b, state, 200);

        fmpz_get_mpz(d, a);
        fmpz_get_mpz(e, b);
        x = z_randtest(state);

        fmpz_addmul_si(b, a, x);
        flint_mpz_init_set_si(xx, x);
        mpz_addmul(e, d, xx);

        fmpz_get_mpz(f, b);

        result = (mpz_cmp(e, f) == 0);

        if (COEFF_IS_MPZ(*b))
        {
            fmpz c = *b;
            _fmpz_demote_val(b);
            if (*b != c)
                result = 0;
        }

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, f = %Zd, xx = %Zd\n", d, e, f, xx);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);

        mpz_clear(xx);
        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
    }

    /* Check aliasing of a and b */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t a;
        mpz_t d, e, xx;
        slong x;

        fmpz_init(a);

        mpz_init(d);
        mpz_init(e);

        fmpz_randtest(a, state, 200);

        fmpz_get_mpz(d, a);
        x = n_randtest(state);

        fmpz_addmul_si(a, a, x);
        flint_mpz_init_set_si(xx, x);
        mpz_addmul(d, d, xx);

        fmpz_get_mpz(e, a);

        result = (mpz_cmp(d, e) == 0);

        if (!result)
        {
            flint_printf("FAIL:\n");
            gmp_printf("d = %Zd, e = %Zd, xx = %Zd\n", d, e, xx);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(a);

        mpz_clear(xx);
        mpz_clear(d);
        mpz_clear(e);
    }

    TEST_FUNCTION_END;
}

