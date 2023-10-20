/*
    Copyright (C) 2009, 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_scalar_submul_fmpz, state)
{
    int i, result;

    /* Compare with fmpz_vec_scalar_submul_si */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *c;
        slong len, n;
        fmpz_t n1;
        len = n_randint(state, 100);
        n = (slong) n_randbits(state, FLINT_BITS - 1);
        if (n_randint(state, 2))
            n = -n;
        fmpz_init(n1);
        fmpz_set_si(n1, n);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        c = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);
        _fmpz_vec_randtest(b, state, len, 200);
        _fmpz_vec_set(c, b, len);

        _fmpz_vec_scalar_submul_fmpz(b, a, len, n1);
        _fmpz_vec_scalar_submul_si(c, a, len, n);

        result = (_fmpz_vec_equal(c, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(c, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n1);
        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        _fmpz_vec_clear(c, len);
    }

    /* Compute a different way */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *c, *d;
        slong len = n_randint(state, 100);
        fmpz_t n1;
        fmpz_init(n1);
        fmpz_randtest(n1, state, 200);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        c = _fmpz_vec_init(len);
        d = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);
        _fmpz_vec_randtest(b, state, len, 200);
        _fmpz_vec_set(c, b, len);

        _fmpz_vec_scalar_submul_fmpz(b, a, len, n1);
        _fmpz_vec_scalar_mul_fmpz(d, a, len, n1);
        _fmpz_vec_sub(c, c, d, len);

        result = (_fmpz_vec_equal(c, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(c, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(n1);
        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        _fmpz_vec_clear(c, len);
        _fmpz_vec_clear(d, len);
    }

    TEST_FUNCTION_END(state);
}
