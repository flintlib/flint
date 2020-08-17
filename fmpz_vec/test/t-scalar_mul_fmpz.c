/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("scalar_mul_fmpz....");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        fmpz_t n;
        slong len = n_randint(state, 100);
        fmpz_init(n);
        fmpz_randtest(n, state, 100);
        if (n_randint(state, 2))
            fmpz_neg(n, n);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_scalar_mul_fmpz(b, a, len, n);
        _fmpz_vec_scalar_mul_fmpz(a, a, len, n);

        result = (_fmpz_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        fmpz_clear(n);
    }

    /* Check that n (a + b) == na + nb */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *lhs, *rhs;
        fmpz_t n;
        slong len = n_randint(state, 100);
        fmpz_init(n);
        fmpz_randtest(n, state, 100);
        if (n_randint(state, 2))
            fmpz_neg(n, n);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        lhs = _fmpz_vec_init(len);
        rhs = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);
        _fmpz_vec_randtest(b, state, len, 200);

        _fmpz_vec_scalar_mul_fmpz(lhs, a, len, n);
        _fmpz_vec_scalar_mul_fmpz(rhs, b, len, n);
        _fmpz_vec_add(rhs, lhs, rhs, len);
        _fmpz_vec_add(lhs, a, b, len);
        _fmpz_vec_scalar_mul_fmpz(lhs, lhs, len, n);

        result = (_fmpz_vec_equal(lhs, rhs, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            _fmpz_vec_print(lhs, len), flint_printf("\n\n");
            _fmpz_vec_print(rhs, len), flint_printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        _fmpz_vec_clear(lhs, len);
        _fmpz_vec_clear(rhs, len);
        fmpz_clear(n);
    }

    /* Check that n2 * (n1 a) == (n1 * n2) a */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        fmpz_t n1, n2, n;
        slong len = n_randint(state, 100);
        fmpz_init(n1);
        fmpz_init(n2);
        fmpz_init(n);
        fmpz_randtest(n1, state, 100);
        fmpz_randtest(n2, state, 100);
        if (n_randint(state, 2))
            fmpz_neg(n1, n1);
        if (n_randint(state, 2))
            fmpz_neg(n2, n2);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_scalar_mul_fmpz(b, a, len, n1);
        _fmpz_vec_scalar_mul_fmpz(b, b, len, n2);
        fmpz_mul(n, n1, n2);
        _fmpz_vec_scalar_mul_fmpz(a, a, len, n);

        result = (_fmpz_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        fmpz_clear(n1);
        fmpz_clear(n2);
        fmpz_clear(n);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
