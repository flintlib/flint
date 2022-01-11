/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"

void
_fmpz_vec_content_naive(fmpz_t res, const fmpz * vec, slong len)
{
    if (len == 0)
        fmpz_zero(res);
    else if (len == 1)
        fmpz_abs(res, vec + 0);

    fmpz_gcd(res, vec + 0, vec + 1);
    for (len -= 2, vec += 2; len > 0; len--, vec++)
        fmpz_gcd(res, res, vec);
}

int
main(void)
{
    int ix, result;
    fmpz_t res_naive, res, multiplier;
    FLINT_TEST_INIT(state);

    flint_printf("content....");
    fflush(stdout);

    fmpz_init(res_naive);
    fmpz_init(res);
    fmpz_init(multiplier);
    
    /* Check that the naive implementation gives the same result */
    for (ix = 0; ix < 1000 * flint_test_multiplier(); ix++)
    {
        fmpz * vec;
        slong len = n_randint(state, 100);

        vec = _fmpz_vec_init(len);

        _fmpz_vec_randtest(vec, state, len, 200);
        fmpz_randtest_not_zero(multiplier, state, 100);

        _fmpz_vec_scalar_mul_fmpz(vec, vec, len, multiplier);
        _fmpz_vec_content(res, vec, len);
        _fmpz_vec_content_naive(res_naive, vec, len);

        result = fmpz_equal(res_naive, res);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("res        = "), fmpz_print(res), flint_printf("\n");
            flint_printf("res_naive  = "), fmpz_print(res_naive), flint_printf("\n");
            flint_printf("multiplier = "), fmpz_print(multiplier), flint_printf("\n");
            flint_printf("vec        = "); _fmpz_vec_scalar_divexact_fmpz(vec, vec, len, multiplier); _fmpz_vec_print(vec, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        _fmpz_vec_clear(vec, len);
    }

    fmpz_clear(res_naive);
    fmpz_clear(res);
    fmpz_clear(multiplier);

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
