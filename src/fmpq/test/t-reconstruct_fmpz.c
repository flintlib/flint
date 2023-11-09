/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"

TEST_FUNCTION_START(fmpq_reconstruct_fmpz, state)
{
    int i;

    for (i = 0; i < 1000*flint_test_multiplier(); i++)
    {
        int result;
        fmpq_t x, y;
        fmpz_t mod;
        fmpz_t res;
        mpz_t tmp;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(mod);
        fmpz_init(res);
        mpz_init(tmp);

        fmpq_randtest(x, state, 1000);

        /* Modulus m > 2*max(|n|,d)^2 */
        if (fmpz_cmpabs(&x->num, &x->den) >= 0)
            fmpz_mul(mod, &x->num, &x->num);
        else
            fmpz_mul(mod, &x->den, &x->den);
        fmpz_mul_2exp(mod, mod, 1);
        do fmpz_add_ui(mod, mod, 1);
        while (!fmpq_mod_fmpz(res, x, mod));

        result = fmpq_reconstruct_fmpz(y, res, mod);

        if (!result || !fmpq_equal(x, y))
        {
            flint_printf("FAIL: reconstruction failed\n");
            flint_printf("input = ");
            fmpq_print(x);
            flint_printf("\nmodulus = ");
            fmpz_print(mod);
            flint_printf("\nresidue = ");
            fmpz_print(res);
            flint_printf("\nreconstructed = ");
            fmpq_print(y);
            flint_printf("\nfmpq_reconstruct_fmpz return value = %d", result);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(mod);
        fmpz_clear(res);
        mpz_clear(tmp);
    }

    TEST_FUNCTION_END(state);
}
