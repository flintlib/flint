/*
    Copyright (C) 2025 Fredrik Johansson, Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

TEST_FUNCTION_START(fmpq_poly_interpolate_multi_mod, state)
{
    slong i;

    ulong adversarial_primes[100];
    for (i = 0; i < 100; i++)
    {
        adversarial_primes[i] = n_nextprime((i == 0) ? (UWORD(1) << (FLINT_BITS - 1))
                                                     : adversarial_primes[i - 1], 1);
    }

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t P, Q;
        fmpq *x, *y;
        fmpz *tmp;
        slong j, k, n, npoints, nbad, bits;

        if (n_randint(state, 10))
        {
            npoints = n_randint(state, 10);
            n = n_randint(state, npoints + 1);
            bits = n_randint(state, 400) + 1;
        }
        else
        {
            npoints = n_randint(state, 6);
            n = n_randint(state, npoints + 1);
            bits = n_randint(state, 2000) + 1;
        }

        x = _fmpq_vec_init(npoints);
        y = _fmpq_vec_init(npoints);

        fmpq_poly_init(P);
        fmpq_poly_init(Q);

        fmpq_poly_randtest(P, state, n, bits);

        for (j = 0; j < npoints; j++)
            fmpq_set_si(x + j, -npoints/2 + j, WORD(1));

        nbad = n_randint(state, 10) ? 0 : n_randint(state, 50);
        if (n_randint(state, 2))
        {

            for (k = 0; k < nbad && npoints; k++)
                fmpz_mul_ui(P->den, P->den, adversarial_primes[n_randint(state, 100)]);
            fmpq_poly_canonicalise(P);
        }
        else
        {
            for (k = 0; k < nbad && npoints; k++)
            {
                j = n_randint(state, npoints);
                tmp = n_randint(state, 2) ? fmpq_numref(x + j) : fmpq_denref(x + j);
                fmpz_mul_ui(tmp, tmp, adversarial_primes[n_randint(state, 100)]);
                fmpq_canonicalise(x + j);
            }
        }

        for (j = 0; j < npoints; j++)
            fmpq_poly_evaluate_fmpq(y + j, P, x + j);

        fmpq_poly_interpolate_multi_mod(Q, x, y, npoints);

        if (!fmpq_poly_equal(P, Q))
        {
            flint_printf("FAIL (P != Q):\n");
            flint_printf("P  %{fmpq_poly}\n\n", P);
            flint_printf("x  %{fmpq*}\n\n", x, npoints);
            flint_printf("y  %{fmpq*}\n\n", y, npoints);
            flint_printf("Q  %{fmpq_poly}\n\n", Q);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(P);
        fmpq_poly_clear(Q);
        _fmpq_vec_clear(x, npoints);
        _fmpq_vec_clear(y, npoints);
    }

    TEST_FUNCTION_END(state);
}
