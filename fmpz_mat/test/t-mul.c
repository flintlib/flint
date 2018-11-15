/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "ulong_extras.h"

void fmpz_mat_mul_1(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);
void fmpz_mat_mul_2a(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);
void fmpz_mat_mul_2b(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);
void fmpz_mat_mul_4(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B);

int main(void)
{
    fmpz_mat_t A, B, C, D;
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("mul....");
    fflush(stdout);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong m, n, k;
        slong abits, bbits, bits;

        if (n_randint(state, 10) == 0)
        {
            m = n_randint(state, 50);
            n = n_randint(state, 50);
            k = n_randint(state, 50);
        }
        else
        {
            m = n_randint(state, 8);
            n = n_randint(state, 8);
            k = n_randint(state, 8);
        }

        fmpz_mat_init(A, m, n);
        fmpz_mat_init(B, n, k);
        fmpz_mat_init(C, m, k);
        fmpz_mat_init(D, m, k);

        fmpz_mat_randtest(A, state, n_randint(state, 200) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, 200) + 1);

        abits = fmpz_mat_max_bits(A);
        bbits = fmpz_mat_max_bits(B);
        abits = FLINT_ABS(abits);
        bbits = FLINT_ABS(bbits);
        bits = abits + bbits + FLINT_BIT_COUNT(n) + 1;

        /* Make sure noise in the output is ok */
        fmpz_mat_randtest(C, state, n_randint(state, 200) + 1);

        fmpz_mat_mul(C, A, B);
        fmpz_mat_mul_classical_inline(D, A, B);

        if (!fmpz_mat_equal(C, D))
        {
            flint_printf("FAIL: results not equal\n\n");
            fmpz_mat_print(A); flint_printf("\n\n");
            fmpz_mat_print(B); flint_printf("\n\n");
            fmpz_mat_print(C); flint_printf("\n\n");
            fmpz_mat_print(D); flint_printf("\n\n");
            flint_abort();
        }

        if (bits <= FLINT_BITS - 2)
        {
            fmpz_mat_mul_1(C, A, B);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal (mul_1)\n\n");
                fmpz_mat_print(A); flint_printf("\n\n");
                fmpz_mat_print(B); flint_printf("\n\n");
                fmpz_mat_print(C); flint_printf("\n\n");
                fmpz_mat_print(D); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (abits <= FLINT_BITS - 2 && bbits <= FLINT_BITS - 2 && bits <= 2 * FLINT_BITS - 1)
        {
            fmpz_mat_mul_2a(C, A, B);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal (mul_2a)\n\n");
                fmpz_mat_print(A); flint_printf("\n\n");
                fmpz_mat_print(B); flint_printf("\n\n");
                fmpz_mat_print(C); flint_printf("\n\n");
                fmpz_mat_print(D); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (abits <= FLINT_BITS - 2 && bbits <= FLINT_BITS - 2)
        {
            fmpz_mat_mul_2b(C, A, B);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal (mul_2b)\n\n");
                fmpz_mat_print(A); flint_printf("\n\n");
                fmpz_mat_print(B); flint_printf("\n\n");
                fmpz_mat_print(C); flint_printf("\n\n");
                fmpz_mat_print(D); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (abits <= 2 * FLINT_BITS && bbits <= 2 * FLINT_BITS && bits <= 4 * FLINT_BITS - 1)
        {
            fmpz_mat_mul_4(C, A, B);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL: results not equal (mul_4)\n\n");
                fmpz_mat_print(A); flint_printf("\n\n");
                fmpz_mat_print(B); flint_printf("\n\n");
                fmpz_mat_print(C); flint_printf("\n\n");
                fmpz_mat_print(D); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (n == k)
        {
            fmpz_mat_mul(A, A, B);

            if (!fmpz_mat_equal(A, C))
            {
                flint_printf("FAIL: aliasing failed\n");
                flint_abort();
            }
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
