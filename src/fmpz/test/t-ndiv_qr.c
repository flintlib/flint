/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"

TEST_FUNCTION_START(fmpz_ndiv_qr, state)
{
    slong i;
    int result;

    /* Check that a = b * nquo + nrem, and that nrem is smallest */
    for (i = 0; i < 30000 * flint_test_multiplier(); i++)
    {
        fmpz_t tmp;
        fmpz_t a, b;
        fmpz_t A, B;
        fmpz_t nquo, nrem;
        fmpz_t fquo, frem;
        fmpz_t cquo, crem;

        fmpz_init(tmp);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(A);
        fmpz_init(B);
        fmpz_init(nquo);
        fmpz_init(nrem);
        fmpz_init(fquo);
        fmpz_init(frem);
        fmpz_init(cquo);
        fmpz_init(crem);

        fmpz_randbits(a, state, n_randint(state, 200));
        fmpz_randbits(b, state, 1 + n_randint(state, 200));

        fmpz_ndiv_qr(nquo, nrem, a, b);
        {
            fmpz_set(A, a);
            fmpz_set(B, b);
            fmpz_ndiv_qr(A, B, A, B);
            if (!fmpz_equal(A, nquo) || !fmpz_equal(B, nrem))
            {
                flint_printf("FAIL: check (A, B, A, B) aliasing\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_set(A, a);
            fmpz_set(B, b);
            fmpz_ndiv_qr(B, A, A, B);
            if (!fmpz_equal(B, nquo) || !fmpz_equal(A, nrem))
            {
                flint_printf("FAIL: check (B, A, A, B) aliasing\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_fdiv_qr(fquo, frem, a, b);
        fmpz_cdiv_qr(cquo, crem, a, b);

        fmpz_set(tmp, nrem);
        fmpz_addmul(tmp, b, nquo);
        result = ( fmpz_cmp(tmp, a) == 0
                && fmpz_cmpabs(nrem, frem) <= 0
                && fmpz_cmpabs(nrem, crem) <= 0)
                && _fmpz_is_canonical(nquo) && _fmpz_is_canonical(nrem);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n");
            flint_printf("q = "); fmpz_print(nquo); flint_printf("\n");
            flint_printf("r = "); fmpz_print(nrem); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(tmp);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(A);
        fmpz_clear(B);
        fmpz_clear(nquo);
        fmpz_clear(nrem);
        fmpz_clear(fquo);
        fmpz_clear(frem);
        fmpz_clear(cquo);
        fmpz_clear(crem);
    }

    /* Check that it rounds towards zero for ties */
    for (i = 0; i < 10000 * flint_test_multiplier(); i++)
    {
        fmpz_t max;
        fmpz_t tmp;
        fmpz_t a, b;
        fmpz_t nquo, nrem;
        fmpz_t tquo, trem;
        ulong multiplier;

        fmpz_init(max);
        fmpz_init(tmp);
        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(nquo);
        fmpz_init(nrem);
        fmpz_init(tquo);
        fmpz_init(trem);

        fmpz_set_d_2exp(max, 1.0, 2*FLINT_BITS);
        fmpz_randm(a, state, max);
        fmpz_set(b, a);

        /* a -> odd * a, b -> 2 * b */
        /* Thus a // b = odd // 2 and q should have to round down */
        multiplier = n_randint(state, 50);
        multiplier += (multiplier % 2 == 0);
        fmpz_mul_ui(a, a, multiplier);
        fmpz_mul_ui(b, b, 2);

        if (n_randint(state, 2))
            fmpz_neg(a, a);
        if (n_randint(state, 2))
            fmpz_neg(b, b);
        if (fmpz_is_zero(b))
            fmpz_one(b);

        fmpz_ndiv_qr(nquo, nrem, a, b);
        fmpz_tdiv_qr(tquo, trem, a, b);

        fmpz_set(tmp, nrem);
        fmpz_addmul(tmp, b, nquo);
        result = ( fmpz_cmp(tmp, a) == 0
                && fmpz_cmp(nquo, tquo) == 0
                && fmpz_cmp(nrem, trem) == 0)
                && _fmpz_is_canonical(nquo) && _fmpz_is_canonical(nrem);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "); fmpz_print(a); flint_printf("\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n");
            flint_printf("q = "); fmpz_print(nquo); flint_printf("\n");
            flint_printf("r = "); fmpz_print(nrem); flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(max);
        fmpz_clear(tmp);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(nquo);
        fmpz_clear(nrem);
        fmpz_clear(tquo);
        fmpz_clear(trem);
    }

    TEST_FUNCTION_END(state);
}
