/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "fmpz.h"

TEST_FUNCTION_START(nmod_poly_bit_pack, state)
{
    int i, result;

    /* Check aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n;
        ulong bits;
        mp_ptr mpn;

        do
        {
            n = n_randtest_not_zero(state);
        } while (n == 1);
        bits = 2 * FLINT_BIT_COUNT(n) + n_randint(state, FLINT_BITS);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        do
        {
            nmod_poly_randtest(a, state, n_randint(state, 100));
        } while (a->length == 0);

        mpn =
            flint_malloc(sizeof(mp_limb_t) *
                   ((bits * a->length - 1) / FLINT_BITS + 1));

        _nmod_poly_bit_pack(mpn, a->coeffs, a->length, bits);
        nmod_poly_fit_length(b, a->length);
        _nmod_poly_bit_unpack(b->coeffs, a->length, mpn, bits, a->mod);
        b->length = a->length;

        result = (nmod_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        flint_free(mpn);
    }

    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        fmpz_t f;
        nmod_poly_t A, B;
        slong b;
        mp_limb_t n;

        do
        {
            n = n_randtest_not_zero(state);
        } while (n == 1);

        fmpz_init(f);
        nmod_poly_init(A, n);
        nmod_poly_init(B, n);

        nmod_poly_randtest(A, state, 1+n_randint(state,100));

        b = FLINT_BIT_COUNT(n) + n_randint(state, FLINT_BITS);

        nmod_poly_bit_pack(f, A, b);
        nmod_poly_bit_unpack(B, f, b);

        if (!nmod_poly_equal(A, B))
        {
            mpz_t zz;
            flint_printf("FAIL:\n");
            flint_printf("INPUT: ");
            nmod_poly_print(A);
            flint_printf("\n");
            mpz_init(zz); fmpz_get_mpz(zz, f);
            flint_printf("PACKED: ");
            mpz_out_str(stdout, 2, zz);
            flint_printf("\n");
            flint_printf("OUTPUT: ");
            nmod_poly_print(B);
            flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(f);
        nmod_poly_clear(A);
        nmod_poly_clear(B);
    }

    TEST_FUNCTION_END(state);
}
