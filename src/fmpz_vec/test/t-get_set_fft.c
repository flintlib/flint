/*
    Copyright (C) 2009, 2011 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fft.h"
#include "fmpz_vec.h"

TEST_FUNCTION_START(fmpz_vec_get_set_fft, state)
{
    int i, result;

     /* convert back and forth and compare */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz * a, * b;
        flint_bitcnt_t bits;
        slong len, limbs;
        mp_limb_t ** ii, * ptr;
        slong i, bt;

        bits = n_randint(state, 300) + 1;
        len = n_randint(state, 300) + 1;
        limbs = 2*((bits - 1)/FLINT_BITS + 1);

        ii = flint_malloc((len + len*(limbs + 1))*sizeof(mp_limb_t));
        ptr = (mp_limb_t *) ii + len;
        for (i = 0; i < len; i++, ptr += (limbs + 1))
           ii[i] = ptr;

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, bits);

        fmpz_one(a + 0);
        fmpz_mul_2exp(a + 0, a + 0, FLINT_BITS*limbs - 1);

        _fmpz_vec_get_fft(ii, a, limbs, len);
        bt = _fmpz_vec_max_bits(a, len);
        for (i = 0; i < len; i++)
           mpn_normmod_2expp1(ii[i], limbs);
        _fmpz_vec_set_fft(b, len, ii, limbs, bt < 0);

        result = (_fmpz_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(ii);
        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
    }

     /* convert back and forth unsigned and compare */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz * a, * b;
        flint_bitcnt_t bits;
        slong len, limbs;
        mp_limb_t ** ii, * ptr;
        slong i, bt;

        bits = n_randint(state, 300) + 1;
        len = n_randint(state, 300) + 1;
        limbs = 2*((bits - 1)/FLINT_BITS + 1);

        ii = flint_malloc((len + len*(limbs + 1))*sizeof(mp_limb_t));
        ptr = (mp_limb_t *) ii + len;
        for (i = 0; i < len; i++, ptr += (limbs + 1))
           ii[i] = ptr;

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest_unsigned(a, state, len, bits);

        _fmpz_vec_get_fft(ii, a, limbs, len);
        bt = _fmpz_vec_max_bits(a, len);
        _fmpz_vec_set_fft(b, len, ii, limbs, bt < 0);

        result = (_fmpz_vec_equal(a, b, len));
        if (!result)
        {
            flint_printf("FAIL:\n");
            _fmpz_vec_print(a, len), flint_printf("\n\n");
            _fmpz_vec_print(b, len), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        flint_free(ii);
        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
    }

    TEST_FUNCTION_END(state);
}
