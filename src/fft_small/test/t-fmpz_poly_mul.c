/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fft_small.h"

TEST_FUNCTION_START(_fmpz_poly_mul_mid_mpn_ctx, state)
{
    mpn_ctx_t R;

    mpn_ctx_init(R, UWORD(0x0003f00000000001));

    {
        fmpz * a, * c, * d;
        ulong an, zn, zl, zh, sz, i, reps, abits;

        for (reps = 0; reps < 1000 * flint_test_multiplier(); reps++)
        {
            flint_set_num_threads(1 + n_randint(state, 10));

            abits = 1 + n_randint(state, 250);

            an = 1 + n_randint(state, 5000);
            an = 1 + n_randint(state, 1 + an);
            zn = an + an - 1;
            zl = n_randint(state, zn+10);
            zh = n_randint(state, zn+20);

            sz = FLINT_MAX(zl, zh);
            sz = FLINT_MAX(sz, zn);

            a = _fmpz_vec_init(an);
            c = _fmpz_vec_init(sz);
            d = _fmpz_vec_init(sz);

            if (n_randint(state, 2))
                _fmpz_vec_randtest_unsigned(a, state, an, abits);
            else
                _fmpz_vec_randtest(a, state, an, abits);

            if (_fmpz_poly_mul_mid_mpn_ctx(d, zl, zh, a, an, a, an, R))
            {
                _fmpz_poly_mul_KS(c, a, an, a, an);

                for (i = zl; i < zh; i++)
                {
                    if (!fmpz_equal(c + i, d + i-zl))
                    {
                        flint_printf("(squaring) mulmid error at index %wu\n", i);
                        flint_printf("abits=%wu\n", abits);
                        flint_printf("zl=%wu, zh=%wu, an=%wu\n", zl, zh, an);
                        flint_abort();
                    }
                }
            }

            _fmpz_vec_clear(a, an);
            _fmpz_vec_clear(c, sz);
            _fmpz_vec_clear(d, sz);
        }
    }

    {
        fmpz * a, * b, * c, * d;
        ulong an, bn, zn, zl, zh, sz, i, reps, abits, bbits;

        for (reps = 0; reps < 1000 * flint_test_multiplier(); reps++)
        {
            flint_set_num_threads(1 + n_randint(state, 10));

            abits = 1 + n_randint(state, 250);
            bbits = 1 + n_randint(state, 250);

            an = 1 + n_randint(state, 5000);
            an = 1 + n_randint(state, 1 + an);
            bn = 1 + n_randint(state, an);
            zn = an + bn - 1;
            zl = n_randint(state, zn+10);
            zh = n_randint(state, zn+20);

            sz = FLINT_MAX(zl, zh);
            sz = FLINT_MAX(sz, zn);

            a = _fmpz_vec_init(an);
            b = _fmpz_vec_init(bn);
            c = _fmpz_vec_init(sz);
            d = _fmpz_vec_init(sz);

            if (n_randint(state, 2))
                _fmpz_vec_randtest_unsigned(a, state, an, abits);
            else
                _fmpz_vec_randtest(a, state, an, abits);

            if (n_randint(state, 2))
                _fmpz_vec_randtest_unsigned(b, state, bn, bbits);
            else
                _fmpz_vec_randtest(b, state, bn, bbits);

            _fmpz_vec_randtest(d, state, sz, bbits);

            if (_fmpz_poly_mul_mid_mpn_ctx(d, zl, zh, a, an, b, bn, R))
            {
                _fmpz_poly_mul_KS(c, a, an, b, bn);

                for (i = zl; i < zh; i++)
                {
                    if (!fmpz_equal(c + i, d + i-zl))
                    {
                        flint_printf("mulmid error at index %wu\n", i);
                        flint_printf("abits=%wu, bbits=%wu\n", abits, bbits);
                        flint_printf("zl=%wu, zh=%wu, an=%wu, bn=%wu\n", zl, zh, an, bn);
                        flint_abort();
                    }
                }
            }

            _fmpz_vec_clear(a, an);
            _fmpz_vec_clear(b, bn);
            _fmpz_vec_clear(c, sz);
            _fmpz_vec_clear(d, sz);
        }
    }

    mpn_ctx_clear(R);

    TEST_FUNCTION_END(state);
}
