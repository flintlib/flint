/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "fft_small.h"

TEST_FUNCTION_START(_nmod_poly_divrem_mpn_ctx, state)
{
    flint_bitcnt_t nbits;
    mpn_ctx_t R;
    nmod_t mod;

    mpn_ctx_init(R, UWORD(0x0003f00000000001));

    {
        ulong * a, * b, * q1, * q2, * q3, * r1, * r2, * r3;
        ulong an, bn, qn, i, reps;
        nmod_poly_divrem_precomp_struct M[1];

        for (reps = 0; reps < 1000 * flint_test_multiplier(); reps++)
        {
            flint_set_num_threads(1 + n_randint(state, 10));

            /* 2 <= nbits <= FLINT_BITS */
            nbits = 2 + n_randint(state, FLINT_BITS - 1);
            nmod_init(&mod, n_randbits(state, nbits));

            bn = 2 + n_randint(state, 5000);
            qn = 1 + n_randint(state, 5000);
            an = bn + qn - 1;

            a = FLINT_ARRAY_ALLOC(an, ulong);
            b = FLINT_ARRAY_ALLOC(bn, ulong);
            q1 = FLINT_ARRAY_ALLOC(qn, ulong);
            q2 = FLINT_ARRAY_ALLOC(qn, ulong);
            q3 = FLINT_ARRAY_ALLOC(qn, ulong);
            r1 = FLINT_ARRAY_ALLOC(bn, ulong);
            r2 = FLINT_ARRAY_ALLOC(bn, ulong);
            r3 = FLINT_ARRAY_ALLOC(bn, ulong);

            for (i = 0; i < an; i++)
                a[i] = n_randint(state, mod.n);

            for (i = 0; i < bn; i++)
                b[i] = n_randint(state, mod.n);

            while (n_gcd(b[bn-1], mod.n) != 1)
                b[bn-1] = n_randint(state, mod.n);

            _nmod_poly_divrem(q1, r1, a, an, b, bn, mod);

            _nmod_poly_divrem_mpn_ctx(q2, r2, a, an, b, bn, mod, R);

            ulong prec = qn + n_randint(state, 200);
            _nmod_poly_divrem_precomp_init(M, b, bn, prec, mod, R);
            _nmod_poly_divrem_precomp(q3, r3, a, an, M, mod, R);
            _nmod_poly_divrem_precomp_clear(M);

            for (i = qn; i > 0; i--)
            {
                if (q1[i-1] != q2[i-1] || q1[i-1] != q3[i-1])
                {
                    flint_printf("quotient error at index %wu\n", i-1);
                    flint_printf("qn=%wu, an=%wu, bn=%wu\n", qn, an, bn);
                    flint_printf("mod: %wu\n", mod.n);
                    flint_abort();
                }
            }

            for (i = bn-1; i > 0; i--)
            {
                if (r1[i-1] != r2[i-1] || r1[i-1] != r3[i-1])
                {
                    flint_printf("remainder error at index %wu\n", i-1);
                    flint_printf("r1[i]=%wu, r2[i]=%wu, bn=%wu\n", r1[i-1], r2[i-1]);
                    flint_printf("qn=%wu, an=%wu, bn=%wu\n", qn, an, bn);
                    flint_printf("mod: %wu\n", mod.n);
                    flint_abort();
                }
            }

            flint_free(a);
            flint_free(b);
            flint_free(q1);
            flint_free(q2);
            flint_free(r1);
            flint_free(r2);
        }
    }

    mpn_ctx_clear(R);

    TEST_FUNCTION_END(state);
}
