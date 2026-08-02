/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fft_small.h"

/* Accumulate K integer products in transform space and compare the
   exported limbs against mpn_mul. Inflating len_bound drives the plan
   across np = 4..8 and onto varying chunk sizes. */
TEST_FUNCTION_START(fft_small_op_mpn, state)
{
    mpn_ctx_t R;

    mpn_ctx_init(R, UWORD(0x0003f00000000001));

    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        int large = n_randint(state, 8) == 0;
        ulong K = 1 + n_randint(state, large ? 2 : 4);
        /* the occasional large iterations cross the internal threading
           thresholds of fft_small_fft_mpn and fft_small_export_mpn,
           including the overhang stitching between crt ranges */
        ulong an = large ? 1100 + n_randint(state, 2000)
                         : 1 + n_randint(state, 150);
        ulong bn = large ? 1100 + n_randint(state, 2000)
                         : 1 + n_randint(state, 150);
        ulong zn;

        /* equal lengths every few iterations so the square op actually
           runs (independent draws almost never coincide); this must
           precede every use of bn: the sizes, the plan and the buffers */
        if (iter % 3 == 0)
            bn = an;
        zn = an + bn;
        /* weighted toward large inflations, where the transform-cost
           score starts preferring more primes with larger chunks */
        ulong inflate = n_randint(state, 2) ? n_randint(state, 60)
                                            : 45 + n_randint(state, 15);
        ulong len_bound;
        ulong i, k;
        int use_sqr, use_sub;
        fft_small_plan_t P;
        fft_small_op_t A, B, Z, T;
        ulong *a, *b, *z, *ref, *t;

        flint_set_num_threads(1 + n_randint(state, 8));

        /* inflate the accumulation bound to sweep prime counts and
           chunk sizes */
        if (inflate >= FLINT_BITS - FLINT_BIT_COUNT(K))
            len_bound = UWORD_MAX;   /* certainly fails: tests the 0 return */
        else
            len_bound = K << inflate;

        if (!fft_small_plan_init_mpn(P, R, an, bn, len_bound))
        {
            continue;
        }

        a = FLINT_ARRAY_ALLOC(K*an, ulong);
        b = FLINT_ARRAY_ALLOC(K*bn, ulong);
        z = FLINT_ARRAY_ALLOC(zn, ulong);
        ref = FLINT_ARRAY_ALLOC(zn, ulong);
        t = FLINT_ARRAY_ALLOC(zn, ulong);

        use_sqr = (an == bn) && n_randint(state, 2) == 0;
        use_sub = (K >= 2) && n_randint(state, 3) == 0;

        for (k = 0; k < K; k++)
        {
            flint_mpn_rrandom(a + k*an, state, an);
            /* keep the K-fold sum below 2^(FLINT_BITS*zn) so that the
               reference accumulation cannot carry out of zn limbs,
               except in the untruncated K = 1 case */
            if (K > 1)
                a[k*an + an - 1] = 0;

            if (use_sqr)
                flint_mpn_copyi(b + k*bn, a + k*an, bn);
            else
                flint_mpn_rrandom(b + k*bn, state, bn);
        }

        fft_small_op_init(A, P);
        fft_small_op_init(B, P);
        fft_small_op_init(Z, P);
        fft_small_op_init(T, P);

        for (k = 0; k < K; k++)
        {
            fft_small_fft_mpn(A, a + k*an, an, P);

            if (use_sqr)
            {
                if (k == 0)
                    fft_small_op_sqr(Z, A, P);
                else
                {
                    fft_small_op_sqr(T, A, P);
                    fft_small_op_add(Z, Z, T, P);
                }
            }
            else
            {
                fft_small_fft_mpn(B, b + k*bn, bn, P);

                if (k == 0)
                    fft_small_op_mul(Z, A, B, P);
                else
                    fft_small_op_addmul(Z, A, B, P);
            }
        }

        if (use_sub)
        {
            /* add and subtract the same product: the true value stays
               nonnegative as the unsigned reconstruction requires */
            fft_small_fft_mpn(A, a + (K-1)*an, an, P);
            fft_small_fft_mpn(B, b + (K-1)*bn, bn, P);
            fft_small_op_addmul(Z, A, B, P);
            fft_small_op_submul(Z, A, B, P);
        }

        fft_small_ifft(Z, P);
        fft_small_export_mpn(z, zn, Z, P);

        /* reference */
        flint_mpn_zero(ref, zn);
        for (k = 0; k < K; k++)
        {
            if (an >= bn)
                flint_mpn_mul(t, a + k*an, an, b + k*bn, bn);
            else
                flint_mpn_mul(t, b + k*bn, bn, a + k*an, an);

            i = mpn_add_n(ref, ref, t, zn);
            FLINT_ASSERT(i == 0);
        }

        for (i = 0; i < zn; i++)
        {
            if (z[i] != ref[i])
            {
                TEST_FUNCTION_FAIL(
                        "limb i = %wu, zn = %wu, an = %wu, bn = %wu\n"
                        "K = %wu, inflate = %wu\n"
                        "np = %wu, bits = %wu, depth = %wu, ztrunc = %wu\n"
                        "use_sqr = %d, use_sub = %d\n"
                        "got %wx, expected %wx\n",
                        i, zn, an, bn, K, inflate,
                        P->np, P->bits, P->depth, P->ztrunc,
                        use_sqr, use_sub,
                        z[i], ref[i]);
            }
        }

        fft_small_op_clear(A);
        fft_small_op_clear(B);
        fft_small_op_clear(Z);
        fft_small_op_clear(T);
        fft_small_plan_clear(P);

        flint_free(a);
        flint_free(b);
        flint_free(z);
        flint_free(ref);
        flint_free(t);
    }

    mpn_ctx_clear(R);

    TEST_FUNCTION_END(state);
}
