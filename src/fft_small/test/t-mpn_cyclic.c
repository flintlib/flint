/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fft_small.h"

/* r := v mod (B^nn - 1), canonical in [0, B^nn - 1) */
static void
_fold_bnm1(nn_ptr r, nn_srcptr v, ulong vn, ulong nn)
{
    ulong cy, m, i;

    m = FLINT_MIN(vn, nn);
    flint_mpn_copyi(r, v, m);
    flint_mpn_zero(r + m, nn - m);
    v += m;
    vn -= m;

    while (vn > 0)
    {
        m = FLINT_MIN(vn, nn);
        cy = mpn_add(r, r, nn, v, m);
        while (cy)
            cy = mpn_add_1(r, r, nn, cy);   /* B^nn ≡ 1 */
        v += m;
        vn -= m;
    }

    for (i = 0; i < nn; i++)
        if (r[i] != UWORD_MAX)
            return;
    flint_mpn_zero(r, nn);                  /* B^nn - 1 ≡ 0 */
}

TEST_FUNCTION_START(fft_small_mpn_cyclic, state)
{
    mpn_ctx_t R;

    mpn_ctx_init(R, UWORD(0x0003f00000000001));

    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        int large = n_randint(state, 10) == 0;
        ulong K = 1 + n_randint(state, 3);
        ulong nn_req = large ? 2000 + n_randint(state, 6000)
                             : 100 + n_randint(state, 1500);
        ulong nn = nn_req;
        ulong len_bound = K;
        ulong an, bn, zn, i, k;
        int use_sqr;
        fft_small_plan_t P;
        fft_small_op_t A, B, Z, T;
        ulong *a, *b, *v, *z, *ref, *w;

        flint_set_num_threads(1 + n_randint(state, 8));

        if (n_randint(state, 30) == 0)
            len_bound = UWORD_MAX;   /* certainly fails: tests the 0 return */

        if (!fft_small_plan_init_mpn_cyclic(P, R, &nn, len_bound))
        {
            if (len_bound != UWORD_MAX)
                TEST_FUNCTION_FAIL("no plan for nn = %wu, len_bound = %wu\n",
                                   nn_req, len_bound);
            continue;
        }

        if (nn < nn_req || P->bits * n_pow2(P->depth) != FLINT_BITS * nn)
            TEST_FUNCTION_FAIL("bad geometry: nn_req = %wu, nn = %wu, "
                               "bits = %wu, depth = %wu\n",
                               nn_req, nn, P->bits, P->depth);

        /* exact representative bound from the plan.c comment */
        zn = nn + (2*P->bits + P->depth + FLINT_BIT_COUNT(len_bound))
                    / FLINT_BITS + 1;

        /* operand lengths: exactly nn most of the time, so the exact
           chunk fit and the wraparound actually engage */
        an = (n_randint(state, 4) != 0) ? nn : 1 + n_randint(state, nn);
        bn = (n_randint(state, 4) != 0) ? nn : 1 + n_randint(state, nn);

        a = FLINT_ARRAY_ALLOC(K*nn, ulong);
        b = FLINT_ARRAY_ALLOC(K*nn, ulong);
        v = FLINT_ARRAY_ALLOC(2*nn, ulong);
        z = FLINT_ARRAY_ALLOC(zn, ulong);
        ref = FLINT_ARRAY_ALLOC(nn, ulong);
        w = FLINT_ARRAY_ALLOC(nn, ulong);

        use_sqr = (an == bn) && n_randint(state, 3) == 0;

        for (k = 0; k < K; k++)
        {
            flint_mpn_rrandom(a + k*nn, state, an);
            if (n_randint(state, 20) == 0)          /* B^nn - 1 ≡ 0 */
                for (i = 0; i < an; i++)
                    a[k*nn + i] = UWORD_MAX;
            if (use_sqr)
                flint_mpn_copyi(b + k*nn, a + k*nn, bn);
            else
                flint_mpn_rrandom(b + k*nn, state, bn);
        }

        fft_small_op_init(A, P);
        fft_small_op_init(B, P);
        fft_small_op_init(Z, P);
        fft_small_op_init(T, P);

        /* the second and later products reuse B's transform from the
           k = 0 round when the operands coincide, exercising the cached
           PRIMAL reuse the plan interface promises */
        for (k = 0; k < K; k++)
        {
            fft_small_fft_mpn(A, a + k*nn, an, P);

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
                if (k == 0 || n_randint(state, 2) == 0)
                    fft_small_fft_mpn(B, b + k*nn, bn, P);
                else
                    flint_mpn_copyi(b + k*nn, b + (k-1)*nn, nn);

                if (k == 0)
                    fft_small_op_mul(Z, A, B, P);
                else
                    fft_small_op_addmul(Z, A, B, P);
            }
        }

        fft_small_ifft(Z, P);
        fft_small_export_mpn(z, zn, Z, P);
        _fold_bnm1(w, z, zn, nn);

        /* reference: sum of full products, folded */
        flint_mpn_zero(ref, nn);
        for (k = 0; k < K; k++)
        {
            ulong cy;

            if (an >= bn)
                mpn_mul(v, a + k*nn, an, b + k*nn, bn);
            else
                mpn_mul(v, b + k*nn, bn, a + k*nn, an);
            _fold_bnm1(z, v, an + bn, nn);   /* reuse z as scratch */
            cy = mpn_add_n(ref, ref, z, nn);
            while (cy)
                cy = mpn_add_1(ref, ref, nn, cy);
        }
        for (i = 0; i < nn; i++)
            if (ref[i] != UWORD_MAX)
                break;
        if (i == nn)
            flint_mpn_zero(ref, nn);

        if (mpn_cmp(w, ref, nn) != 0)
            TEST_FUNCTION_FAIL("nn_req = %wu, nn = %wu, an = %wu, bn = %wu, "
                               "K = %wu, np = %wu, bits = %wu, depth = %wu, "
                               "sqr = %d\n",
                               nn_req, nn, an, bn, K,
                               P->np, P->bits, P->depth, use_sqr);

        fft_small_op_clear(A);
        fft_small_op_clear(B);
        fft_small_op_clear(Z);
        fft_small_op_clear(T);
        fft_small_plan_clear(P);
        flint_free(a);
        flint_free(b);
        flint_free(v);
        flint_free(z);
        flint_free(ref);
        flint_free(w);
    }

    mpn_ctx_clear(R);

    TEST_FUNCTION_END(state);
}
