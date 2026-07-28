/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fft_small.h"

/* Accumulate K products in transform space (as a bilinear application
   like matrix multiplication would) and compare the exported window
   against classical multiplications. Inflating prod_bits beyond the
   2*modbits an ordinary multiplication needs also drives the plan onto
   5 to 8 primes, exercising the template-generated nmod chinese
   remaindering. */
TEST_FUNCTION_START(fft_small_op_nmod, state)
{
    mpn_ctx_t R;

    mpn_ctx_init(R, UWORD(0x0003f00000000001));

    for (slong iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        int large = n_randint(state, 8) == 0;
        ulong K = 1 + n_randint(state, large ? 2 : 5);
        /* the occasional large iterations cross the internal threading
           thresholds of fft_small_fft_nmod, fft_small_ifft and
           fft_small_export_nmod */
        ulong an = large ? 1500 + n_randint(state, 3500)
                         : 1 + n_randint(state, 700);
        ulong bn = large ? 1500 + n_randint(state, 3500)
                         : 1 + n_randint(state, 700);
        ulong zn = an + bn - 1;
        ulong zl = n_randint(state, zn);
        ulong zh = zl + 1 + n_randint(state, zn - zl);
        ulong atrunc = n_round_up(an, BLK_SZ);
        ulong btrunc = n_round_up(bn, BLK_SZ);
        ulong extra_bits = n_randint(state, 300);
        ulong modbits = 2 + n_randint(state, 49);
        ulong i, k;
        int use_sqr, use_sub, use_add;
        fft_small_plan_t P;
        fft_small_op_t A, B, Z, T;
        ulong *a, *b, *z, *ref, *t, *ta;

        flint_set_num_threads(1 + n_randint(state, 8));

        nmod_init(&mod, n_randbits(state, modbits));

        /* z = sum_k a_k * b_k needs len_bound = K*min(an, bn); inflate
           prod_bits to force more primes than 2*modbits would need */
        if (!fft_small_plan_init_nmod(P, R, zl, zh, zn,
                        FLINT_MAX(atrunc, btrunc), K*FLINT_MIN(an, bn),
                        2*modbits + extra_bits, mod, 0))
        {
            continue;   /* bound needs more than MPN_CTX_NCRTS primes */
        }

        a = FLINT_ARRAY_ALLOC(K*an, ulong);
        b = FLINT_ARRAY_ALLOC(K*bn, ulong);
        z = FLINT_ARRAY_ALLOC(zh - zl, ulong);
        ref = FLINT_ARRAY_ALLOC(zn, ulong);
        t = FLINT_ARRAY_ALLOC(zn, ulong);
        ta = FLINT_ARRAY_ALLOC(an, ulong);

        use_sqr = (an == bn) && n_randint(state, 4) == 0;
        /* an extra pair with subtracted product keeps the true
           coefficients nonnegative, which the unsigned nmod export
           requires */
        use_sub = (K >= 2) && n_randint(state, 3) == 0;
        /* sometimes fold a_1 into a_0 with a pointwise add instead of a
           separate product; the sums are unreduced, so this needs a bit
           of bound headroom */
        use_add = !use_sqr && !use_sub && (K >= 2) && extra_bits >= 2 &&
                  n_randint(state, 3) == 0;

        for (k = 0; k < K; k++)
        {
            _nmod_vec_randtest(a + k*an, state, an, mod);
            if (use_sqr)
                flint_mpn_copyi(b + k*bn, a + k*an, bn);
            else
                _nmod_vec_randtest(b + k*bn, state, bn, mod);
        }

        fft_small_op_init(A, P);
        fft_small_op_init(B, P);
        fft_small_op_init(Z, P);
        fft_small_op_init(T, P);

        /* transform space: Z = sum_k A_k B_k (with the use_* variations) */
        for (k = 0; k < K; k++)
        {
            fft_small_fft_nmod(A, a + k*an, an, atrunc, mod, P);

            if (use_add && k == 0)
            {
                /* A = A_0 + A_1 pointwise, then skip k == 1 */
                fft_small_fft_nmod(T, a + 1*an, an, atrunc, mod, P);
                fft_small_op_add(A, A, T, P);
            }

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
                fft_small_fft_nmod(B, b + k*bn, bn, btrunc, mod, P);

                if (k == 0)
                    fft_small_op_mul(Z, A, B, P);
                else
                    fft_small_op_addmul(Z, A, B, P);
            }

            if (use_add && k == 0)
                k++;   /* a_1 was already folded in */
        }

        if (use_sub)
        {
            /* subtract the last product again: Z += A_{K-1} B_{K-1}
               was done above, now Z -= A_{K-1} B_{K-1} */
            fft_small_fft_nmod(A, a + (K-1)*an, an, atrunc, mod, P);
            fft_small_fft_nmod(B, b + (K-1)*bn, bn, btrunc, mod, P);
            fft_small_op_submul(Z, A, B, P);
        }

        fft_small_ifft(Z, P);
        fft_small_export_nmod(z, Z, mod, P);

        /* reference via classical multiplication */
        _nmod_vec_zero(ref, zn);
        for (k = 0; k < K; k++)
        {
            if (use_sub && k == K - 1)
                continue;

            if (use_add && k == 0)
            {
                _nmod_vec_add(ta, a + 0*an, a + 1*an, an, mod);
                if (an >= bn)
                    _nmod_poly_mul(t, ta, an, b + 0*bn, bn, mod);
                else
                    _nmod_poly_mul(t, b + 0*bn, bn, ta, an, mod);
            }
            else if (an >= bn)
                _nmod_poly_mul(t, a + k*an, an, b + k*bn, bn, mod);
            else
                _nmod_poly_mul(t, b + k*bn, bn, a + k*an, an, mod);

            _nmod_vec_add(ref, ref, t, zn, mod);

            if (use_add && k == 0)
                k++;
        }

        for (i = zl; i < zh; i++)
        {
            if (z[i - zl] != ref[i])
            {
                TEST_FUNCTION_FAIL(
                        "i = %wu, zl = %wu, zh = %wu, zn = %wu\n"
                        "an = %wu, bn = %wu, K = %wu\n"
                        "np = %wu, depth = %wu, ztrunc = %wu\n"
                        "modulus = %wu, extra_bits = %wu\n"
                        "use_sqr = %d, use_sub = %d, use_add = %d\n"
                        "got %wu, expected %wu\n",
                        i, zl, zh, zn, an, bn, K,
                        P->np, P->depth, P->ztrunc,
                        mod.n, extra_bits,
                        use_sqr, use_sub, use_add,
                        z[i - zl], ref[i]);
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
        flint_free(ta);
    }

    mpn_ctx_clear(R);

    TEST_FUNCTION_END(state);
}
