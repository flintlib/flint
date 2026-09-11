/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "nmod_vec.h"
#include "n_fft.h"
#include "n_fft/impl.h"

#define MAX_EVAL_DEPTH 11
#define NB_TESTS 1000

TEST_FUNCTION_START(n_fft_tft, state)
{
    int i;

    for (i = 0; i < NB_TESTS * flint_test_multiplier(); i++)
    {
        // take some FFT prime p with max_depth >= MAX_EVAL_DEPTH
        ulong max_depth, prime;

        // half of tests == fixed large prime, close to limit
        // 62 bits: prime = 4611686018427322369 == 2**62 - 2**16 + 1
        // 30 bits: prime = 1073479681 == 2**30 - 2**18 + 1
        if (2*i > NB_TESTS * flint_test_multiplier())
#if FLINT_BITS == 64
            prime = UWORD(4611686018427322369);
#else // FLINT_BITS == 32
            prime = UWORD(1073479681);
#endif
        else
        {
            max_depth = MAX_EVAL_DEPTH + n_randint(state, 6);
            prime = 1 + (UWORD(1) << max_depth);
            while (! n_is_prime(prime))
                prime += (UWORD(1) << max_depth);
        }
        max_depth = flint_ctz(prime-1);

        nmod_t mod;
        nmod_init(&mod, prime);

        // init FFT root tables
        n_fft_ctx_t F;
        n_fft_ctx_init2(F, MAX_EVAL_DEPTH, prime);

        for (ulong depth = 4; depth <= MAX_EVAL_DEPTH; depth++)
        {
            const ulong len = (UWORD(1) << depth);
            ulong plen, ilen, ilen4, olen, olen4;
            ilen = 5 + n_randint(state, 4*len);  /* [5, 4*len] */
            /* FIXME handle smaller ilen */
            olen = 1+len/2 + n_randint(state, len/2);  /* (len/2, len] */
            plen = n_fft_tft_prepare(&ilen4, &olen4, ilen, olen, F);

            /* flint_printf("---\n" */
            /*         "prime = %wu\n" */
            /*         "ilen = %wu\n" */
            /*         "olen = %wu\n" */
            /*         "ilen4 = %wu\n" */
            /*         "olen4 = %wu\n" */
            /*         "plen = %wu\n" */
            /*         "root of unity = %wu\n" */
            /*         "max_depth = %wu\n" */
            /*         "depth = %wu\n", */
            /*         prime, ilen, olen, */
            /*         ilen4, olen4, plen, */
            /*         F->tab_w2[2*(max_depth-2)], max_depth, depth); */

            // choose random poly of degree < ilen
            nmod_poly_t pol;
            nmod_poly_init(pol, mod.n);
            nmod_poly_randtest(pol, state, ilen);

            // find evals via full-length DFT
            nn_ptr evals_br = _nmod_vec_init(plen);
            _nmod_vec_set(evals_br, pol->coeffs, pol->length);
            _nmod_vec_zero(evals_br + pol->length, plen - pol->length);
            if (pol->length > (slong)len)
                _nmod_poly_divrem_circulant1(evals_br, plen, len, prime);
            n_fft_dft(evals_br, depth, F);

            // find evals by TFT
            nn_ptr p = _nmod_vec_init(plen);
            _nmod_vec_set(p, pol->coeffs, pol->length);
            _nmod_vec_zero(p + pol->length, plen - pol->length);
            n_fft_args_t Fargs;
            n_fft_set_args(Fargs, F->mod, F->tab_w);
            n_fft_tft(p, ilen4, olen4, F);

            /* tft_lazy_1_4(p, ilen4, olen4, Fargs); */
            /* for (ulong k = 0; k < olen; k++) */
            /* { */
            /*     if (p[k] >= 2*prime) */
            /*         p[k] -= 2*prime; */
            /*     if (p[k] >= prime) */
            /*         p[k] -= prime; */
            /* } */

            int res = _nmod_vec_equal(evals_br, p, olen);

            if (!res)
            {
                _nmod_vec_print(evals_br, olen, mod);
                _nmod_vec_print(p, olen, mod);
                TEST_FUNCTION_FAIL(
                    "prime = %wu\n"
                    "ilen = %wu\n"
                    "olen = %wu\n"
                    "root of unity = %wu\n"
                    "max_depth = %wu\n"
                    "depth = %wu\n"
                    "failed equality test\n",
                    prime, ilen, olen, F->tab_w2[2*(max_depth-2)], max_depth, depth);
            }

            _nmod_vec_clear(p);
            nmod_poly_clear(pol);
            _nmod_vec_clear(evals_br);
        }

        n_fft_ctx_clear(F);
    }

    TEST_FUNCTION_END(state);
}

#undef MAX_EVAL_DEPTH
#undef NB_TESTS
