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
#define NB_TESTS 500

TEST_FUNCTION_START(n_fft_itft, state)
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
            ulong iolen, iolen4;
            if (n_randint(state, 4) == 0)
                iolen = 1 + n_randint(state, 12);  /* small, including <= 2 */
            else
                iolen = 1 + n_randint(state, len);  /* [1, len] */
            ulong plen = n_fft_itft_prepare(&iolen4, iolen, F);
            /* p below needs len entries for the full-length DFT */
            const ulong plen2 = FLINT_MAX(plen, len);

            // choose random poly of degree < iolen4
            nmod_poly_t pol;
            nmod_poly_init(pol, mod.n);
            nmod_poly_randtest(pol, state, iolen4);

            // find evals via full-length DFT
            // (the first iolen4 evaluation points at depth are the same as those
            // at any smaller depth >= ceil(log_2(iolen4)), by the bit-reversed order)
            nn_ptr p = _nmod_vec_init(plen2);
            _nmod_vec_set(p, pol->coeffs, pol->length);
            _nmod_vec_zero(p + pol->length, plen2 - pol->length);
            n_fft_dft(p, depth, F);

            /* make sure entries beyond iolen4 do not matter */
            _nmod_vec_randtest(p+iolen4, state, plen2 - iolen4, mod);

            /* apply itft */
            n_fft_itft(p, iolen4, F);

            int res = _nmod_vec_equal(p, pol->coeffs, pol->length) && _nmod_vec_is_zero(p + pol->length, iolen4 - pol->length);

            if (!res)
            {
                _nmod_vec_print(p, iolen4, mod);
                _nmod_vec_print(pol->coeffs, pol->length, mod);
                TEST_FUNCTION_FAIL(
                    "prime = %wu\n"
                    "iolen4 = %wu\n"
                    "pol->length = %wu\n"
                    "root of unity = %wu\n"
                    "max_depth = %wu\n"
                    "depth = %wu\n"
                    "failed equality test\n",
                    prime, iolen4, pol->length, F->tab_w2[2*(max_depth-2)], max_depth, depth);
            }

            _nmod_vec_clear(p);
            nmod_poly_clear(pol);
        }

        n_fft_ctx_clear(F);
    }

    TEST_FUNCTION_END(state);
}

#undef MAX_EVAL_DEPTH
#undef NB_TESTS
