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

#define MAX_EVAL_DEPTH 11

TEST_FUNCTION_START(n_fft_tft, state)
{
    int i;
    flint_rand_set_seed(state, rand(), rand()+12);
    ulong seed1, seed2;
    flint_rand_get_seed(&seed1, &seed2, state);
    flint_rand_set_seed(state, 846930886L, 1804289395L);
    flint_printf("seeds: %wu, %wu\n", seed1, seed2);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        // take some FFT prime p with max_depth >= 10
        ulong max_depth, prime;

        // half of tests == fixed large prime, close to limit
        // 62 bits: prime = 4611686018427322369 == 2**62 - 2**16 + 1
        // 30 bits: prime = 1073479681 == 2**30 - 2**18 + 1
        if (i > 100)
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

        // retrieve roots, used later for multipoint evaluation
        nn_ptr roots = flint_malloc((UWORD(1) << MAX_EVAL_DEPTH) * sizeof(ulong));
        for (ulong k = 0; k < (UWORD(1) << (MAX_EVAL_DEPTH-1)); k++)
        {
            roots[2*k] = F->tab_w[2*k];
            roots[2*k+1] = prime - F->tab_w[2*k];  // < prime since F->tab_w[2*k] != 0
        }

        for (ulong depth = 4; depth <= MAX_EVAL_DEPTH; depth++)
        {
            const ulong len = (UWORD(1) << depth);
            /* FOR V1: */
            /* const ulong ilen = 8 + n_randint(state, len); */
            /* const ulong olen = 8 * ((8 + n_randint(state, len)) / 8); */
            /* FOR V2: */
            const ulong ilen = len;
            ulong olen = len/2 + 4 * (n_randint(state, len/2) / 4);  /* (len/2, len) */
            if (olen == len/2) olen += 4;

            flint_printf("---\n"
                    "prime = %wu\n"
                    "ilen = %wu\n"
                    "olen = %wu\n"
                    "root of unity = %wu\n"
                    "max_depth = %wu\n"
                    "depth = %wu\n",
                    prime, ilen, olen, F->tab_w2[2*(max_depth-2)], max_depth, depth);

            // choose random poly of degree < ilen
            nmod_poly_t pol;
            nmod_poly_init(pol, mod.n);
            nmod_poly_randtest(pol, state, ilen);
            // copy it for DFT
            nn_ptr p = _nmod_vec_init(FLINT_MAX(len, ilen));
            _nmod_vec_set(p, pol->coeffs, ilen);
            /* TODO need to fill ilen..len with zeros? */

            // evals via general multipoint evaluation
            nn_ptr evals_br = _nmod_vec_init(olen);
            nmod_poly_evaluate_nmod_vec(evals_br, pol, roots, olen);

            // evals by TFT
            n_fft_args_t Fargs;
            n_fft_set_args(Fargs, F->mod, F->tab_w);
            /* tft_node_lazy_4_4_v1(p, ilen, olen, depth, 0, Fargs); */
            tft_node_lazy_4_4_v2_olen(p, olen, depth, 0, Fargs);

            for (ulong k = 0; k < olen; k++)
            {
                if (p[k] >= 2*prime)
                    p[k] -= 2*prime;
                if (p[k] >= prime)
                    p[k] -= prime;
            }
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

        flint_free(roots);
        n_fft_ctx_clear(F);
    }

    TEST_FUNCTION_END(state);
}

#undef MAX_EVAL_DEPTH
