/*
    Copyright (C) 2024 Vincent Neiger

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
#include "nmod_vec.h"
#include "n_fft.h"

#define MAX_EVAL_DEPTH 9

/** computes the weighted power sums
 *      q == [PowerSum(p, w**j) for 0 <= j < len]
 * where PowerSum(p, w**j) == sum(p[i] * w[i]**j for 0 <= i < len)
 * and where roots == [w[i] for 0 <= i < len]
 */
static void t_dft_t_weighted_power_sums(nn_ptr q, nn_srcptr p, nn_ptr roots, ulong len, nmod_t mod)
{
    // initially w**0 == [1,..,1]:
    nn_ptr w_pow_j = _nmod_vec_init(len);
    for (ulong i = 0; i < len; i++)
        w_pow_j[i] = 1;

    for (ulong j = 0; j < len; j++)
    {
        // at this stage, w_pow_j holds [w[i]**j for 0 <= i < len]
        q[j] = 0;
        for (ulong i = 0; i < len; i++)
        {
            q[j] = nmod_add(q[j], 
                            nmod_mul(p[i], w_pow_j[i], mod),
                            mod);
            w_pow_j[i] = nmod_mul(w_pow_j[i], roots[i], mod);
        }
    }
    _nmod_vec_clear(w_pow_j);
}

TEST_FUNCTION_START(n_fft_dft_t, state)
{
    int i;

    for (i = 0; i < 200 * flint_test_multiplier(); i++)
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

        // retrieve roots
        nn_ptr roots = flint_malloc((UWORD(1) << MAX_EVAL_DEPTH) * sizeof(ulong));
        for (ulong k = 0; k < (UWORD(1) << (MAX_EVAL_DEPTH-1)); k++)
        {
            roots[2*k] = F->tab_w[2*k];
            roots[2*k+1] = prime - F->tab_w[2*k];  // < prime since F->tab_w[2*k] != 0
        }

        for (ulong depth = 1; depth <= MAX_EVAL_DEPTH; depth++)
        {
            const ulong len = (UWORD(1) << depth);

            // construct random array of length len
            nn_ptr p = _nmod_vec_init(len);
            for (ulong k = 0; k < len; k++)
                p[k] = n_randint(state, prime);
            // copy it before in-place transform
            ulong * q = _nmod_vec_init(len);
            _nmod_vec_set(q, p, len);

            // naive weighted power sums
            t_dft_t_weighted_power_sums(q, p, roots, len, mod);

            // transposed DFT
            n_fft_dft_t(p, depth, F);

            int res = _nmod_vec_equal(p, q, len);

            if (!res)
                TEST_FUNCTION_FAIL(
                        "prime = %wu\n"
                        "root of unity = %wu\n"
                        "max_depth = %wu\n"
                        "depth = %wu\n"
                        "failed equality test\n",
                        prime, F->tab_w2[2*(max_depth-2)], max_depth, depth);

            _nmod_vec_clear(p);
            _nmod_vec_clear(q);
        }

        flint_free(roots);
        n_fft_ctx_clear(F);
    }

    TEST_FUNCTION_END(state);
}

#undef MAX_EVAL_DEPTH
