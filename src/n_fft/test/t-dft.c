/*
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "nmod_vec.h"
#include "n_fft.h"

#define MAX_EVAL_DEPTH 10

// vector equality up to reduction mod
static inline int nmod_vec_red_equal(nn_srcptr vec1, nn_srcptr vec2, ulong len, nmod_t mod)
{
    for (ulong k = 0; k < len; k++)
    {
        ulong v1;
        ulong v2;
        NMOD_RED(v1, vec1[k], mod);
        NMOD_RED(v2, vec2[k], mod);
        if (v1 != v2)
            return 0;
    }

    return 1;
}

// testing that all elements of "vec" are less than "bound"
static inline int nmod_vec_range(nn_srcptr vec, ulong len, ulong bound)
{
    for (ulong k = 0; k < len; k++)
        if (vec[k] >= bound)
            return 0;

    return 1;
}


TEST_FUNCTION_START(n_fft_dft, state)
{
    int i;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        // take some FFT prime p with max_depth >= 12
        ulong max_depth = 12 + n_randint(state, 10);
        ulong p = 1 + (UWORD(1) << max_depth);
        while (! n_is_prime(p))
            p += (UWORD(1) << max_depth);
        max_depth = flint_ctz(p-1);

        nmod_t mod;
        nmod_init(&mod, p);

        // init FFT root tables
        n_fft_ctx_t F;
        n_fft_ctx_init2(F, MAX_EVAL_DEPTH, p);

        for (ulong depth = 3; depth <= MAX_EVAL_DEPTH; depth++)
        {
            const ulong len = (1UL<<depth);

            // choose random poly of degree < len
            nmod_poly_t pol;
            nmod_poly_init(pol, mod.n);
            nmod_poly_randtest(pol, state, len);

            // naive evals by Horner, in bit reversed order
            nn_ptr evals_br = _nmod_vec_init(len);
            for (ulong k = 0; k < len/2; k++)
            {
                ulong point = F->tab_w[2*k];
                evals_br[2*k] = nmod_poly_evaluate_nmod(pol, point);
                evals_br[2*k+1] = nmod_poly_evaluate_nmod(pol, nmod_neg(point, mod));
            }

            // evals by DFT
            ulong * p = _nmod_vec_init(len);
            _nmod_vec_set(p, pol->coeffs, len);

            n_fft_dft(p, len, depth, F);

            int res = nmod_vec_red_equal(evals_br, p, len, mod);

            if (!res)
                TEST_FUNCTION_FAIL(
                        "prime = %wu\n"
                        "root of unity = %wu\n"
                        "max_depth = %wu\n"
                        "depth = %wu\n"
                        "failed equality test\n",
                        p, F->tab_w2[2*(max_depth-2)], max_depth, depth);

            res = nmod_vec_range(p, len, 4*mod.n);

            if (!res)
                TEST_FUNCTION_FAIL(
                        "prime = %wu\n"
                        "root of unity = %wu\n"
                        "max_depth = %wu\n"
                        "depth = %wu\n"
                        "failed range test\n",
                        p, F->tab_w2[2*(max_depth-2)], max_depth, depth);

            _nmod_vec_clear(p);
            nmod_poly_clear(pol);
            _nmod_vec_clear(evals_br);
        }

        n_fft_ctx_clear(F);
    }

    TEST_FUNCTION_END(state);
}
