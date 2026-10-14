/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#if FLINT_HAVE_FFT_SMALL
#include "fft_small.h"
#endif

/* Geometric progressions built with stored transforms (bit 3 of function)
   give the same results as without, for evaluation, interpolation and
   extrapolation with any lengths up to len, and can be shared by threads. */

/* random vector with some zeros at both ends */
static void
_rand_vec(nn_ptr v, flint_rand_t state, slong len, nmod_t mod)
{
    slong i, lo, hi;

    _nmod_vec_rand(v, state, len, mod);

    if (len > 0 && n_randint(state, 3) == 0)
    {
        lo = n_randint(state, len);
        hi = lo + n_randint(state, len - lo + 1);
        for (i = 0; i < lo; i++)
            v[i] = 0;
        for (i = hi; i < len; i++)
            v[i] = 0;
    }
}

typedef struct
{
    const nmod_geometric_progression_struct * G;
    nn_srcptr in;
    nn_ptr out;
    slong len;
}
work_t;

static void
_eval_worker(slong i, void * args)
{
    work_t * w = (work_t *) args + i;
    _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(w->out, w->in, w->len,
                                                        w->G, w->len, w->G->mod);
}

#define NTHREADS 3

TEST_FUNCTION_START(nmod_geometric_progression_fft, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        nmod_geometric_progression_t G0, G1;
        nmod_t mod;
        ulong p, r, function;
        slong len, k, i;
        int prime = 1;
        nn_ptr a, b, c;

        len = 1 + n_randint(state, n_randint(state, 4) ? 1100 : 3000);
        if (n_randint(state, 4) == 0)
            len = (WORD(1) << n_randint(state, 12)) + n_randint(state, 3) - 1;
        len = FLINT_MAX(len, 1);

        switch (n_randint(state, 5))
        {
            case 0:  /* one of the primes of the fft_small context */
#if FLINT_HAVE_FFT_SMALL
                p = get_default_mpn_ctx()->ffts[n_randint(state, 4)].mod.n;
#else
                p = UWORD(1108307720798209);
#endif
                break;
            case 1:  /* primes transformed directly by fft_small */
                p = (n_randint(state, 2)) ? UWORD(998244353) : UWORD(2013265921);
                break;
            case 2:  /* a composite modulus, for evaluation only */
                p = 2 * len + n_randint(state, n_randint(state, 2) ? 1000 : UWORD(1) << 62);
                prime = n_is_prime(p);
                break;
            default:
                do
                    p = n_randprime(state, 12 + n_randint(state, 53), 0);
                while (p <= (ulong) 2 * len + 1);
        }

        nmod_init(&mod, p);

        if (prime)
        {
            r = n_primitive_root_prime(p);
            function = 1 + n_randint(state, 7);
        }
        else
        {
            do
                r = n_randint(state, p);
            while (n_gcd(r, p) != 1);
            function = 1;
        }

        _nmod_geometric_progression_init_function(G0, r, len, mod, function);
        _nmod_geometric_progression_init_function(G1, r, len, mod, function | 8);

        a = _nmod_vec_init(len);
        b = _nmod_vec_init(len);
        c = _nmod_vec_init(len);

        for (k = 0; k < 10; k++)
        {
            /* evaluation, the full shape first */
            if (function & 1)
            {
                slong ilen = (k == 0) ? len : (slong) n_randint(state, len + 1);
                slong olen = (k == 0) ? len : (slong) n_randint(state, len + 1);

                _rand_vec(a, state, ilen, mod);
                _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(b, a, ilen, G0, olen, mod);
                _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(c, a, ilen, G1, olen, mod);

                if (!_nmod_vec_equal(b, c, olen))
                    TEST_FUNCTION_FAIL("evaluate: p = %wu, r = %wu, len = %wd, ilen = %wd, olen = %wd\n",
                                       p, r, len, ilen, olen);
            }

            /* interpolation */
            if (function & 2)
            {
                slong m = (k == 0) ? len : 1 + (slong) n_randint(state, len);

                _rand_vec(a, state, m, mod);
                _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(b, a, G0, m, mod);
                _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(c, a, G1, m, mod);

                if (!_nmod_vec_equal(b, c, m))
                    TEST_FUNCTION_FAIL("interpolate: p = %wu, r = %wu, len = %wd, m = %wd\n",
                                       p, r, len, m);
            }

            /* extrapolation, forward and backward: ilen + olen <= len */
            if ((function & 4) && len >= 3)
            {
                slong olen = (k == 0) ? len / 2 : 1 + (slong) n_randint(state, len - 2);
                slong ilen = (k == 0) ? len - len / 2 : 2 + (slong) n_randint(state, len - olen - 1);
                slong offset;

                if (n_randint(state, 2))
                    offset = ilen + n_randint(state, len - olen - ilen + 1);
                else
                    offset = - olen - (slong) n_randint(state, len - olen - ilen + 1);

                _rand_vec(a, state, ilen, mod);
                nmod_poly_extrapolate_geometric_precomp(b, olen, a, ilen, offset, G0);
                nmod_poly_extrapolate_geometric_precomp(c, olen, a, ilen, offset, G1);

                if (!_nmod_vec_equal(b, c, olen))
                    TEST_FUNCTION_FAIL("extrapolate: p = %wu, r = %wu, len = %wd, ilen = %wd, olen = %wd, offset = %wd\n",
                                       p, r, len, ilen, olen, offset);
            }
        }

        /* several threads evaluating with the same G1 */
        if ((function & 1) && iter % 10 == 0)
        {
            work_t w[NTHREADS];
            nn_ptr in = _nmod_vec_init(NTHREADS * len);
            nn_ptr out = _nmod_vec_init(NTHREADS * len);

            _nmod_vec_rand(in, state, NTHREADS * len, mod);

            for (i = 0; i < NTHREADS; i++)
            {
                w[i].G = G1;
                w[i].in = in + i * len;
                w[i].out = out + i * len;
                w[i].len = len;
            }

            k = flint_get_num_threads();
            flint_set_num_threads(NTHREADS);
            flint_parallel_do(_eval_worker, w, NTHREADS, NTHREADS, FLINT_PARALLEL_STRIDED);
            flint_set_num_threads(k);

            for (i = 0; i < NTHREADS; i++)
            {
                _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(b, in + i * len, len, G0, len, mod);
                if (!_nmod_vec_equal(b, out + i * len, len))
                    TEST_FUNCTION_FAIL("threads: p = %wu, r = %wu, len = %wd\n", p, r, len);
            }

            _nmod_vec_clear(in);
            _nmod_vec_clear(out);
        }

        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
        _nmod_vec_clear(c);
        nmod_geometric_progression_clear(G0);
        nmod_geometric_progression_clear(G1);
    }

    TEST_FUNCTION_END(state);
}

#undef NTHREADS
