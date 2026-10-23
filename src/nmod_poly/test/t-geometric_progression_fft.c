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

/* Geometric progressions built with stored transforms (bit 3 of function,
   or _nmod_geometric_progression_init_fft with any bounds) give the same
   results as without, for evaluation, interpolation and extrapolation with
   any lengths up to len, and can be shared by threads. */

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
        slong len, ev_fft, int_fft, ext_fft, k, i;
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
#elif FLINT_BITS == 64
                p = UWORD(1108307720798209);
#else
                p = UWORD(2013265921);
#endif
                break;
            case 1:  /* primes transformed directly by fft_small */
                p = (n_randint(state, 2)) ? UWORD(998244353) : UWORD(2013265921);
                break;
            case 2:  /* a composite modulus, for evaluation only */
                p = 2 * len + n_randint(state, n_randint(state, 2) ? 1000 : UWORD(1) << (FLINT_BITS - 2));
                prime = n_is_prime(p);
                break;
            default:
                do
                    p = n_randprime(state, 12 + n_randint(state, FLINT_BITS - 11), 0);
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

        /* stored transforms through bit 3, or through init_fft with bounds
           for each functionality, 0 (no transform) included (possibly
           replacing the transforms stored through bit 3) */
        ev_fft = n_randint(state, 4) ? 1 + n_randint(state, len) : 0;
        int_fft = n_randint(state, 4) ? 1 + n_randint(state, len) : 0;
        ext_fft = n_randint(state, 4) ? 1 + n_randint(state, len) : 0;
        switch (n_randint(state, 3))
        {
            case 0:
                _nmod_geometric_progression_init_function(G1, r, len, mod, function | 8);
                ev_fft = int_fft = ext_fft = len;
                break;
            case 1:
                _nmod_geometric_progression_init_function(G1, r, len, mod, function);
                _nmod_geometric_progression_init_fft(G1, ev_fft, int_fft, ext_fft);
                break;
            default:
                _nmod_geometric_progression_init_function(G1, r, len, mod, function | 8);
                _nmod_geometric_progression_init_fft(G1, ev_fft, int_fft, ext_fft);
        }

        if ((ev_fft == 0 && G1->ev_f_fft != NULL) || (int_fft == 0 && G1->int_f_fft != NULL)
                || (ext_fft == 0 && (G1->ext_ff_fft != NULL || G1->ext_fb_fft != NULL)))
            TEST_FUNCTION_FAIL("transform stored with bound 0: len = %wd\n", len);

        /* one time in two, new bounds for some functionalities, the
           transforms of the others being kept (negative bound) */
        if (n_randint(state, 2))
        {
            const void * ev_old = G1->ev_f_fft;
            const void * int_old = G1->int_f_fft;
            const void * ext_old = G1->ext_ff_fft;
            slong ev_new = n_randint(state, 2) ? -1 - (slong) n_randint(state, 3) : (slong) n_randint(state, len + 1);
            slong int_new = n_randint(state, 2) ? -1 : (slong) n_randint(state, len + 1);
            slong ext_new = n_randint(state, 2) ? -1 : (slong) n_randint(state, len + 1);

            _nmod_geometric_progression_init_fft(G1, ev_new, int_new, ext_new);

            if ((ev_new < 0 && (const void *) G1->ev_f_fft != ev_old)
                    || (int_new < 0 && (const void *) G1->int_f_fft != int_old)
                    || (ext_new < 0 && (const void *) G1->ext_ff_fft != ext_old))
                TEST_FUNCTION_FAIL("transform not kept with a negative bound: len = %wd\n", len);

            if (ev_new >= 0)
                ev_fft = ev_new;
            if (int_new >= 0)
                int_fft = int_new;
            if (ext_new >= 0)
                ext_fft = ext_new;
        }

        a = _nmod_vec_init(len);
        b = _nmod_vec_init(len);
        c = _nmod_vec_init(len);

        for (k = 0; k < 10; k++)
        {
            /* evaluation, the full shape first, then a polynomial longer
               than ev_fft with ev_fft nonzero coefficients */
            if (function & 1)
            {
                slong ilen = (k <= 1) ? len : (slong) n_randint(state, len + 1);
                slong olen = (k <= 1) ? len : (slong) n_randint(state, len + 1);

                _rand_vec(a, state, ilen, mod);
                if (k == 1)
                    _nmod_vec_zero(a, len - ev_fft);
                _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(b, a, ilen, G0, olen, mod);
                _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(c, a, ilen, G1, olen, mod);

                if (!_nmod_vec_equal(b, c, olen))
                    TEST_FUNCTION_FAIL("evaluate: p = %wu, r = %wu, len = %wd, ilen = %wd, olen = %wd\n",
                                       p, r, len, ilen, olen);
            }

            /* interpolation, from all the points, then from the int_fft
               first ones, then from random numbers of points */
            if (function & 2)
            {
                slong m = (k == 0) ? len : (k == 1 && int_fft > 0) ? int_fft
                                         : 1 + (slong) n_randint(state, len);

                _rand_vec(a, state, m, mod);
                _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(b, a, G0, m, mod);
                _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(c, a, G1, m, mod);

                if (!_nmod_vec_equal(b, c, m))
                    TEST_FUNCTION_FAIL("interpolate: p = %wu, r = %wu, len = %wd, m = %wd\n",
                                       p, r, len, m);
            }

            /* extrapolation, forward and backward: ilen + olen <= len; the
               largest shape first, then from at most ext_fft values */
            if ((function & 4) && len >= 3)
            {
                slong olen = (k == 0) ? len / 2 : 1 + (slong) n_randint(state, len - 2);
                slong ilen = (k == 0) ? len - len / 2 : 2 + (slong) n_randint(state, len - olen - 1);
                if (k == 1 && ext_fft >= 2)
                {
                    ilen = FLINT_MIN(ext_fft, len - 1);
                    olen = 1 + n_randint(state, len - ilen);
                }
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
