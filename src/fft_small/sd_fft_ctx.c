/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fft_small.h"
#include "nmod.h"

void sd_fft_ctx_clear(sd_fft_ctx_t Q)
{
    ulong k;
    flint_aligned_free(Q->w2tab[0]);
    for (k = SD_FFT_CTX_W2TAB_INIT; k < SD_FFT_CTX_W2TAB_SIZE; k++)
        flint_aligned_free(Q->w2tab[k]);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&Q->mutex);
#endif
}

/*
    Return a primitive 2^depth-th root modulo the prime pp.
    Requires depth == valuation(pp - 1, 2).
*/
static ulong
sd_fft_ctx_primitive_2power_root(ulong pp, ulong depth, nmod_t mod)
{
    ulong a = n_quadratic_nonresidue(pp);
    return nmod_pow_ui(a, (pp - 1) >> depth, mod);
}

/*
    Return the primitive 2^(k+1)-th root used to generate w2tab[k].
    Requires depth == valuation(Q->mod.n - 1, 2).
*/
static ulong
sd_fft_ctx_w2tab_root(const sd_fft_ctx_t Q, ulong depth, ulong k)
{
    FLINT_ASSERT(k + 1 <= depth);
    return nmod_pow_ui(Q->primitive_2power_root, UWORD(1) << (depth - k - 1), Q->mod);
}

/*
    Initialize FFT context.
    pp is a prime with at most ~ 50 bits (exactly representable with a `double`)
    such that pp - 1 has sufficiently high 2-valuation.
    Used in sd_fft_trunc, sd_ifft_trunc, sd_fft_ctx_point_mul, etc.
*/
void sd_fft_ctx_init_prime(sd_fft_ctx_t Q, ulong pp)
{
    ulong N, i, k, l, init_depth, two_power_depth;
    double * t;
    double n, ninv;

    if (!fft_small_mulmod_satisfies_bounds(pp))
        flint_throw(FLINT_ERROR, "FFT prime %wu does not satisfy bounds for arithmetic", pp);

    Q->p = pp;
    Q->pinv = 1.0/Q->p;
    nmod_init(&Q->mod, pp);
    two_power_depth = n_trailing_zeros(pp - 1);
    if (two_power_depth == 0)
        flint_throw(FLINT_ERROR, "Input %wu is either 2 or not a prime", pp);
    Q->primitive_2power_root = sd_fft_ctx_primitive_2power_root(pp, two_power_depth, Q->mod);
    init_depth = n_min(two_power_depth, SD_FFT_CTX_W2TAB_INIT);

    n = Q->p;
    ninv = Q->pinv;

    /*
        fill wtab to a depth of init_depth:
        2^(init_depth-1) entries: 1, e(1/4), e(1/8), e(3/8), ...

        Q->w2tab[j] is itself a table of length 2^(j-1) containing 2^(j+1) st
        roots of unity. More documentation on the layout of w2tab can be found
        before the definition of SD_FFT_CTX_W2TAB_SIZE.

        All entries in w2tab are exactly-representable integers modulo pp, but
        they're stored as `double` to make use of the vectorized functions in
        machine_vectors.h.
    */
    N = n_pow2(init_depth - 1);
    t = (double*) flint_aligned_alloc(4096, n_round_up(N*sizeof(double), 4096));

    Q->w2tab[0] = t;
    t[0] = 1;
    for (k = 1, l = 1; k < init_depth; k++, l *= 2)
    {
        ulong ww = sd_fft_ctx_w2tab_root(Q, two_power_depth, k);
        double w = vec1d_set_d(vec1d_reduce_0n_to_pmhn(ww, n));
        double* curr = t + l;
        Q->w2tab[k] = curr;
        i = 0; do {
            vec1d x = vec1d_load(t + i);
            x = vec1d_mulmod(x, w, n, ninv);
            x = vec1d_reduce_pm1n_to_pmhn(x, n);
            vec1d_store(curr + i, x);
        } while (i += 1, i < l);
    }

#if FLINT_USES_PTHREAD
    atomic_init(&Q->w2tab_depth, (unsigned int)k);
#else
    Q->w2tab_depth = (unsigned int)k;
#endif

    /* the rest of the tables are uninitialized */
    for ( ; k < SD_FFT_CTX_W2TAB_SIZE; k++)
        Q->w2tab[k] = NULL;

#if FLINT_USES_PTHREAD
    pthread_mutex_init(&Q->mutex, NULL);
#endif

#if FLINT_WANT_ASSERT
    for (k = 1; k < init_depth; k++)
    {
        ulong ww = sd_fft_ctx_w2tab_root(Q, two_power_depth, k);
        for (i = 0; i < n_pow2(k-1); i++)
        {
            ulong www = nmod_pow_ui(ww, n_revbin(i+n_pow2(k-1), k), Q->mod);
            FLINT_ASSERT(Q->w2tab[k][i] == vec1d_reduce_0n_to_pmhn(www, n));
        }
    }
#endif
}

void sd_fft_ctx_fit_depth_with_lock(sd_fft_ctx_t Q, ulong depth)
{
    ulong two_power_depth = n_trailing_zeros(Q->mod.n - 1);

    if (depth > two_power_depth)
        flint_throw(FLINT_ERROR, "FFT prime %wu does not support depth %wu", Q->mod.n, depth);

#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&Q->mutex);
#endif

#if FLINT_USES_PTHREAD
    ulong k = (ulong)atomic_load(&Q->w2tab_depth);
#else
    ulong k = (ulong)Q->w2tab_depth;
#endif

    while (k < depth)
    {
        ulong i, j, l, off;
        ulong ww = sd_fft_ctx_w2tab_root(Q, two_power_depth, k);
        vec8d w    = vec8d_set_d(vec1d_reduce_0n_to_pmhn(ww, Q->p));
        vec8d n    = vec8d_set_d(Q->p);
        vec8d ninv = vec8d_set_d(Q->pinv);
        ulong N = n_pow2(k - 1);
        double* curr = (double*) flint_aligned_alloc(4096, n_round_up(N*sizeof(double), 4096));
        double* t = Q->w2tab[0];
        Q->w2tab[k] = curr;

        /* The first few tables are stored consecutively, so vec16 is ok. */
        off = 0;
        l = n_pow2(SD_FFT_CTX_W2TAB_INIT - 1);
        for (j = SD_FFT_CTX_W2TAB_INIT - 1; j < k; j++)
        {
            i = 0; do {
                vec8d x0 = vec8d_load_aligned(t + i + 0);
                vec8d x1 = vec8d_load_aligned(t + i + 8);
                x0 = vec8d_mulmod(x0, w, n, ninv);
                x1 = vec8d_mulmod(x1, w, n, ninv);
                x0 = vec8d_reduce_pm1n_to_pmhn(x0, n);
                x1 = vec8d_reduce_pm1n_to_pmhn(x1, n);
                vec8d_store_aligned(curr + off + i + 0, x0);
                vec8d_store_aligned(curr + off + i + 8, x1);
            } while (i += 16, i < l);
            FLINT_ASSERT(i == l);
            t = Q->w2tab[j + 1];
            l += off;
            off = l;
        }

#if FLINT_WANT_ASSERT
        {
            ulong ww = sd_fft_ctx_w2tab_root(Q, two_power_depth, k);
            for (i = 0; i < n_pow2(k-1); i++)
            {
                ulong www = nmod_pow_ui(ww, n_revbin(i+n_pow2(k-1), k), Q->mod);
                FLINT_ASSERT(Q->w2tab[k][i] == vec1d_reduce_0n_to_pmhn(www, Q->p));
            }
        }
#endif

        k++;
#if FLINT_USES_PTHREAD
        atomic_store(&Q->w2tab_depth, (unsigned int)k);
#else
        Q->w2tab_depth = (unsigned int)k;
#endif
    }


#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&Q->mutex);
#endif
}
