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

void sd_fft_ctx_init_prime(sd_fft_ctx_t Q, ulong pp)
{
    ulong N, i, k, l;
    double * t;
    double n, ninv;

    if (!fft_small_mulmod_satisfies_bounds(pp))
        flint_throw(FLINT_ERROR, "FFT prime %wu does not satisfy bounds for arithmetic", pp);

    Q->p = pp;
    Q->pinv = 1.0/Q->p;
    nmod_init(&Q->mod, pp);
    Q->primitive_root = n_primitive_root_prime(pp);

    n = Q->p;
    ninv = Q->pinv;

    /*
        fill wtab to a depth of SD_FFT_CTX_W2TAB_INIT:
        2^(SD_FFT_CTX_W2TAB_INIT-1) entries: 1, e(1/4), e(1/8), e(3/8), ...

        Q->w2tab[j] is itself a table of length 2^(j-1) containing 2^(j+1) st
        roots of unity.
    */
    N = n_pow2(SD_FFT_CTX_W2TAB_INIT - 1);
    t = (double*) flint_aligned_alloc(4096, n_round_up(N*sizeof(double), 4096));

    Q->w2tab[0] = t;
    t[0] = 1;
    for (k = 1, l = 1; k < SD_FFT_CTX_W2TAB_INIT; k++, l *= 2)
    {
        ulong ww = nmod_pow_ui(Q->primitive_root, (Q->mod.n - 1)>>(k + 1), Q->mod);
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
    for (k = 1; k < SD_FFT_CTX_W2TAB_INIT; k++)
    {
        ulong ww = nmod_pow_ui(Q->primitive_root, (Q->mod.n - 1)>>(k + 1), Q->mod);
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
        ulong ww = nmod_pow_ui(Q->primitive_root, (Q->mod.n - 1)>>(k + 1), Q->mod);
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
            ulong ww = nmod_pow_ui(Q->primitive_root, (Q->mod.n - 1)>>(k + 1), Q->mod);
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
