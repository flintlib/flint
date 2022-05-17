/* 
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <string.h> /* for memcpy */
#undef ulong

#include "fft_small.h"
#include "machine_vectors.h"

void* flint_aligned_alloc(ulong alignment, ulong size)
{
    void* p = NULL;
    if (posix_memalign(&p, alignment, size))
    {
        flint_printf("Exception (FLINT memory_manager). Unable to allocate "
                     "%wu bytes with alignment %wu.\n", size, alignment);
        fflush(stdout);
        flint_abort();
    }
    return p;
}

void flint_aligned_free(void* p)
{
    free(p);
}

void sd_fft_ctx_init_prime(sd_fft_ctx_t Q, ulong pp)
{
    Q->blk_sz = BLK_SZ;
    Q->data = NULL;
    Q->p = pp;
    Q->pinv = 1.0/Q->p;
    nmod_init(&Q->mod, pp);
    Q->primitive_root = n_primitive_root_prime(pp);

    /* fill wtab to a depth of 10 (512 entries: 1, e(1/4), e(1/8), e(3/8), ...) */
    Q->wtab_depth = 10;
    ulong N = n_pow2(Q->wtab_depth-1);

    Q->w2s = (double*)flint_aligned_alloc(4096, n_round_up(N*sizeof(double), 4096));

    ulong w = nmod_pow_ui(Q->primitive_root, (Q->mod.n - 1)>>Q->wtab_depth, Q->mod);
    ulong wi = 1;
    for (ulong j = 0; j < N; j++)
    {
        ulong jr = n_revbin(j, Q->wtab_depth)/2;
        Q->w2s[jr] = vec1d_reduce_0n_to_pmhn(wi, Q->p);
        wi = nmod_mul(wi, w, Q->mod);
    }
}

void sd_fft_ctx_fit_wtab(sd_fft_ctx_t Q, ulong k)
{
    ulong N = n_pow2(Q->wtab_depth);
    double* oldw2s = Q->w2s;

    FLINT_ASSERT(oldw2s != NULL);

    if (Q->wtab_depth >= k)
        return;

    Q->w2s = (double*) flint_aligned_alloc(4096, n_pow2(k)*sizeof(double));
    memcpy(Q->w2s, oldw2s, N/2*sizeof(double));
    flint_aligned_free(oldw2s);

    while (Q->wtab_depth < k)
    {
        slong ww = nmod_pow_ui(Q->primitive_root, (Q->mod.n - 1)>>(Q->wtab_depth+1), Q->mod);
        vec8d w = vec8d_set_d(vec1d_reduce_0n_to_pmhn(ww, Q->p));
        vec8d n    = vec8d_set_d(Q->p);
        vec8d ninv = vec8d_set_d(Q->pinv);
        double* wptr = Q->w2s;
        ulong i = 0; do {
            vec8d x = vec8d_load_aligned(wptr + i);
            x = vec8d_mulmod2(x, w, n, ninv);
            x = vec8d_reduce_pm1n_to_pmhn(x, n);
            vec8d_store_aligned(wptr + N/2 + i, x);
        } while (i += 8, i < N/2);
        FLINT_ASSERT(i == N/2);
        Q->wtab_depth++;
        N *= 2;
    }
}

void sd_fft_ctx_set_depth(sd_fft_ctx_t Q, ulong l)
{
    Q->depth = l;
    if (l < LG_BLK_SZ)
    {
        flint_printf("depth %wu should be at least %d\n", l, LG_BLK_SZ);
        flint_abort();
    }
    sd_fft_ctx_fit_wtab(Q, l);
}
