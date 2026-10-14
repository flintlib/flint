/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fft_small.h"

/*
    Operand storage pools (see fft_small.h). The block is either a
    reservation in the creating thread's mpn_ctx buffer -- warm memory,
    which matters at a few hundred kB per operand, where fresh pages cost
    about a tenth of the transform time -- or one aligned allocation.
    Only the bookkeeping (the free-slot stack and the staging flag) is
    shared state; it is guarded by a mutex, so elements of one pool may
    be created and released on several threads. The memory itself is
    owned by the pool until it is cleared.
*/

#if FLINT_USES_PTHREAD
# define POOL_LOCK(pool) pthread_mutex_lock(&(pool)->mutex)
# define POOL_UNLOCK(pool) pthread_mutex_unlock(&(pool)->mutex)
#else
# define POOL_LOCK(pool)
# define POOL_UNLOCK(pool)
#endif

void
fft_small_op_pool_init(fft_small_op_pool_t pool)
{
    pool->R = NULL;
    pool->base = NULL;
    pool->slab_size = 0;
    pool->count = 0;
    pool->free_slots = NULL;
    pool->navail = 0;
    pool->stage = NULL;
    pool->stage_size = 0;
    pool->stage_busy = 0;
    pool->allocated = 0;
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&pool->mutex, NULL);
#endif
}

int
fft_small_op_pool_fit(fft_small_op_pool_t pool, mpn_ctx_t R,
                      ulong slab_size, slong count, ulong stage_bytes,
                      int flags)
{
    ulong total;
    slong i;

    FLINT_ASSERT(pool->base == NULL);
    FLINT_ASSERT(slab_size % FLINT_FFT_SMALL_ALIGNMENT == 0);

    if (count < 0)
        count = 0;

    stage_bytes = n_round_up(stage_bytes, FLINT_FFT_SMALL_ALIGNMENT);
    total = (ulong) count * slab_size + stage_bytes;

    if (total == 0)
        return 0;

    if (flags & FFT_SMALL_OP_POOL_RESERVE)
    {
        pool->base = (char *) mpn_ctx_fit_buffer_reserve(R, 0, total);
        if (pool->base != NULL)
            pool->R = R;
    }

    if (pool->base == NULL && (flags & FFT_SMALL_OP_POOL_ALLOC))
    {
        pool->base = flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT, total);
        pool->allocated = 1;
    }

    if (pool->base == NULL)
        return 0;

    pool->slab_size = slab_size;
    pool->count = count;
    pool->free_slots = flint_malloc(FLINT_MAX(count, 1) * sizeof(slong));
    /* hand out the lowest slabs first */
    for (i = 0; i < count; i++)
        pool->free_slots[i] = count - 1 - i;
    pool->navail = count;
    pool->stage = (stage_bytes != 0) ? pool->base + (ulong) count * slab_size : NULL;
    pool->stage_size = stage_bytes;
    pool->stage_busy = 0;
    return 1;
}

void
fft_small_op_pool_clear(fft_small_op_pool_t pool)
{
    FLINT_ASSERT(pool->navail == pool->count);
    FLINT_ASSERT(!pool->stage_busy);

    if (pool->R != NULL)
        mpn_ctx_fit_buffer_release(pool->R);
    else if (pool->allocated)
        flint_aligned_free(pool->base);

    flint_free(pool->free_slots);
#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&pool->mutex);
#endif
}

int
fft_small_op_pool_owns(const fft_small_op_pool_t pool, const void * d)
{
    const char * c = (const char *) d;
    return pool->base != NULL && c >= pool->base
        && c < pool->base + (ulong) pool->count * pool->slab_size
               + pool->stage_size;
}

double *
fft_small_op_pool_take(fft_small_op_pool_t pool)
{
    double * d = NULL;

    if (pool->base == NULL)
        return NULL;

    POOL_LOCK(pool);
    if (pool->navail > 0)
    {
        slong slot = pool->free_slots[--pool->navail];
        d = (double *) (pool->base + (ulong) slot * pool->slab_size);
    }
    POOL_UNLOCK(pool);

    return d;
}

void
fft_small_op_pool_give(fft_small_op_pool_t pool, double * d)
{
    slong slot;

    FLINT_ASSERT(fft_small_op_pool_owns(pool, d));
    slot = (slong) ((ulong) ((char *) d - pool->base) / pool->slab_size);

    POOL_LOCK(pool);
    FLINT_ASSERT(pool->navail < pool->count);
    pool->free_slots[pool->navail++] = slot;
    POOL_UNLOCK(pool);
}

ulong *
fft_small_op_pool_scratch_take(fft_small_op_pool_t pool, ulong limbs,
                               int use_slabs)
{
    ulong bytes = limbs * sizeof(ulong);
    ulong * t = NULL;

    if (pool->base != NULL)
    {
        POOL_LOCK(pool);
        if (pool->stage != NULL && !pool->stage_busy && bytes <= pool->stage_size)
        {
            pool->stage_busy = 1;
            t = (ulong *) pool->stage;
        }
        else if (use_slabs && bytes <= pool->slab_size && pool->navail > 0)
        {
            slong slot = pool->free_slots[--pool->navail];
            t = (ulong *) (pool->base + (ulong) slot * pool->slab_size);
        }
        POOL_UNLOCK(pool);
    }

    if (t == NULL)
        t = flint_malloc(FLINT_MAX(bytes, sizeof(ulong)));

    return t;
}

void
fft_small_op_pool_scratch_give(fft_small_op_pool_t pool, ulong * t)
{
    if (pool->stage != NULL && (char *) t == pool->stage)
    {
        POOL_LOCK(pool);
        pool->stage_busy = 0;
        POOL_UNLOCK(pool);
    }
    else if (fft_small_op_pool_owns(pool, t))
        fft_small_op_pool_give(pool, (double *) t);
    else
        flint_free(t);
}
