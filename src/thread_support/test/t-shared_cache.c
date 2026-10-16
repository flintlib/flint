/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "ulong_extras.h"

/* a version: the squares below n */
typedef struct
{
    slong n;
    ulong * sq;
}
sq_struct;

static void
sq_free(void * v)
{
    sq_struct * S = v;
    flint_free(S->sq);
    flint_free(S);
}

static int
sq_sufficient(const void * v, void * args)
{
    return ((const sq_struct *) v)->n >= *(slong *) args;
}

/* geometric growth, copying the old version's entries */
static void *
sq_build(const void * old, void * args)
{
    const sq_struct * O = old;
    sq_struct * S = flint_malloc(sizeof(sq_struct));
    slong i, n = *(slong *) args;

    if (O != NULL)
        n = FLINT_MAX(n, O->n + O->n / 2);
    S->n = n;
    S->sq = flint_malloc(n * sizeof(ulong));
    i = 0;
    if (O != NULL)
    {
        for (; i < O->n; i++)
            S->sq[i] = O->sq[i];
    }
    for (; i < n; i++)
        S->sq[i] = (ulong) i * (ulong) i;
    return S;
}

static flint_shared_cache_struct sq_cache = FLINT_SHARED_CACHE_INIT(sq_free);

typedef struct
{
    const sq_struct ** seen;    /* a version each worker saw early */
    int * fail;
    ulong seed;
}
work_struct;

static void
worker(slong i, void * args)
{
    work_struct * W = args;
    flint_rand_t state;
    slong iter, k, n;
    const sq_struct * S, * first = NULL;

    flint_rand_init(state);
    flint_rand_set_seed(state, W->seed + i, 0);

    for (iter = 0; iter < 200; iter++)
    {
        n = 1 + n_randint(state, 1000);
        S = flint_shared_cache_ensure(&sq_cache, sq_sufficient, sq_build, &n);
        if (S->n < n)
            W->fail[i] |= 1;
        for (k = n_randint(state, n); k < n; k += 1 + n_randint(state, 50))
            if (S->sq[k] != (ulong) k * (ulong) k)
                W->fail[i] |= 2;
        if (first == NULL)
            first = S;
        /* the version seen first stays valid: retired, not freed */
        if (first->sq[first->n - 1] != (ulong) (first->n - 1) * (ulong) (first->n - 1))
            W->fail[i] |= 4;
        /* get returns the latest, at least as large */
        S = flint_shared_cache_get(&sq_cache);
        if (S == NULL || S->n < n)
            W->fail[i] |= 8;
    }

    W->seen[i] = first;
    flint_rand_clear(state);
}

TEST_FUNCTION_START(thread_support_shared_cache, state)
{
    slong iter, i, nthreads = 8;
    int old_threads = flint_get_num_threads();

    for (iter = 0; iter < 3 * flint_test_multiplier(); iter++)
    {
        const sq_struct * seen[8];
        int fail[8] = { 0 };
        work_struct W;
        slong n;

        /* a fresh cache, re-registering itself */
        flint_shared_cache_clear(&sq_cache);
        if (flint_shared_cache_get(&sq_cache) != NULL)
            TEST_FUNCTION_FAIL("clear\n");

        W.seen = seen;
        W.fail = fail;
        W.seed = n_randtest(state);

        flint_set_num_threads(nthreads);
        flint_parallel_do(worker, &W, nthreads, nthreads, FLINT_PARALLEL_UNIFORM);
        flint_set_num_threads(old_threads);

        for (i = 0; i < nthreads; i++)
            if (fail[i])
                TEST_FUNCTION_FAIL("worker %wd: %d\n", i, fail[i]);

        /* the primitives */
        n = 5000;
        flint_shared_cache_lock(&sq_cache);
        if (!sq_sufficient(flint_shared_cache_get(&sq_cache), &n))
            flint_shared_cache_publish(&sq_cache,
                sq_build(flint_shared_cache_get(&sq_cache), &n));
        flint_shared_cache_unlock(&sq_cache);
        if (((const sq_struct *) flint_shared_cache_get(&sq_cache))->n < n)
            TEST_FUNCTION_FAIL("publish\n");

        /* the versions seen by the workers are still intact */
        for (i = 0; i < nthreads; i++)
        {
            const sq_struct * S = seen[i];
            if (S != NULL && S->sq[S->n - 1] != (ulong) (S->n - 1) * (ulong) (S->n - 1))
                TEST_FUNCTION_FAIL("retired version\n");
        }

        /* flint_cleanup_master frees everything */
        flint_cleanup_master();
        if (flint_shared_cache_get(&sq_cache) != NULL || sq_cache.retired != NULL
            || sq_cache.registered)
            TEST_FUNCTION_FAIL("cleanup\n");
    }

    TEST_FUNCTION_END(state);
}
