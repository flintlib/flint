/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "ulong_extras.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"

#if FLINT_USES_PTHREAD
#include <pthread.h>

/*
    Several threads share one lazy tower context (and a few elements of
    it): the operations of the context are serialized by its mutex, so
    the registry of towers is never mutated concurrently. Each thread
    checks radical identities of its own and of the shared elements.
*/

typedef struct
{
    gr_ctx_struct * K;
    gr_ptr * shared;
    slong num_shared;
    ulong seed;
    int failed;
}
thread_arg_struct;

static void *
_worker(void * argp)
{
    thread_arg_struct * arg = (thread_arg_struct *) argp;
    gr_ctx_struct * K = arg->K;
    flint_rand_t state;
    slong iter;

    flint_rand_init(state);
    flint_rand_set_seed(state, arg->seed, arg->seed * 7 + 1);

    for (iter = 0; iter < 30 && !arg->failed; iter++)
    {
        gr_ptr a, b, s, t, u;
        ulong p = 2 + n_randint(state, 30), q = 2 + n_randint(state, 30);
        fmpq_t e;

        a = gr_heap_init(K);
        b = gr_heap_init(K);
        s = gr_heap_init(K);
        t = gr_heap_init(K);
        u = gr_heap_init(K);
        fmpq_init(e);

        /* (sqrt p + sqrt q)^2 == p + q + 2 sqrt(p q) */
        GR_MUST_SUCCEED(gr_set_ui(a, p, K));
        GR_MUST_SUCCEED(gr_sqrt(a, a, K));
        GR_MUST_SUCCEED(gr_set_ui(b, q, K));
        GR_MUST_SUCCEED(gr_sqrt(b, b, K));
        GR_MUST_SUCCEED(gr_add(s, a, b, K));
        GR_MUST_SUCCEED(gr_sqr(s, s, K));
        GR_MUST_SUCCEED(gr_set_ui(t, p * q, K));
        GR_MUST_SUCCEED(gr_sqrt(t, t, K));
        GR_MUST_SUCCEED(gr_mul_ui(t, t, 2, K));
        GR_MUST_SUCCEED(gr_add_ui(t, t, p + q, K));
        if (gr_equal(s, t, K) != T_TRUE)
            arg->failed = 1;

        /* cbrt(x)^3 == x for a shared element x plus a random integer */
        {
            slong i = n_randint(state, arg->num_shared);
            GR_MUST_SUCCEED(gr_add_ui(u, arg->shared[i], n_randint(state, 5), K));
            fmpq_set_si(e, 1, 3);
            GR_MUST_SUCCEED(gr_pow_fmpq(t, u, e, K));
            GR_MUST_SUCCEED(gr_pow_ui(t, t, 3, K));
            if (gr_equal(t, u, K) != T_TRUE)
                arg->failed = 1;
        }

        /* the shared elements in a product with the thread's radicals */
        {
            slong i = n_randint(state, arg->num_shared);
            GR_MUST_SUCCEED(gr_mul(t, arg->shared[i], a, K));
            GR_MUST_SUCCEED(gr_div(t, t, a, K));
            if (gr_equal(t, arg->shared[i], K) != T_TRUE)
                arg->failed = 1;
        }

        /* repeated products with a shared element (lock-free with the
           dense forms, while other threads use it) */
        {
            slong i = n_randint(state, arg->num_shared), k, n = 2 + n_randint(state, 8);
            GR_MUST_SUCCEED(gr_set(u, arg->shared[i], K));
            for (k = 1; k < n; k++)
            {
                GR_MUST_SUCCEED(gr_mul(u, u, arg->shared[i], K));
                GR_MUST_SUCCEED(gr_add_si(u, u, 1, K));
                GR_MUST_SUCCEED(gr_sub_si(u, u, 1, K));
            }
            GR_MUST_SUCCEED(gr_pow_ui(t, arg->shared[i], n, K));
            if (gr_equal(t, u, K) != T_TRUE)
                arg->failed = 1;
            GR_MUST_SUCCEED(gr_div(u, u, arg->shared[i], K));
            GR_MUST_SUCCEED(gr_mul(u, u, arg->shared[i], K));
            if (gr_equal(t, u, K) != T_TRUE)
                arg->failed = 1;
        }

        fmpq_clear(e);
        gr_heap_clear(a, K);
        gr_heap_clear(b, K);
        gr_heap_clear(s, K);
        gr_heap_clear(t, K);
        gr_heap_clear(u, K);
    }

    flint_rand_clear(state);
    /* (the thread's caches) */
    flint_cleanup();
    return NULL;
}

TEST_FUNCTION_START(gr_tower_threads, state)
{
    gr_ctx_t QQ, K;
    slong iter;

    gr_ctx_init_fmpq(QQ);

    for (iter = 0; iter < 2 * flint_test_multiplier(); iter++)
    {
        pthread_t threads[4];
        thread_arg_struct args[4];
        gr_ptr shared[3];
        slong i, num_threads = 2 + n_randint(state, 3);

        gr_ctx_init_tower_lazy(K, QQ, n_randint(state, 2) ? GR_TOWER_MERGE_EXPRESS : 0);

        if (gr_ctx_is_threadsafe(K) != T_TRUE)
        {
            flint_printf("FAIL: context not threadsafe\n");
            flint_abort();
        }

        /* shared: sqrt 2, 1 + cbrt 5, pi + sqrt 3 */
        for (i = 0; i < 3; i++)
            shared[i] = gr_heap_init(K);
        GR_MUST_SUCCEED(gr_set_ui(shared[0], 2, K));
        GR_MUST_SUCCEED(gr_sqrt(shared[0], shared[0], K));
        {
            fmpq_t e;
            fmpq_init(e);
            fmpq_set_si(e, 1, 3);
            GR_MUST_SUCCEED(gr_set_ui(shared[1], 5, K));
            GR_MUST_SUCCEED(gr_pow_fmpq(shared[1], shared[1], e, K));
            GR_MUST_SUCCEED(gr_add_ui(shared[1], shared[1], 1, K));
            fmpq_clear(e);
        }
        GR_MUST_SUCCEED(gr_set_ui(shared[2], 3, K));
        GR_MUST_SUCCEED(gr_sqrt(shared[2], shared[2], K));
        GR_MUST_SUCCEED(gr_pi(shared[1], K));   /* (reused below) */
        GR_MUST_SUCCEED(gr_add(shared[2], shared[2], shared[1], K));
        GR_MUST_SUCCEED(gr_set_ui(shared[1], 5, K));
        {
            fmpq_t e;
            fmpq_init(e);
            fmpq_set_si(e, 1, 3);
            GR_MUST_SUCCEED(gr_pow_fmpq(shared[1], shared[1], e, K));
            GR_MUST_SUCCEED(gr_add_ui(shared[1], shared[1], 1, K));
            fmpq_clear(e);
        }

        {
            /* (an explicit stack size: some C libraries, such as musl,
               give threads small default stacks, and zero tests recurse) */
            pthread_attr_t attr;
            pthread_attr_init(&attr);
            pthread_attr_setstacksize(&attr, (size_t) 8 << 20);
            for (i = 0; i < num_threads; i++)
            {
                args[i].K = K;
                args[i].shared = shared;
                args[i].num_shared = 3;
                args[i].seed = n_randtest(state);
                args[i].failed = 0;
                if (pthread_create(threads + i, &attr, _worker, args + i) != 0)
                {
                    flint_printf("FAIL: pthread_create\n");
                    flint_abort();
                }
            }
            pthread_attr_destroy(&attr);
        }

        for (i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], NULL);
            if (args[i].failed)
            {
                flint_printf("FAIL: thread %wd\n", i);
                flint_abort();
            }
        }

        for (i = 0; i < 3; i++)
            gr_heap_clear(shared[i], K);
        gr_ctx_clear(K);
    }

    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}

#else

TEST_FUNCTION_START(gr_tower_threads, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}

#endif
