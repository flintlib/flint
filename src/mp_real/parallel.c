/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_pool.h"
#include "thread_support.h"
#include "mp_real.h"
#include "impl.h"

/* Thread helpers for the mp_real module's series.

   _mp_real_parallel_pair runs f1(a1) and f2(a2), the second on a pool
   thread when one is free, splitting the calling thread's budget (the
   flint_get_num_threads of the caller) between the two halves as
   flint_parallel_binary_splitting does, so that nested forks divide it
   further.  Returns 1 if a thread was used.  The pool threads are given
   back before returning, so that a merge following the call gets them
   for its multiplications.

   _mp_real_parallel_tasks runs f(i, args) for 0 <= i < n on
   w = min(n, threads) threads that take the tasks in order of i as they
   become free (tasks listed by decreasing cost balance best), each
   thread with a budget of threads / w for nested parallelism. */

int
_mp_real_parallel_pair(void (* f1)(void *), void * a1,
    void (* f2)(void *), void * a2)
{
    thread_pool_handle * threads;
    slong nt = flint_get_num_threads(), nw = 0;

    if (nt >= 2)
        nw = flint_request_threads(&threads, 2);

    if (nw == 0)
    {
        if (nt >= 2)
            flint_give_back_threads(threads, nw);
        f1(a1);
        f2(a2);
        return 0;
    }
    else
    {
        int save = flint_set_num_workers(nt - nt / 2 - 1);
        thread_pool_wake(global_thread_pool, threads[0], nt / 2 - 1, f2, a2);
        f1(a1);
        flint_reset_num_workers(save);
        thread_pool_wait(global_thread_pool, threads[0]);
        flint_give_back_threads(threads, nw);
        return 1;
    }
}

typedef struct
{
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
    slong next, n;
    void (* f)(slong, void *);
    void * args;
}
_tasks_struct;

static void
_tasks_worker(void * arg)
{
    _tasks_struct * S = (_tasks_struct *) arg;
    slong i;

    for (;;)
    {
#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&S->mutex);
#endif
        i = S->next++;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&S->mutex);
#endif
        if (i >= S->n)
            return;
        S->f(i, S->args);
    }
}

void
_mp_real_parallel_tasks(void (* f)(slong, void *), void * args, slong n)
{
    thread_pool_handle * handles;
    slong nt = flint_get_num_threads(), nw, i, budget;
    _tasks_struct S;

    if (n <= 0)
        return;

    nw = (nt >= 2 && n >= 2) ? flint_request_threads(&handles,
        FLINT_MIN(n, nt)) : 0;

    if (nw == 0)
    {
        if (nt >= 2 && n >= 2)
            flint_give_back_threads(handles, 0);
        for (i = 0; i < n; i++)
            f(i, args);
        return;
    }

    S.next = 0;
    S.n = n;
    S.f = f;
    S.args = args;
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&S.mutex, NULL);
#endif

    budget = nt / (nw + 1);
    {
        int save = flint_set_num_workers(budget - 1);
        for (i = 0; i < nw; i++)
            thread_pool_wake(global_thread_pool, handles[i], budget - 1,
                _tasks_worker, &S);
        _tasks_worker(&S);
        flint_reset_num_workers(save);
    }
    for (i = 0; i < nw; i++)
        thread_pool_wait(global_thread_pool, handles[i]);
    flint_give_back_threads(handles, nw);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&S.mutex);
#endif
}

/* the P, Q, T merge T = T1 Q2 + P1 T2, Q = Q1 Q2, P = P1 P2 (need_p),
   in place in (P, Q, T), destroying T2; with par the two independent
   halves {T1 Q2, Q1 Q2} and {P1 T2, P1 P2} on two threads (for the
   merges above the truncation, whose subtrees run one after the other
   with the whole thread budget) */
typedef struct
{
    mp_real_struct * P, * Q, * T, * P2, * Q2, * T2;
    int need_p;
    slong n;
}
_pqt_struct;

static void
_pqt_left(void * arg)
{
    _pqt_struct * A = (_pqt_struct *) arg;
    mp_real_mul(A->T, A->T, A->Q2, A->n);
    mp_real_mul(A->Q, A->Q, A->Q2, A->n);
}

static void
_pqt_right(void * arg)
{
    _pqt_struct * A = (_pqt_struct *) arg;
    mp_real_mul(A->T2, A->T2, A->P, A->n);
    if (A->need_p)
        mp_real_mul(A->P, A->P, A->P2, A->n);
}

void
_mp_real_pqt_merge(mp_real_t P, mp_real_t Q, mp_real_t T, mp_real_t P2, mp_real_t Q2,
    mp_real_t T2, int need_p, slong n, int par)
{
    _pqt_struct A;

    A.P = P; A.Q = Q; A.T = T; A.P2 = P2; A.Q2 = Q2; A.T2 = T2;
    A.need_p = need_p;
    A.n = n;

    if (par)
        _mp_real_parallel_pair(_pqt_left, &A, _pqt_right, &A);
    else
    {
        _pqt_left(&A);
        _pqt_right(&A);
    }
    mp_real_add(T, T, T2, n);
}
