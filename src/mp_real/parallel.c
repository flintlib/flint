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
   thread with a budget of threads / w for nested parallelism; the
   _max version caps w at max_workers (the remaining budget going to
   the tasks), for tasks whose memory scales with their concurrency,
   and the _cost version takes the tasks' estimated costs and sizes the
   budgets by them (see below).

   _mp_real_vec_prod multiplies the entries of a vector (destroyed)
   by a product tree whose halves fork while they hold at least two
   entries each (and n reaches MP_REAL_PAR_MIN_LIMBS);
   _mp_real_vec_prod_complex does the same for a vector of complex
   numbers given as separate real and imaginary parts. */

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

/* a worker: its assigned first task (or none), then the queue */
typedef struct
{
    _tasks_struct * S;
    slong first;
}
_tasks_worker_struct;

#define _TASK_CALL(S, i) ((S)->f(i, (S)->args))

static void
_tasks_worker(void * arg)
{
    _tasks_worker_struct * W = (_tasks_worker_struct *) arg;
    _tasks_struct * S = W->S;
    slong i;

    if (W->first >= 0)
        _TASK_CALL(S, W->first);

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
        _TASK_CALL(S, i);
    }
}

/* Runs the n tasks on w = min(n, max_workers, threads) threads, the
   calling one included.  With cost NULL, the threads take the tasks
   in order of index as they become free, each with a budget of
   threads / w for nested parallelism.  With cost (the tasks' estimated
   costs, listed in decreasing order), the first w tasks are assigned
   one per thread with budgets in proportion to their costs -- so that
   a task much costlier than the others, such as the series remainder
   of a cascade, gets most of the threads for the parallelism inside
   it rather than a single one while the others idle -- and the
   remaining tasks go to whichever thread frees up first. */
static void
_parallel_tasks_core(void (* f)(slong, void *), void * args,
    slong n, slong max_workers, const double * cost)
{
    thread_pool_handle * handles;
    slong nt = flint_get_num_threads(), nw, i, want, w;
    slong budget[FLINT_BITS * 4] = { 0 };   /* all read entries are set below */
    _tasks_struct S;
    _tasks_worker_struct W[FLINT_BITS * 4];

    if (n <= 0)
        return;

    want = FLINT_MIN(n, max_workers);
    want = FLINT_MIN(want, FLINT_BITS * 4);
    nw = (nt >= 2 && want >= 2) ? flint_request_threads(&handles,
        FLINT_MIN(want, nt)) : 0;

    if (nw == 0)
    {
        if (nt >= 2 && want >= 2)
            flint_give_back_threads(handles, 0);
        for (i = 0; i < n; i++)
            f(i, args);
        return;
    }

    w = nw + 1;
    S.n = n;
    S.f = f;
    S.args = args;
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&S.mutex, NULL);
#endif

    if (cost == NULL)
    {
        for (i = 0; i < w; i++)
        {
            budget[i] = nt / w;
            W[i].first = -1;
        }
        S.next = 0;
    }
    else
    {
        /* budgets proportional to the costs of the first w tasks, at
           least one each, summing to nt (the rounding leftovers to the
           costliest) */
        double total = 0.0;
        slong sum = 0;
        for (i = 0; i < w; i++)
            total += FLINT_MAX(cost[i], 0.0);
        for (i = 0; i < w; i++)
        {
            budget[i] = (total > 0.0)
                ? (slong) ((double) (nt - w) * FLINT_MAX(cost[i], 0.0) / total) + 1
                : nt / w;
            sum += budget[i];
            W[i].first = i;
        }
        budget[0] += nt - sum;
        S.next = w;
    }

    for (i = 0; i < w; i++)
        W[i].S = &S;

    {
        int save = flint_set_num_workers((int) budget[0] - 1);
        for (i = 0; i < nw; i++)
            thread_pool_wake(global_thread_pool, handles[i],
                (int) budget[i + 1] - 1, _tasks_worker, &W[i + 1]);
        _tasks_worker(&W[0]);
        flint_reset_num_workers(save);
    }
    for (i = 0; i < nw; i++)
        thread_pool_wait(global_thread_pool, handles[i]);
    flint_give_back_threads(handles, nw);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&S.mutex);
#endif
}

void
_mp_real_parallel_tasks_cost(void (* f)(slong, void *), void * args,
    slong n, slong max_workers, const double * cost)
{
    _parallel_tasks_core(f, args, n, max_workers, cost);
}

/* Lanes: the n tasks (costs in decreasing order) are assigned to
   nlanes lanes by the greedy longest-processing-time rule, each task to
   the lane with the least cost so far (ties to the lowest index), and
   the lanes run as parallel tasks, each calling f(i, lane, args) for
   its tasks in increasing order of i.  The assignment depends only on
   n, nlanes and the costs -- not on how many threads the pool grants or
   on timing -- so per-lane state such as a running product gives the
   same result on every run. */
typedef struct
{
    void (* f)(slong, slong, void *);
    void * args;
    const slong * lane_of;
    slong n;
}
_lanes_struct;

static void
_lane_worker(slong lane, void * arg)
{
    _lanes_struct * L = (_lanes_struct *) arg;
    slong i;
    for (i = 0; i < L->n; i++)
        if (L->lane_of[i] == lane)
            L->f(i, lane, L->args);
}

void
_mp_real_parallel_lanes(void (* f)(slong, slong, void *), void * args,
    slong n, slong nlanes, const double * cost)
{
    slong lane_of[FLINT_BITS * 4];
    double load[FLINT_BITS * 4], lcost[FLINT_BITS * 4];
    slong order[FLINT_BITS * 4], rank[FLINT_BITS * 4];
    slong i, j, k;
    _lanes_struct L;

    if (n <= 0)
        return;
    FLINT_ASSERT(n <= FLINT_BITS * 4);
    nlanes = FLINT_MAX(1, FLINT_MIN(nlanes, n));

    for (j = 0; j < nlanes; j++)
        load[j] = 0.0;
    for (i = 0; i < n; i++)
    {
        k = 0;
        for (j = 1; j < nlanes; j++)
            if (load[j] < load[k])
                k = j;
        lane_of[i] = k;
        load[k] += (cost != NULL) ? cost[i] : 1.0;
    }

    /* the lanes as tasks, costliest first (a stable sort by load),
       renumbered so that lane j is the j-th costliest */
    for (j = 0; j < nlanes; j++)
        order[j] = j;
    for (j = 1; j < nlanes; j++)
        for (k = j; k > 0 && load[order[k]] > load[order[k - 1]]; k--)
            FLINT_SWAP(slong, order[k], order[k - 1]);
    for (j = 0; j < nlanes; j++)
    {
        rank[order[j]] = j;
        lcost[j] = load[order[j]];
    }
    for (i = 0; i < n; i++)
        lane_of[i] = rank[lane_of[i]];

    L.f = f;
    L.args = args;
    L.lane_of = lane_of;
    L.n = n;
    _parallel_tasks_core(_lane_worker, &L, nlanes, nlanes, lcost);
}

/* the number of lanes to use for n tasks under max_workers: bounded by
   the thread setting, which (unlike the threads the pool grants at the
   moment) is fixed for a given configuration */
slong
_mp_real_parallel_lanes_count(slong n, slong max_workers)
{
    slong w = FLINT_MIN(n, max_workers);
    w = FLINT_MIN(w, FLINT_BITS * 4);
    w = FLINT_MIN(w, flint_get_num_threads());
    return FLINT_MAX(w, 1);
}

void
_mp_real_parallel_tasks_max(void (* f)(slong, void *), void * args, slong n,
    slong max_workers)
{
    _mp_real_parallel_tasks_cost(f, args, n, max_workers, NULL);
}

void
_mp_real_parallel_tasks(void (* f)(slong, void *), void * args, slong n)
{
    _mp_real_parallel_tasks_cost(f, args, n, n, NULL);
}

/* the product trees */
typedef struct
{
    mp_real_struct * vec;
    mp_real_struct * vec2;
    slong len, n;
}
_prod_struct;

static void _vec_prod(mp_real_struct * vec, slong len, slong n);

static void
_prod_job(void * arg)
{
    _prod_struct * A = (_prod_struct *) arg;
    _vec_prod(A->vec, A->len, A->n);
}

/* vec[0] = prod vec[0..len), len >= 1 */
static void
_vec_prod(mp_real_struct * vec, slong len, slong n)
{
    slong m;

    if (len == 1)
        return;
    if (len == 2)
    {
        mp_real_mul(vec, vec, vec + 1, n);
        return;
    }
    m = len / 2;
    if (m >= 2 && n >= MP_REAL_PAR_MIN_LIMBS && flint_get_num_threads() >= 2)
    {
        _prod_struct L, R;
        L.vec = vec; L.len = m; L.n = n;
        R.vec = vec + m; R.len = len - m; R.n = n;
        _mp_real_parallel_pair(_prod_job, &L, _prod_job, &R);
    }
    else
    {
        _vec_prod(vec, m, n);
        _vec_prod(vec + m, len - m, n);
    }
    mp_real_mul(vec, vec, vec + m, n);
    /* release the consumed entry at once (the factors of a bit-burst
       evaluation are each a full-precision number) */
    mp_real_clear(vec + m);
    mp_real_init(vec + m);
}

void
_mp_real_vec_prod(mp_real_t res, mp_real_struct * vec, slong len, slong n)
{
    if (len == 0)
    {
        mp_real_set_ui(res, 1);
        return;
    }
    _vec_prod(vec, len, n);
    mp_real_swap(res, vec);
}

static void _vec_prod_complex(mp_real_struct * re, mp_real_struct * im, slong len, slong n);

static void
_prod_complex_job(void * arg)
{
    _prod_struct * A = (_prod_struct *) arg;
    _vec_prod_complex(A->vec, A->vec2, A->len, A->n);
}

static void
_vec_prod_complex(mp_real_struct * re, mp_real_struct * im, slong len, slong n)
{
    slong m;

    if (len == 1)
        return;
    if (len == 2)
    {
        mp_real_mul_complex(re, im, re, im, re + 1, im + 1, n);
        return;
    }
    m = len / 2;
    if (m >= 2 && n >= MP_REAL_PAR_MIN_LIMBS && flint_get_num_threads() >= 2)
    {
        _prod_struct L, R;
        L.vec = re; L.vec2 = im; L.len = m; L.n = n;
        R.vec = re + m; R.vec2 = im + m; R.len = len - m; R.n = n;
        _mp_real_parallel_pair(_prod_complex_job, &L, _prod_complex_job, &R);
    }
    else
    {
        _vec_prod_complex(re, im, m, n);
        _vec_prod_complex(re + m, im + m, len - m, n);
    }
    mp_real_mul_complex(re, im, re, im, re + m, im + m, n);
    mp_real_clear(re + m);
    mp_real_init(re + m);
    mp_real_clear(im + m);
    mp_real_init(im + m);
}

void
_mp_real_vec_prod_complex(mp_real_t rre, mp_real_t rim, mp_real_struct * re,
    mp_real_struct * im, slong len, slong n)
{
    if (len == 0)
    {
        mp_real_set_ui(rre, 1);
        mp_real_zero(rim);
        return;
    }
    _vec_prod_complex(re, im, len, n);
    mp_real_swap(rre, re);
    mp_real_swap(rim, im);
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
