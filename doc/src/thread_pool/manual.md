---
title: Using the thread pool directly
summary: The request, wake, wait and give-back pattern behind all parallel code in FLINT, how to divide a thread budget between nested workers, and a complete example.
order: 1
---

Parallel code in FLINT follows one pattern: obtain workers, wake each with
a function and its argument, do a share of the work on the calling thread,
wait for the workers, and give them back. `flint_parallel_do` and
`flint_parallel_binary_splitting` implement this pattern; use the pool
directly when the work does not fit either of them, for example when the
threads must do different things. This guide describes the pattern for
the global pool, which is the pool that all FLINT functions use.

## The pattern {#pattern}

1. Call `flint_request_threads(&handles, limit)` with the number of threads
   the computation may use, counting the calling thread. It returns the
   number $w$ of workers obtained, which may be smaller than requested and
   may be 0.
2. Split the work into $w + 1$ shares. Wake worker $i$ with
   `thread_pool_wake(global_thread_pool, handles[i], 0, f, &share[i])`.
   The call returns at once, and `f` runs on the worker.
3. Process the last share on the calling thread.
4. Call `thread_pool_wait(global_thread_pool, handles[i])` for each
   worker, then `flint_give_back_threads(handles, w)`. Call
   `flint_give_back_threads` even if $w = 0$: the array of handles may
   have been allocated anyway.

Never assume that $w$ equals the number requested: the pool is shared, and
FLINT may have been built without thread support. Split the work according
to the $w$ actually obtained, so that the result does not depend on it.

::: summary algorithms/parallel/thread-pool :::

## Nested parallelism {#nested}

The third argument of `thread_pool_wake` is the worker budget of the woken
thread: while the function runs, `flint_get_num_threads()` returns that
value plus one on the worker, so FLINT functions called there may start
that many further workers. Pass 0 unless the function is itself parallel.

The pool does not check budgets. To keep a computation within the limit
$n$ = `flint_get_num_threads()` of the calling thread, the $w$ workers,
their budgets and the budget the calling thread keeps for its own share
should add up to at most $n - 1$. Lower the calling thread's budget with
`flint_set_num_workers` before processing its share, and restore it with
`flint_reset_num_workers` afterwards. For example,
`flint_parallel_binary_splitting` with a limit of $n$ threads wakes one
worker with budget $\lfloor n/2 \rfloor - 1$ for the right half and keeps
$n - \lfloor n/2 \rfloor - 1$ for the left half.

## Example {#example}

This sums $0 + 1 + \cdots + 65535$ with up to four threads. Each thread
sums a contiguous range, and the calling thread takes the last range.

```c example name=thread_pool_sum
#include <flint/flint.h>
#include <flint/thread_pool.h>
#include <flint/thread_support.h>

typedef struct
{
    ulong start;
    ulong stop;
    ulong sum;
}
work_t;

static void worker(void * arg)
{
    work_t * w = arg;
    ulong i, s = 0;

    for (i = w->start; i < w->stop; i++)
        s += i;

    w->sum = s;
}

int main(void)
{
    thread_pool_handle * handles;
    work_t work[4];
    slong i, num_workers, num_threads;
    ulong n = 65536, total = 0;

    flint_set_num_threads(4);

    /* Ask for workers for a computation using at most 4 threads in total.
       We may get fewer than 3, for example if FLINT has no pthreads. */
    num_workers = flint_request_threads(&handles, 4);
    num_threads = num_workers + 1;

    /* Thread i sums [i n / num_threads, (i + 1) n / num_threads). */
    for (i = 0; i < num_threads; i++)
    {
        work[i].start = (i * n) / num_threads;
        work[i].stop = ((i + 1) * n) / num_threads;
    }

    for (i = 0; i < num_workers; i++)
        thread_pool_wake(global_thread_pool, handles[i], 0, worker, &work[i]);

    worker(&work[num_workers]);     /* the last share runs here */

    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);

    for (i = 0; i < num_threads; i++)
        total += work[i].sum;

    flint_printf("0 + 1 + ... + %wu = %wu\n", n - 1, total);
    flint_printf("expected:          %wu\n", (n / 2) * (n - 1));

    flint_cleanup_master();
    return 0;
}
```
```output
0 + 1 + ... + 65535 = 2147450880
expected:          2147450880
```
