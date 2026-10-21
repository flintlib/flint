---
title: thread_pool.h
header: thread_pool.h
description: The pool of worker threads that all parallel code in FLINT draws from.
extra-headers: flint
category: Foundations
order: 3
---

FLINT keeps one process-wide pool of sleeping worker threads,
`global_thread_pool`, created by `flint_set_num_threads`. A parallel
function borrows threads from the pool, gives each one a function to run,
waits for them, and returns them. Most code should not call this module
directly: `flint_parallel_do` and `flint_parallel_binary_splitting` in
[](thread_support) cover the common patterns, and `flint_request_threads`
handles the bookkeeping for the global pool. The guide
[](thread_pool/manual) shows how to run work on the pool's threads when
neither pattern fits. The rules for using threads with FLINT in general
are in [](guide/conventions/threading).

::: summary algorithms/parallel/thread-pool :::

## Conventions {#conventions}

A thread in the pool is *available* or *unavailable* (handed out by
`thread_pool_request` and not yet given back), and independently
*sleeping* or *working* (between `thread_pool_wake` and the end of its
function). A handle may only be woken while unavailable and sleeping, and
only given back while sleeping. None of this is checked except by
`FLINT_ASSERT` in builds configured with `--enable-assert`.

Without pthreads (FLINT configured with `--disable-pthread`), no thread
can run the work handed to a worker, so every pool has size 0:
`thread_pool_init` and `thread_pool_set_size` treat any size as 0, and
`thread_pool_request` returns no handles. In particular
`flint_request_threads` obtains no workers and parallel code runs
serially.

## Types {#types}

::: type
thread_pool_t

A thread pool: an array of thread entries and a mutex protecting their
states. The pool owns its threads; they are created by `thread_pool_init`
or `thread_pool_set_size` and joined by `thread_pool_clear`.
:::

::: type
thread_pool_handle

An index identifying a thread of a pool (an `int`, declared in `flint.h`).
Handles are obtained from `thread_pool_request`.
:::

::: variable
thread_pool_t global_thread_pool
int global_thread_pool_initialized

The pool used by all of FLINT, and whether it has been initialised. It is
created by the first call to `flint_set_num_threads` and resized by later
calls.
:::

## Creating and resizing {#lifecycle}

::: function
void thread_pool_init(thread_pool_t T, slong size)

Initialises `T` with `size` sleeping, available threads, or no threads if
$\mathit{size} \le 0$ or FLINT was built without pthreads.
:::

::: function
slong thread_pool_get_size(thread_pool_t T)

Returns the number of threads in `T`.
:::

::: function
int thread_pool_set_size(thread_pool_t T, slong new_size)

If every thread of `T` is available, replaces them by `new_size` new
threads (none if $\mathit{new\_size} \le 0$ or FLINT was built without
pthreads) and returns 1. Otherwise
changes nothing and returns 0.

::: details
All old threads are told to exit and joined, even if the size does not
change, and new threads are created. Resizing therefore costs a thread
creation per thread and should not be done in a loop.
:::
:::

::: function
void thread_pool_clear(thread_pool_t T)

Joins all threads of `T` and frees its memory. All threads must have been
given back.
:::

## Running work {#work}

A computation on the pool's threads calls these four functions in order:
it requests threads, wakes each of them with a function to run, waits for
each of them, and gives them back. See [](thread_pool/manual) for the
complete pattern and an example.

::: function
slong thread_pool_request(thread_pool_t T, thread_pool_handle * out, slong requested)

Marks up to `requested` available threads as unavailable, writes their
handles to `out` and returns how many were obtained, possibly 0. The
threads must later be returned with `thread_pool_give_back`.

::: details
`thread_pool_request` scans the pool from the first thread under the pool
mutex, so the call costs time linear in the pool size. Callers should not
assume that they obtain as many threads as requested, since other threads
may hold some of the pool. The request is not limited by the calling
thread's worker budget (`flint_get_num_threads()` minus 1); for the global
pool, `flint_request_threads` applies that limit.
:::
:::

::: function
void thread_pool_wake(thread_pool_t T, thread_pool_handle i, int max_workers, void (*f)(void*), void * a)

Makes the sleeping thread `i` run `f(a)`. While running, that thread may
itself start at most `max_workers` further worker threads (it sees
`flint_get_num_threads()` as $\mathit{max\_workers} + 1$). Pass 0 unless
`f` is itself parallel.

::: details
`thread_pool_wake` returns without waiting; `f(a)` runs concurrently with
the caller until `thread_pool_wait` is called for `i`. The thread must
have been obtained with `thread_pool_request` and must be sleeping.
`max_workers` is not checked against the caller's own budget; how to
divide a budget between the caller and its workers is explained in
[](thread_pool/manual).
:::
:::

::: function
void thread_pool_wait(thread_pool_t T, thread_pool_handle i)

Blocks until thread `i` has finished the function it was woken with and
is sleeping again.
:::

::: function
void thread_pool_give_back(thread_pool_t T, thread_pool_handle i)

Marks the sleeping thread `i` as available again, so that later calls
of `thread_pool_request` can hand it out. If the thread was woken, call
`thread_pool_wait` for it first.
:::

## Processor affinity {#affinity}

::: function
int thread_pool_set_affinity(thread_pool_t T, int * cpus, slong length)
int thread_pool_restore_affinity(thread_pool_t T)

Pins the calling thread to processor `cpus[0]` and the threads of `T` to
`cpus[1]`, ..., `cpus[length - 1]` (threads beyond `length - 1` are left
alone), or restores the affinity that the calling thread had when `T` was
created. Return 0 on success and an error code from
`pthread_setaffinity_np` otherwise. Without CPU set support both functions
do nothing and return 0.

::: details
The affinity restored by `thread_pool_restore_affinity` is that of the
thread that called `thread_pool_init`, recorded at that point. It is
applied both to the pool's threads and to the calling thread, so the
function should be called from the thread that created the pool (for the
global pool, the master thread, which creates it in
`flint_set_num_threads`); any other thread receives the creator's
affinity rather than its own original one.
:::
:::

## Internal helpers {#internal}

::: function
void * thread_pool_idle_loop(void * varg)

The function executed by every pool thread: sleep on the thread's
condition variable, run the assigned function when woken, signal
completion, repeat until told to exit. Not for direct use.
:::

::: function
void _thread_pool_distribute_work_2(slong start, slong stop, slong * Astart, slong * Astop, slong Alen, slong * Bstart, slong * Bstop, slong Blen)

Splits the index range $[\mathit{start}, \mathit{stop})$ of the
concatenation of two arrays of lengths `Alen` and `Blen` into the
corresponding ranges $[\mathit{Astart}, \mathit{Astop})$ of the first array
and $[\mathit{Bstart}, \mathit{Bstop})$ of the second; empty ranges are
returned as $[0, 0)$. Requires
$0 \le \mathit{start} \le \mathit{stop} \le \mathit{Alen} + \mathit{Blen}$.
:::

::: function
ulong _thread_pool_find_work_2(ulong a, ulong alpha, ulong b, ulong beta, ulong yn, ulong yd)

Load balancing for work with two cost rates: if the first $a$ items cost
$\alpha$ each and the next $b$ items cost $\beta$ each, returns the number
$x$ of items whose total cost is the fraction $y_n/y_d$ of the whole,
rounded down. Requires $y_n \le y_d$. The intermediate product
$y_n (a \alpha + b \beta)$ must not overflow a word.
:::
