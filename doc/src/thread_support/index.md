---
title: thread_support.h
header: thread_support.h
description: Thread limits, parallel loops, parallel binary splitting and lock-free shared caches.
extra-headers: flint thread_pool
category: Foundations
order: 4
---

FLINT runs on one thread unless told otherwise. `flint_set_num_threads(n)`
allows FLINT functions to use up to $n$ threads, and many functions then
parallelise themselves: large integer and polynomial multiplication,
integer factoring, multivariate polynomial arithmetic and others. This
module provides the helpers that such functions are built from: thread
limits, the requesting of worker threads from the global pool, parallel
loops, parallel binary splitting, and shared caches for data computed
once and read by all threads. The
parallel helpers are experimental and their interfaces may change. The
pool itself is described in [](thread_pool); the guide
[](thread_pool/manual) shows how to use requested threads directly. The
user guide page [](guide/conventions/threading) explains which
operations are parallelised and the rules for calling FLINT from several
threads.

::: note Threads and workers
A function running on the *master* thread that may use $n$ threads in
total may start $n - 1$ *workers*. `flint_get_num_threads` returns the
number of threads the current thread may use, including itself, so its
value differs between threads. Library code should request at most
`flint_get_num_threads()` threads.
:::

::: summary algorithms/parallel :::

## Thread limits {#limits}

These functions set and query how many threads FLINT may use. All of them
except `flint_get_num_available_threads` and `FLINT_DEFAULT_THREAD_LIMIT`
are declared in `flint.h`.

::: function
void flint_set_num_threads(int num_threads)

Allows FLINT to use up to `num_threads` threads, resizing the global thread
pool to `num_threads - 1` workers. Must be called from the master thread,
outside any parallel region. Throws an exception if threads of the pool are
in use.

::: details
The first call creates `global_thread_pool`. Later calls resize it, which
joins all existing worker threads and creates new ones, even if the size is
unchanged. Without pthreads the value is ignored and FLINT uses one thread.
:::
:::

::: function
int flint_get_num_threads(void)

Returns the number of threads the calling thread may use, counting itself.
On the master thread this is the value set by `flint_set_num_threads`
(default 1); on a worker it is one more than the `max_workers` value it was
woken with.
:::

::: function
int flint_set_num_workers(int num_workers)
void flint_reset_num_workers(int num_workers)

Temporarily lowers the number of workers the calling thread may start to
`num_workers` (a larger value has no effect), and restores it.
`flint_set_num_workers` returns the previous value, which must be passed to
`flint_reset_num_workers` before the thread finishes its work.
:::

::: function
void _flint_set_num_workers(int num_workers)

Sets the number of workers the calling thread may start, without the
check of `flint_set_num_workers`. Used by the pool when waking a thread.
:::

::: function
slong flint_get_num_available_threads(void)

Returns one more than the number of threads of the global pool that are
currently available, or `flint_get_num_threads()` if the pool has not been
created.

::: details
The count is a snapshot and ignores the worker limit of the calling thread:
it can exceed `flint_get_num_threads()` on a worker thread. The extra 1
counts the calling thread, so `flint_get_num_available_threads` returns 1
when every thread of the pool is in use.
:::
:::

::: function
int flint_set_thread_affinity(int * cpus, slong length)
int flint_restore_thread_affinity(void)

Pin the master thread to `cpus[0]` and the global pool's threads to
`cpus[1]`, ..., `cpus[length - 1]`, or restore the original affinity.
Return 0 on success and nonzero on failure, including when the global pool
does not exist yet.
:::

::: macro
FLINT_DEFAULT_THREAD_LIMIT

A thread limit larger than any useful number of threads (99999), for
passing to functions that take a limit when no limit is wanted.
:::

## Requesting threads {#request}

::: function
slong flint_request_threads(thread_pool_handle ** handles, slong thread_limit)
void flint_give_back_threads(thread_pool_handle * handles, slong num_handles)

`flint_request_threads` obtains workers from the global pool for a
computation that may use `thread_limit` threads in total, counting the
caller. It stores an array of handles in `*handles` and returns the
number $w$ of workers obtained: 0 if $\mathit{thread\_limit} \le 1$, and
otherwise at most
$\min(\mathit{thread\_limit}, \texttt{flint\_get\_num\_threads()}) - 1$.
`flint_give_back_threads` returns the workers and frees the array.

::: details
These are wrappers of `thread_pool_request` and `thread_pool_give_back`
for `global_thread_pool`. `flint_request_threads` may obtain fewer workers
than allowed, or none, when other threads hold part of the pool. It may
allocate the array even when it returns $w = 0$, so always call
`flint_give_back_threads(handles, w)`; it accepts a `NULL` array. The
guide [](thread_pool/manual) shows how to wake, wait for and give back
the workers.
:::
:::

## Parallel loops {#parallel-do}

::: type
do_func_t

The type `void (* do_func_t)(slong i, void * args)` of a loop body.
:::

::: function
void flint_parallel_do(do_func_t f, void * args, slong n, int thread_limit, int flags)

Calls `f(i, args)` for $0 \le i < n$, distributing the calls over at most
`thread_limit` threads including the caller (`flint_get_num_threads()` if
$\mathit{thread\_limit} \le 0$). The calls for different $i$ may run in
any order and concurrently. `flags` chooses the scheduling.
:::

::: macro
FLINT_PARALLEL_UNIFORM
FLINT_PARALLEL_STRIDED
FLINT_PARALLEL_DYNAMIC
FLINT_PARALLEL_VERBOSE

Flags for `flint_parallel_do`. `UNIFORM` gives each thread one contiguous
block, which suits calls of equal cost. `STRIDED` gives thread $t$ of $T$
the indices $i \equiv t \pmod T$, which balances costs that grow or shrink
with $i$. `DYNAMIC` is not implemented and currently behaves like
`UNIFORM`. `VERBOSE` prints the schedule.
:::

## Parallel binary splitting {#bsplit}

::: type
bsplit_basecase_func_t
bsplit_merge_func_t
bsplit_init_func_t
bsplit_clear_func_t

Callback types for `flint_parallel_binary_splitting`:
`void basecase(void * res, slong a, slong b, void * args)`,
`void merge(void * res, void * x, void * y, void * args)`,
`void init(void * x, void * args)` and `void clear(void * x, void * args)`.
:::

::: function
void flint_parallel_binary_splitting(void * res, bsplit_basecase_func_t basecase, bsplit_merge_func_t merge, size_t sizeof_res, bsplit_init_func_t init, bsplit_clear_func_t clear, void * args, slong a, slong b, slong basecase_cutoff, int thread_limit, int flags)

Sets `res` to $f(a) \circ f(a+1) \circ \cdots \circ f(b-1)$ for an
associative operation $\circ$, by recursive halving of $[a, b)$ with the
two halves computed in parallel when workers are available. The intended
number of threads running base cases at the same time, the calling
thread included, is at most `thread_limit` (`flint_get_num_threads()` if
$\mathit{thread\_limit} \le 0$).
:::

::: macro
FLINT_PARALLEL_BSPLIT_LEFT_INPLACE

Flag for `flint_parallel_binary_splitting`: the merge function supports
`res` aliasing its left operand.
:::

## Shared caches {#shared-cache}

A shared cache holds one process-wide, read-only *version* of some lazily
computed data, such as a table of constants. Any thread reads the current
version without locking. A thread that needs more data builds a larger
version under the cache's lock and publishes it. Superseded versions are
kept until `flint_cleanup_master`, since other threads may still be
reading them. The `mp_real` module keeps its reduction tables in shared
caches.

::: summary algorithms/parallel/shared-cache :::

::: type
flint_shared_cache_struct

A cache: its lock, the current version, the list of retired versions, the
function freeing a version, and a link in the list of caches freed by
`flint_cleanup_master`. Caches are static objects initialised with
`FLINT_SHARED_CACHE_INIT`.
:::

::: macro
FLINT_SHARED_CACHE_INIT(free_version)

Static initialiser of a cache whose versions are freed by
`void free_version(void * version)`, for example
`static flint_shared_cache_struct cache = FLINT_SHARED_CACHE_INIT(my_free);`.
:::

::: type
flint_shared_cache_free_t
flint_shared_cache_sufficient_t
flint_shared_cache_build_t

Callback types: `void free(void * version)`;
`int sufficient(const void * version, void * args)`, whether a version
satisfies the request `args`; and
`void * build(const void * version, void * args)`, which returns a new
version satisfying the request, given the current version (or `NULL`),
which it must not modify.
:::

::: function
void * flint_shared_cache_ensure(flint_shared_cache_struct * cache, flint_shared_cache_sufficient_t sufficient, flint_shared_cache_build_t build, void * args)

Returns a version satisfying the request `args`: the current version if it
is sufficient, and otherwise one built and published under the lock after
checking the current version again. The returned version stays valid until
`flint_cleanup_master` or `flint_shared_cache_clear`.
:::

::: function
void * flint_shared_cache_get(flint_shared_cache_struct * cache)

Returns the current version, or `NULL` if none has been published. Costs
one atomic load with acquire semantics.
:::

::: function
void flint_shared_cache_lock(flint_shared_cache_struct * cache)
void flint_shared_cache_unlock(flint_shared_cache_struct * cache)
void flint_shared_cache_publish(flint_shared_cache_struct * cache, void * version)

The primitives that `flint_shared_cache_ensure` is built from. Under the
lock, `flint_shared_cache_publish` makes `version` current, retires the
previous version and registers the cache for `flint_cleanup_master`.
:::

::: function
void flint_shared_cache_clear(flint_shared_cache_struct * cache)

Frees all versions of the cache, leaving it empty. No other thread may be
using the cache. Intended for tests and tuning programs that measure cold
builds.
:::

::: function
void _flint_shared_cache_cleanup_all(void)
slong _flint_shared_cache_count(void)

Free every registered cache (called by `flint_cleanup_master`), and return
the number of registered caches (for tests).
:::
