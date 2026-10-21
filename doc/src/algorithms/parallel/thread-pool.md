---
title: The thread pool
summary: A fixed set of sleeping threads, handed out on request and woken with a function to run, with per-thread worker budgets that keep nested parallel code within the limit set by the user.
complexity: Waking or waiting for a thread costs one mutex and one condition-variable operation; requesting threads scans the pool under a lock.
order: 1
---

Creating a thread is expensive compared with many of the subcomputations
that FLINT runs in parallel. FLINT therefore creates its threads once,
in `flint_set_num_threads`, and keeps them asleep in a pool between uses.
This is the classical thread-pool design for POSIX threads [@But1997].

## Thread states {#states}

Each thread of a pool has an entry with a mutex, two condition variables
(`sleep1`, on which the thread waits for work, and `sleep2`, on which the
owner waits for completion), the function to run with its argument, a
`working` flag, an `available` flag and an `exit` flag.

| State | Set by | Meaning |
|---|---|---|
| available | `thread_pool_init`, `thread_pool_give_back` | free to be handed out |
| unavailable, sleeping | `thread_pool_request`, end of a function | owned by a caller, idle |
| unavailable, working | `thread_pool_wake` | running a function |

Every thread runs `thread_pool_idle_loop`:

::: algorithm Idle loop of a pool thread {#alg-idle}
1. Lock the entry. Set `working` to 0.
2. If `exit` is set, unlock and terminate.
3. Signal `sleep2` (waking a caller waiting in `thread_pool_wait`) and wait
   on `sleep1`, releasing the lock while waiting.
4. If `working` is still 0, the wake-up was spurious: go to step 2.
5. Unlock. Set this thread's worker budget to the `max_workers` value it
   was woken with, run the function, and go to step 1.
:::

`thread_pool_wake` sets the function, its argument and the budget, sets
`working` to 1 and signals `sleep1`, all under the entry's lock.
`thread_pool_wait` waits on `sleep2` until `working` is 0. Testing the flag
after every wake-up, rather than trusting the signal, is what makes the
protocol immune to spurious wake-ups and to a signal arriving before the
other side waits [@But1997, Chapter 3].

## Worker budgets {#budgets}

The pool itself does not know how many threads a computation may use.
Instead every thread has a thread-local *worker budget*: the number of
further workers it may start. `flint_get_num_threads()` returns the budget
plus one.

- The master thread's budget is set by `flint_set_num_threads(n)` to
  $n - 1$.
- A thread woken with `thread_pool_wake(..., max_workers, ...)` gets the
  budget `max_workers` while it runs.
- `flint_set_num_workers` lowers the budget of the calling thread
  temporarily, and `flint_reset_num_workers` restores it.

A parallel function that may use $t$ threads requests up to $t - 1$
workers. If it obtains $w$ of them, it divides its remaining budget among
itself and the workers it wakes. When every function follows this rule,
the number of threads working on one top-level call never exceeds the
budget of the caller plus one. Independently of the budgets, the number of
running threads never exceeds the pool size plus one, since threads are
only ever borrowed from the pool.

## Requesting threads {#request}

`thread_pool_request` scans the entries in order under the pool's mutex
and takes the first available ones. A caller may receive fewer threads
than it asked for, or none, when other threads hold parts of the pool.
FLINT's parallel functions treat the number of threads as a hint: the
work is split according to the threads actually obtained, and the result
does not depend on it.

## Resizing and affinity {#resize}

`thread_pool_set_size` succeeds only when every thread is available. It
tells all threads to exit, joins them, and creates the new set, so its cost
is that of creating the threads anew. `thread_pool_set_affinity` pins the
calling thread and the pool's threads to given processors with
`pthread_setaffinity_np`, where available; the original affinity of the
creating thread is recorded at initialisation and can be restored.

::: functions thread_pool_request thread_pool_wake thread_pool_wait flint_set_num_threads flint_request_threads :::
