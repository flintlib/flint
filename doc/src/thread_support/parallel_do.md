---
title: Running a loop on several threads
functions: flint_parallel_do do_func_t FLINT_PARALLEL_UNIFORM
---

`flint_parallel_do` is the simplest way to parallelise a computation that
splits into $n$ independent pieces of work. The loop body receives the
index $i$ and a pointer to arguments shared by all calls, and typically
writes its result to slot $i$ of an output array, so no locking is needed.
The work is divided among the threads before any of them starts, according
to the scheduling flag.

## Thread limit {#limit}

The loop uses at most `thread_limit` threads, counting the caller, and
never more than $n$. With $\mathit{thread\_limit} \le 0$ the limit is
`flint_get_num_threads()`, which is what library code should use. The
threads actually obtained may be fewer, when other computations hold some
of the pool; the result is the same, only slower. If no worker is
available, the loop runs serially on the calling thread.

## Scheduling {#scheduling}

Let $T$ be the number of threads obtained, including the calling thread,
which takes the last share ($t = T - 1$).

- With `FLINT_PARALLEL_UNIFORM` (and currently also with
  `FLINT_PARALLEL_DYNAMIC` or no flag), thread $t$ gets the contiguous
  block of indices $t s \le i < \min((t + 1) s, n)$, where
  $s = \lceil n/T \rceil$. The last blocks may be shorter or empty.
- With `FLINT_PARALLEL_STRIDED`, thread $t$ gets the indices
  $t, t + T, t + 2T, \ldots$ below $n$.
- `FLINT_PARALLEL_VERBOSE` may be combined with either flag. It prints
  $T$ and the share of each thread, unless the limit or $n$ is at most 1,
  in which case the loop runs serially and nothing is printed.

Strided scheduling is the right choice when the cost of `f(i, ...)` grows
or shrinks with $i$, for example when $i$ indexes coefficients of
increasing size: the work of any two threads then differs by at most the
cost of the most expensive call. Blocks are better when costs are equal, since each
thread touches a contiguous part of memory.

::: summary algorithms/parallel/parallel-loops :::

## Rules for the loop body {#rules}

- Calls for different indices may run concurrently and in any order.
  Writing to distinct memory is safe; anything shared must be read-only or
  protected by the caller.
- The body runs with a worker budget of zero on worker threads, so FLINT
  functions it calls run serially there. On the calling thread the budget
  is not lowered, so FLINT functions called there may take threads of the
  pool that are still free. Nested parallelism needs
  `flint_parallel_binary_splitting` or explicit use of the thread pool
  (see [](thread_pool/manual)).

## Example {#example}

This counts the primes below $10^6$ in 100 blocks of $10^4$ numbers.

```c example name=parallel_do_primes
#include <flint/flint.h>
#include <flint/thread_support.h>
#include <flint/ulong_extras.h>

#define BLOCKS 100
#define BLOCK_SIZE 10000

/* Count the primes in [i * BLOCK_SIZE, (i + 1) * BLOCK_SIZE). Each call
   writes only its own entry of the array, so no locking is needed. */
static void count_block(slong i, void * args)
{
    ulong * counts = args;
    ulong n, c = 0;

    for (n = i * BLOCK_SIZE; n < (i + 1) * BLOCK_SIZE; n++)
        c += n_is_prime(n);

    counts[i] = c;
}

int main(void)
{
    ulong counts[BLOCKS], total = 0;
    slong i;

    flint_set_num_threads(4);

    flint_parallel_do(count_block, counts, BLOCKS, 0, FLINT_PARALLEL_UNIFORM);

    for (i = 0; i < BLOCKS; i++)
        total += counts[i];

    flint_printf("primes below 10^6: %wu\n", total);
    flint_printf("in the last block: %wu\n", counts[BLOCKS - 1]);

    flint_cleanup_master();
    return 0;
}
```
```output
primes below 10^6: 78498
in the last block: 721
```

`flint_cleanup_master` at the end frees the thread pool and FLINT's
caches; it is optional but keeps memory checkers quiet.
