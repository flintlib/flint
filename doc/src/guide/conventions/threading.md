---
title: Threads
summary: How to let FLINT use several threads, which operations are parallelised, the rules for calling FLINT from several threads, and the conventions for writing threaded functions.
order: 3
---

FLINT can use several threads for large computations, and it can be
called from several threads of a multithreaded program. The functions
that control threads are documented
in [](thread_support), the thread pool behind them in [](thread_pool),
and their design in [](algorithms/parallel).

## Letting FLINT use threads {#num-threads}

By default FLINT uses only the thread that calls it. After
`flint_set_num_threads(n)`, functions that are parallelised may use up
to $n$ threads in total, the calling thread included; the call creates
a global pool of $n - 1$ worker threads, which sleep while they are not
needed. `flint_get_num_threads()` returns the current limit. The limit
only applies to threads started by FLINT: an external library called by
FLINT, such as a multithreaded BLAS, may start threads of its own.

No function needs a special name or argument to run in parallel: when
the limit allows it and the input is large enough, the ordinary
functions use the thread pool by themselves. Among the operations that
are parallelised are:

- multiplication of large integers and polynomials with the FFTs of
  [](fft_small) and [](fft), and therefore everything built on them;
- matrix multiplication over $\Z$ and $\Z/n\Z$ (`fmpz_mat_mul`,
  `nmod_mat_mul`), and multimodular algorithms for integer matrices and
  polynomials (determinants, solving, characteristic polynomials, GCDs,
  Taylor shifts);
- multiplication and division of multivariate polynomials
  (`fmpz_mpoly_mul`, `fmpz_mpoly_divides`, `nmod_mpoly_mul`, ...);
- integer factoring with the quadratic sieve, Bernoulli numbers, Bell
  numbers and partition numbers;
- many functions of [](arb) and [](mp_real) that evaluate constants and
  elementary functions by binary splitting, and the generic DFTs of
  [](gr_dft).

Threads only pay off for large inputs; small operations run on one
thread regardless of the limit.

## Calling FLINT from several threads {#thread-safety}

FLINT is thread-safe in the usual sense: different threads may call
FLINT functions at the same time, provided that no object is written by
one thread while another thread reads or writes it. Objects that are only
read (inputs) may be shared between threads; outputs must be private to
one thread. Each thread needs its own random state (`flint_rand_t`).
Context objects are only read by most functions and may be shared, with
exceptions: for example the contexts of [](ca) cache data and must not
be shared. For the generic rings, `gr_ctx_is_threadsafe` tells whether a
context may be used from several threads.

This assumes the default build with thread-local storage
(`--enable-thread-safe`), in which the caches of FLINT are per thread or
are protected by locks. See [](guide/getting-started/building#threads)
for the other build modes. FLINT also calls MPFR in a few places (some
square roots and short products in the `arf` and `arb` modules, the
conversions to and from `mpfr_t`), so MPFR must be thread-safe too. It
is when it is built with thread-local storage, as it normally is
(`mpfr_buildopt_tls_p` returns nonzero).

Because most caches are per thread, each thread that needs a cached
value computes its own copy: if many threads evaluate $\pi$ to a million
digits, each computes it and keeps a copy in memory. A thread that has
used FLINT should free its caches with `flint_cleanup` before it exits,
or their memory is lost; the main thread calls `flint_cleanup_master` at
the end, after the other threads have finished (see
[](guide/conventions/memory#cleanup)). With OpenMP and other frameworks
that manage threads themselves, there is no simple way to run code when
a worker thread exits. Possible solutions are to arrange the parallel
work so that each thread calls `flint_cleanup` when it is done (for
example with OpenMP sections), or, in C++, to create a thread-local
object whose destructor calls `flint_cleanup`.

User threads and FLINT's thread pool combine: when several user threads
call parallelised FLINT functions, they share the workers of the global
pool, and each call gets as many workers as are free at that moment, up
to its limit.

## Writing threaded functions {#writing}

FLINT functions that use threads follow a common pattern, built on the
thread pool:

1. Ask for workers with `flint_request_threads`, which returns handles
   to at most `flint_get_num_threads() - 1` idle workers of the global
   pool (possibly none).
2. Divide the work into one part per worker plus one for the calling
   thread. Wake each worker with `thread_pool_wake`, passing a worker
   function and its argument; do the remaining part in the calling
   thread.
3. Wait for each worker with `thread_pool_wait`, and give the workers
   back with `flint_give_back_threads`.

The guide [](thread_pool/manual) shows this in detail. Simpler cases are
covered by `flint_parallel_do`, a parallel loop over independent
iterations, and `flint_parallel_binary_splitting`, which evaluates a
product tree in parallel; see [](thread_support/parallel_do) and
[](thread_support/bsplit).

The number of threads used in total must not exceed the limit set by the
user, also when threaded functions call each other. To this end,
"threads" and "workers" are distinguished: if `flint_get_num_threads()`
returns $n$, a function may start $n - 1$ workers in addition to the
thread it runs in. A woken worker by default may start no workers of its
own; the `max_workers` argument of `thread_pool_wake` allows it some.
To restrict the workers available to a called function temporarily, set
a lower limit with `flint_set_num_workers` and restore it afterwards
with `flint_reset_num_workers`.

The naming conventions are:

- User-facing functions do not have `_threaded` in their name. Either the
  function itself is threaded, or it calls a threaded implementation when
  that is the best strategy. Versions with `_threaded` in the name exist
  where it is useful to call the threaded algorithm directly, for testing
  or comparison, such as `nmod_mat_mul_classical_threaded` and
  `fmpz_mpoly_mul_heap_threaded`.
- Functions that take threads already obtained by the caller, instead of
  requesting them from the pool, take an array of thread pool handles and
  their number, and have names ending in `_threaded_pool`, such as
  `_fmpz_mod_poly_powers_mod_preinv_threaded_pool`. The usual distinction between
  underscore and non-underscore functions still applies.
