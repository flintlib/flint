---
title: Memory management
summary: Initialising and clearing objects, FLINT's allocation functions, temporary allocation, the caches FLINT keeps and how to free them, and custom allocators.
order: 2
---

FLINT objects manage their own memory: every object is initialised
before use and cleared afterwards, and functions grow their outputs as
needed. Raw memory is allocated with FLINT's own versions of `malloc`
and `free`. FLINT also keeps some caches between calls, which can be
freed explicitly. The functions are documented in [](flint).

## Initialising and clearing objects {#init-clear}

Every variable of a FLINT type is initialised with `foo_init` (or a
variant such as `fmpz_init_set_ui`, `fmpz_poly_init2` with a given
allocation, or `nmod_poly_init(f, n)` with a modulus) before it is used,
and cleared with `foo_clear` exactly once when it is no longer needed.
Initialisation sets the object to zero and usually allocates nothing;
memory is allocated when a value is stored and is reused when the object
is overwritten. A cleared object may be initialised again.

Missing or duplicated `init` and `clear` calls are the most common cause
of memory leaks and crashes in programs using FLINT. A useful habit is to
write the `clear` calls together with the declarations, in the same
order. Context objects (`fmpz_mod_ctx_t`, `fq_ctx_t`, `gr_ctx_t`, ...) are
initialised and cleared in the same way, and must be cleared after all
objects that use them.

Vectors have their own functions: `_fmpz_vec_init(n)` allocates and
initialises $n$ integers and `_fmpz_vec_clear(v, n)` clears them, and
similarly `_arb_vec_init`, `_nmod_vec_init` and others. Arrays of objects
allocated with `flint_malloc` must have each entry initialised and
cleared individually.

The `fmpz` type deserves a remark. An `fmpz` holding an integer of at
most $\bits - 2$ bits stores it in the word itself and allocates nothing,
so most small integers cost no memory management at all. Larger values
use a GMP integer, which in the default build comes from a per-thread
pool: `fmpz_clear` returns it to the pool instead of freeing it.

## Allocation functions {#allocation}

FLINT allocates memory with `flint_malloc`, `flint_calloc`,
`flint_realloc` and `flint_free`, which have the same interface as the C
functions, and `flint_aligned_alloc` and `flint_aligned_free` for
aligned blocks. They differ from the C functions in one respect: they
never return `NULL`. When memory is exhausted they throw an exception,
which by default prints a message and aborts. Memory obtained from these
functions must be freed with `flint_free` (or `flint_aligned_free`), and
memory returned by FLINT functions that the caller must free, such as the
strings of `fmpz_get_str(NULL, 10, x)`, is freed with `flint_free` as
well.

GMP allocates memory through its own functions, which are `malloc`,
`realloc` and `free` by default and can be changed with GMP's
`mp_set_memory_functions`; MPFR uses GMP's functions. A large `fmpz`
integer is a GMP integer: the `mpz` structure is allocated by FLINT,
but its digits are allocated, resized and freed by GMP. FLINT does not
pass memory allocated by one set of functions to the other (see
[](#custom)).

## Temporary allocation {#temporary}

Inside a function, scratch space can be allocated with the `TMP` macros,
which take small blocks (at most 8192 bytes) from the stack with `alloca`
and larger ones from `flint_malloc`:

- `TMP_INIT;` with the declarations of the function;
- `TMP_START;` before the first allocation;
- `p = TMP_ALLOC(size)` or `p = TMP_ARRAY_ALLOC(n, type)` for each block;
- `TMP_END;` to free all blocks at once, before every return.

Stack space is only released when the function returns, so `TMP_ALLOC`
should not be called in a loop with many iterations or in a deeply
recursive function: each call consumes stack, and the stack is small
(typically 8 MB for the main thread and less for other threads). In such
cases, allocate once outside the loop, or use `flint_malloc`. In builds
with assertions enabled, all `TMP` blocks come from `flint_malloc`, so
that memory checkers can see them.

```c example name=guide_memory_tmp
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>

/* the sum of the squares of the nonzero entries of a, and their number */
static slong sum_squares(fmpz_t res, const fmpz * a, slong n)
{
    slong * idx;
    fmpz * t;
    slong i, k = 0;
    TMP_INIT;

    TMP_START;
    idx = TMP_ARRAY_ALLOC(n, slong);     /* plain scratch memory */
    for (i = 0; i < n; i++)
        if (!fmpz_is_zero(a + i))
            idx[k++] = i;

    t = _fmpz_vec_init(k);               /* fmpz entries need init and clear */
    for (i = 0; i < k; i++)
        fmpz_mul(t + i, a + idx[i], a + idx[i]);
    _fmpz_vec_sum(res, t, k);
    _fmpz_vec_clear(t, k);

    TMP_END;                             /* frees idx */
    return k;
}

int main(void)
{
    fmpz * a;
    fmpz_t s;
    slong i, n = 5, k;
    ulong * buf;

    a = _fmpz_vec_init(n);
    fmpz_init(s);
    for (i = 0; i < n; i++)
        fmpz_set_si(a + i, i - 2);

    k = sum_squares(s, a, n);
    flint_printf("%wd nonzero entries, sum of squares %{fmpz}\n", k, s);

    buf = flint_malloc(10 * sizeof(ulong));
    buf = flint_realloc(buf, 20 * sizeof(ulong));
    flint_free(buf);

    fmpz_clear(s);
    _fmpz_vec_clear(a, n);
    flint_cleanup_master();
    return 0;
}
```

```output
4 nonzero entries, sum of squares 10
```

## Caches and cleanup {#cleanup}

To save time, FLINT keeps some data after the functions that computed it
return. The main caches are:

- the pool of GMP integers that `fmpz` reuses (in the default build);
- tables of primes and of their inverses (`n_primes_arr_readonly`), and
  the product tree used for trial division;
- constants of ball arithmetic such as $\pi$, $\log 2$ at precisions
  beyond the built-in tables, the logarithms of
  small primes used by the elementary functions, Bernoulli numbers and
  quadrature nodes (see [](guide/types/ball-arithmetic#caches));
- precomputed tables of [](fft_small), of radix conversion, and of some
  matrix multiplication kernels;
- the constants and tables of [](mp_real), which are shared by all
  threads (see [](thread_support));
- the threads of the global thread pool.

This memory is not leaked, since it remains reachable and is reused,
but tools such as Valgrind report it as "still reachable".

Most caches belong to the thread that created them (when FLINT is built
with thread-local storage, the default). `flint_cleanup()` frees the
caches of the calling thread: it calls the cleanup functions registered
in this thread, returns the `fmpz` pool, and calls MPFR's
`mpfr_free_cache`. A thread that used FLINT should call it before it
exits. The shared caches of [](thread_support) are not freed by
`flint_cleanup`, since other threads may be using them.
`flint_cleanup_master()` frees everything: the thread pool, the shared
caches and the caches of the calling thread. Call it once, from the main
thread, when no other thread uses FLINT any more, typically at the end
of `main`. Neither function is needed for correctness. In a build with
`--enable-reentrant` and without thread-local storage, `flint_cleanup`
does nothing.

Code that builds its own caches can register a function `void f(void)`
with `flint_register_cleanup_function`; it is then called by the next
`flint_cleanup` or `flint_cleanup_master` in the same thread.

## Custom allocators {#custom}

The functions behind `flint_malloc`, `flint_calloc`, `flint_realloc` and
`flint_free` can be replaced with
`__flint_set_memory_functions(alloc, calloc, realloc, free)`, for example
to count allocations or to use a special allocator, and retrieved with
`__flint_get_memory_functions`. After this call `flint_aligned_alloc`
and `flint_aligned_free` over-allocate with the new `alloc` and `free`;
`__flint_set_all_memory_functions` takes two more functions and replaces
the aligned allocation functions as well.

The replacement functions must behave like their C counterparts, with a
few precise requirements:

- They must be installed before FLINT allocates any memory, in practice
  at the start of the program before any other FLINT function is
  called. FLINT keeps memory in its caches across calls (the `fmpz`
  pool, prime tables, cached constants), and a block allocated by the
  old functions would later be freed or resized by the new ones.
- `alloc`, `calloc` and `realloc` return `NULL` on failure; `flint_malloc`
  and its relatives then throw. The memory they return must be aligned
  for any type, as with `malloc`.
- `free` must accept `NULL`: `flint_free` passes its argument on
  unchanged, and FLINT does call it with `NULL`. `realloc` is never
  called with `NULL`, since `flint_realloc(NULL, size)` calls `alloc`.
- If FLINT is used from several threads, or uses its own worker threads
  (`flint_set_num_threads`), the functions must be thread-safe.

FLINT does not require its allocation functions to be compatible with
GMP's. Memory allocated by GMP (the digits of large `fmpz` integers and
of the GMP and MPFR numbers used internally) is resized and freed only by
GMP, and memory allocated by FLINT only by FLINT; FLINT never calls
`mp_set_memory_functions` itself, also not in a garbage-collected build.
To route all memory through custom functions, set GMP's functions too,
with `mp_set_memory_functions`, under the same condition: before the
first FLINT or GMP call. This matters especially for FLINT, because the
`fmpz` pool keeps GMP integers, with their GMP-allocated digits, after
`fmpz_clear`, until `flint_cleanup` returns them to GMP.

FLINT can also be built to use the Boehm–Demers–Weiser garbage collector
(`--with-gc`), in which case `flint_free` does nothing and memory is
reclaimed by the collector. This is mainly useful for language bindings
with garbage collection.
