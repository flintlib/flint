---
title: Binary splitting on several threads
functions: flint_parallel_binary_splitting bsplit_basecase_func_t FLINT_PARALLEL_BSPLIT_LEFT_INPLACE
---

Binary splitting computes a product $f(a) \circ f(a+1) \circ \cdots \circ
f(b-1)$ of many factors by splitting $[a, b)$ in the middle, computing both
halves recursively and combining them [@HP1998; @Ber2008]. When the
operation is multiplication of integers, polynomials or matrices whose size
grows with the number of factors, this balances the operand sizes and lets
fast multiplication do most of the work. It is the standard method for
factorials, binomial coefficients, and the evaluation of hypergeometric
series.

`flint_parallel_binary_splitting` runs the two halves of each split on
different threads, as long as threads are available, so the recursion tree
is processed in parallel near its root and serially near its leaves.

## The callbacks {#callbacks}

The results of subranges are opaque objects of `sizeof_res` bytes, which
the function allocates itself and manages with the callbacks:

- `init(x, args)` and `clear(x, args)` create and destroy a temporary;
- `basecase(res, a, b, args)` computes the result for $[a, b)$ directly,
  and is used when $b - a$ is at most `basecase_cutoff` (a cutoff below 1
  is treated as 1);
- `merge(res, x, y, args)` sets `res` to $x \circ y$. The operation must
  be associative; it need not be commutative, and $x$ always belongs to
  the left subrange.

With `FLINT_PARALLEL_BSPLIT_LEFT_INPLACE`, the left result is computed
directly in `res` and `merge` is called with `res` equal to `x`, which
saves one temporary at each level of the recursion. Use it when the merge
function supports this aliasing, as `fmpz_mul` does.

## Threads {#threads}

Let $n$ be `thread_limit`, or `flint_get_num_threads()` if `thread_limit`
is not positive. At each split with $n \ge 2$, the function asks the
global pool for one worker. If it gets one, the right half runs on the
worker and the left half on the calling thread, and the budget is divided
between them: the worker may start $\lfloor n/2 \rfloor - 1$ further
workers and the calling thread at most $n - \lfloor n/2 \rfloor - 1$.
Otherwise both halves run on the calling thread, one after the other.
The recursive calls receive these shares as their `thread_limit`, so at
most $n$ threads, the calling thread included, run base cases at the same
time.

::: summary algorithms/parallel/parallel-loops :::

## Example {#example}

This computes $20000!$ as a product of the integers in $[1, 20001)$ and
checks it against `fmpz_fac_ui`.

```c example name=bsplit_factorial
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/thread_support.h>

/* The product a * (a + 1) * ... * (b - 1), by binary splitting. */
static void basecase(void * res, slong a, slong b, void * args)
{
    slong k;
    fmpz_one(res);
    for (k = a; k < b; k++)
        fmpz_mul_ui(res, res, k);
}

static void merge(void * res, void * x, void * y, void * args)
{
    fmpz_mul(res, x, y);
}

static void init(void * x, void * args) { fmpz_init(x); }
static void clear(void * x, void * args) { fmpz_clear(x); }

int main(void)
{
    fmpz_t f, g;

    flint_set_num_threads(4);
    fmpz_init(f);
    fmpz_init(g);

    /* 20000! = 1 * 2 * ... * 20000; the leaves multiply 64 numbers each */
    flint_parallel_binary_splitting(f, basecase, merge, sizeof(fmpz),
        init, clear, NULL, 1, 20001, 64, -1, 0);

    fmpz_fac_ui(g, 20000);
    flint_printf("20000! has %wu digits\n", (ulong) fmpz_sizeinbase(f, 10));
    flint_printf("agrees with fmpz_fac_ui: %s\n", fmpz_equal(f, g) ? "yes" : "no");

    fmpz_clear(f);
    fmpz_clear(g);
    flint_cleanup_master();
    return 0;
}
```
```output
20000! has 77338 digits
agrees with fmpz_fac_ui: yes
```

`fmpz_fac_ui` is faster for this particular product; the example only
shows the mechanics.
