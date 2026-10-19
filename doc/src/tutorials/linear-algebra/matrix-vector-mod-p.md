---
title: Matrix–vector products modulo p
summary: Why reducing once per dot product instead of once per term pays off, how `nmod_vec` chooses the accumulator size, and how to parallelise over rows.
order: 2
---

Many algorithms over $\Z/p\Z$, from linear system solving to sparse
linear algebra and polynomial evaluation, spend most of their time in dot
products. This tutorial computes $y = Ax$ for an $N \times N$ matrix $A$
three times: with a modular reduction after every operation, with the dot
products of [](nmod_vec), and in parallel.

## Delaying the reduction {#delayed}

The obvious loop reduces after every multiplication and every addition:
$$
s \gets (s + (a_{ij} x_j \bmod p)) \bmod p.
$$
But the exact sum $\sum_j a_{ij} x_j$ of $N$ products of numbers below
$p$ is less than $N (p-1)^2$, which fits in at most three words. It is
therefore cheaper to accumulate the exact products in one, two or three
words and reduce once at the end [@DGP2008]. How many words are needed
depends on $N$ and $p$; `_nmod_vec_dot_params(N, mod)` works this out once
and selects a method, and `_nmod_vec_dot` then uses it for any vectors of
length at most $N$.

::: summary algorithms/linear-algebra/dot-products :::

For a 20-bit prime and $N = 400$, the sum is below $400 \cdot 2^{40} <
2^{49}$, so it fits in one word and the dot product is a plain loop of
multiply-adds, which is vectorised. For a 62-bit prime it needs
three words.

## The program {#program}

The program fills $A$ and $x$ with simple deterministic entries, computes
$y$ in the three ways, and checks that they agree. It also prints the
number of words `_nmod_vec_dot_bound_limbs` reports for the unreduced
sum.

```c example name=tutorial_matvec requires=64
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/thread_support.h>

#define N 400

/* Fill an N x N matrix and a vector with simple deterministic entries. */
static void fill(nn_ptr A, nn_ptr x, nmod_t mod)
{
    slong i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            A[i * N + j] = nmod_set_ui(i * i + 3 * j + 1, mod);
    for (j = 0; j < N; j++)
        x[j] = nmod_set_ui(j + 1, mod);
}

/* Version 1: reduce after every multiplication and addition. */
static void matvec_naive(nn_ptr y, nn_srcptr A, nn_srcptr x, nmod_t mod)
{
    slong i, j;
    for (i = 0; i < N; i++)
    {
        ulong s = 0;
        for (j = 0; j < N; j++)
            s = nmod_add(s, nmod_mul(A[i * N + j], x[j], mod), mod);
        y[i] = s;
    }
}

/* Version 2: one dot product per row, reducing once per row. */
static void matvec_dot(nn_ptr y, nn_srcptr A, nn_srcptr x, nmod_t mod)
{
    slong i;
    dot_params_t params = _nmod_vec_dot_params(N, mod);
    for (i = 0; i < N; i++)
        y[i] = _nmod_vec_dot(A + i * N, x, N, mod, params);
}

/* Version 3: the rows of version 2, distributed over threads. */
typedef struct
{
    nn_ptr y;
    nn_srcptr A;
    nn_srcptr x;
    nmod_t mod;
    dot_params_t params;
}
matvec_args_t;

static void matvec_row(slong i, void * vargs)
{
    matvec_args_t * args = vargs;
    args->y[i] = _nmod_vec_dot(args->A + i * N, args->x, N, args->mod, args->params);
}

static void matvec_parallel(nn_ptr y, nn_srcptr A, nn_srcptr x, nmod_t mod)
{
    matvec_args_t args;
    args.y = y;
    args.A = A;
    args.x = x;
    args.mod = mod;
    args.params = _nmod_vec_dot_params(N, mod);
    flint_parallel_do(matvec_row, &args, N, 0, FLINT_PARALLEL_UNIFORM);
}

static void run(ulong p)
{
    nmod_t mod;
    nn_ptr A = _nmod_vec_init(N * N), x = _nmod_vec_init(N);
    nn_ptr y1 = _nmod_vec_init(N), y2 = _nmod_vec_init(N), y3 = _nmod_vec_init(N);

    nmod_init(&mod, p);
    fill(A, x, mod);

    matvec_naive(y1, A, x, mod);
    matvec_dot(y2, A, x, mod);
    matvec_parallel(y3, A, x, mod);

    flint_printf("%2wu-bit modulus: dot products agree: %s, parallel agrees: %s, limbs needed: %d\n",
        (ulong) NMOD_BITS(mod),
        _nmod_vec_equal(y1, y2, N) ? "yes" : "no",
        _nmod_vec_equal(y1, y3, N) ? "yes" : "no",
        _nmod_vec_dot_bound_limbs(N, mod));

    _nmod_vec_clear(A);
    _nmod_vec_clear(x);
    _nmod_vec_clear(y1);
    _nmod_vec_clear(y2);
    _nmod_vec_clear(y3);
}

int main(void)
{
    flint_set_num_threads(4);

    run(n_nextprime(UWORD(1) << 20, 1));
    run(n_nextprime(UWORD(1) << 28, 1));
    run(n_nextprime(UWORD(1) << 40, 1));
    run(n_nextprime(UWORD(1) << 62, 1));

    flint_cleanup_master();
    return 0;
}
```
```output
21-bit modulus: dot products agree: yes, parallel agrees: yes, limbs needed: 1
29-bit modulus: dot products agree: yes, parallel agrees: yes, limbs needed: 2
41-bit modulus: dot products agree: yes, parallel agrees: yes, limbs needed: 2
63-bit modulus: dot products agree: yes, parallel agrees: yes, limbs needed: 3
```

The program needs a 64-bit build because of its 41- and 63-bit moduli.

## How much faster? {#timing}

On the machine used for these timings (x86-64, FLINT compiled for
AVX-512), with $N = 2000$ and one thread, the dot-product version was
about 5.5 times faster than the naive loop for the 21-bit prime, 4.3 times
for the 29-bit prime and 3 times for the 41- and 63-bit primes. The gain
is largest when the sum fits in one word, because then nothing but
multiplications and additions remains in the inner loop. Timings depend
on the processor; the `nmod_vec` profiling programs measure them on yours.

## Parallelism {#parallelism}

The rows are independent, so `flint_parallel_do` with one call per row
parallelises the product without any locking: call $i$ writes only
$y_i$. The dot-product parameters are computed once and passed to all
calls through the argument structure. For small $N$ the overhead of waking
threads dominates; parallelism pays off when each thread gets at least
some tens of microseconds of work.

For products of matrices rather than matrix–vector products, use
`nmod_mat_mul`, which chooses between classical multiplication with
delayed reduction, SIMD kernels, Strassen's algorithm, kernels for tiny
moduli and floating-point matrix multiplication depending on the size
and the modulus (see [](nmod_mat/mul#dispatch) and
[](algorithms/linear-algebra/matrix-multiplication)).
