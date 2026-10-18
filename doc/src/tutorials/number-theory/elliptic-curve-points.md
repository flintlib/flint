---
title: Counting points on an elliptic curve
summary: From a brute-force double loop to a parallel $O(p)$ count, using `nmod` arithmetic, Legendre symbols, a table of squares and `flint_parallel_do`.
order: 1
---

This tutorial counts the points of the elliptic curve
$$
E: y^2 = x^3 + ax + b
$$
over the field $\mathbb{F}_p$ for a prime $p$, including the point at
infinity. The number of points $N$ is the most basic invariant of the
curve, and computing it is the first step in many applications. We start
with the obvious method and improve it in three steps, each of which uses
a different part of FLINT:

1. a brute-force count with the modular arithmetic of [](nmod);
2. a count with Legendre symbols from [](ulong_extras), which is
   asymptotically faster;
3. a table of squares and a parallel loop from [](thread_support).

Throughout we use the curve $a = 2$, $b = 3$, which is nonsingular for all
the primes below because its discriminant $-16(4a^3 + 27b^2) = -16 \cdot
275 = -2^4 \cdot 5^2 \cdot 11$ is not divisible by them.

## Step 1: brute force {#brute-force}

For each $x \in \mathbb{F}_p$ we compute $r = x^3 + ax + b$ and count the
$y$ with $y^2 = r$. The arithmetic modulo $p$ uses an `nmod_t`, which holds
$p$ with a precomputed inverse so that `nmod_mul` needs no division
instruction (see [](nmod/arithmetic)).

```c example name=tutorial_ec_naive
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>

/* Number of points on y^2 = x^3 + a x + b over GF(p), p an odd prime,
   including the point at infinity, by brute force: for each x, count
   the y with y^2 = x^3 + a x + b. */
ulong count_points_naive(ulong a, ulong b, ulong p)
{
    nmod_t mod;
    ulong x, y, rhs, count = 1;     /* the point at infinity */

    nmod_init(&mod, p);

    for (x = 0; x < p; x++)
    {
        rhs = nmod_add(nmod_mul(nmod_mul(x, x, mod), x, mod),
                       nmod_add(nmod_mul(a, x, mod), b, mod), mod);
        for (y = 0; y < p; y++)
            if (nmod_mul(y, y, mod) == rhs)
                count++;
    }

    return count;
}

int main(void)
{
    flint_printf("p = 101: %wu points\n", count_points_naive(2, 3, 101));
    return 0;
}
```
```output
p = 101: 96 points
```

This takes $p^2$ multiplications, which is already slow for $p$ around
$10^5$.

## Step 2: Legendre symbols {#legendre}

The number of $y$ with $y^2 = r$ is 1 if $r = 0$, and otherwise 2 or 0
depending on whether $r$ is a square modulo $p$. In terms of the Legendre
symbol $(r/p)$, which is 1 for nonzero squares, $-1$ for nonsquares and 0
for $r = 0$, the count is $1 + (r/p)$, so
$$
N = 1 + \sum_{x \in \mathbb{F}_p} \left(1 + \left(\frac{x^3 + ax + b}{p}\right)\right)
  = p + 1 + \sum_{x \in \mathbb{F}_p} \left(\frac{x^3 + ax + b}{p}\right).
$$
For prime $p$, the Legendre symbol equals the Jacobi symbol, which
`n_jacobi_unsigned` computes with a binary algorithm that uses only
shifts, subtractions and bit operations [@ShaSor1993]. We also evaluate $r$ by Horner's rule.

The difference $t = p + 1 - N$ is the *trace of Frobenius*. By Hasse's
theorem [@Sil2009, Theorem V.1.1], $\lvert t \rvert \le 2\sqrt p$, which
gives a useful check.

```c example name=tutorial_ec_jacobi
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>

/* The same count with Legendre symbols: x contributes 1 + (r/p) points,
   where r = x^3 + a x + b. Cost O(p log p) instead of O(p^2). */
ulong count_points(ulong a, ulong b, ulong p)
{
    nmod_t mod;
    ulong x, r;
    slong sum = 0;

    nmod_init(&mod, p);

    for (x = 0; x < p; x++)
    {
        /* Horner: r = (x^2 + a) x + b */
        r = nmod_add(nmod_mul(nmod_add(nmod_mul(x, x, mod), a, mod), x, mod), b, mod);
        sum += n_jacobi_unsigned(r, p);
    }

    return p + 1 + sum;
}

int main(void)
{
    ulong p, a = 2, b = 3, N;
    slong t;
    int i;

    flint_printf("p = 101: %wu points\n", count_points(a, b, 101));

    /* Hasse's theorem: |p + 1 - N| <= 2 sqrt(p) */
    p = 1000000;
    for (i = 0; i < 5; i++)
    {
        p = n_nextprime(p, 1);
        N = count_points(a, b, p);
        t = (slong) (p + 1) - (slong) N;
        flint_printf("p = %wu: N = %wu, trace t = %wd, t^2 <= 4p: %s\n",
                     p, N, t, (ulong) (t * t) <= 4 * p ? "yes" : "no");
    }

    return 0;
}
```
```output
p = 101: 96 points
p = 1000003: N = 999708, trace t = 296, t^2 <= 4p: yes
p = 1000033: N = 1000284, trace t = -250, t^2 <= 4p: yes
p = 1000037: N = 999316, trace t = 722, t^2 <= 4p: yes
p = 1000039: N = 1000564, trace t = -524, t^2 <= 4p: yes
p = 1000081: N = 1000684, trace t = -602, t^2 <= 4p: yes
```

The count for $p = 101$ agrees with step 1. The cost is now about
$p \log p$ bit operations.

## Step 3: a table of squares, in parallel {#parallel}

The Jacobi symbol is the most expensive part of the loop. When there is
enough memory for a table of $p$ bytes, we can instead mark the nonzero
squares $y^2$ for $1 \le y \le (p-1)/2$ (every nonzero square arises from
exactly two values $\pm y$) and look them up. That makes the count
$O(p)$.

The loop over $x$ is now a sum of independent terms, which
`flint_parallel_do` can distribute over threads (see
[](thread_support/parallel_do)). We split $\mathbb{F}_p$ into 64 blocks;
each call of the loop body sums one block and writes its result to its own
slot of an array, so no locking is needed. The table is filled on the
calling thread and only read by the workers.

```c example name=tutorial_ec_parallel
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/thread_support.h>

typedef struct
{
    nmod_t mod;
    ulong a, b;
    const unsigned char * is_square;    /* is_square[r] = 1 iff r is a nonzero square */
    slong * sums;                       /* one partial sum per block */
    slong blocks;
}
count_args_t;

/* Sum of (r/p) over the x in block i, using the table of squares. */
static void count_block(slong i, void * vargs)
{
    count_args_t * args = vargs;
    nmod_t mod = args->mod;
    ulong p = mod.n, x, r;
    ulong start = (p * i) / args->blocks, stop = (p * (i + 1)) / args->blocks;
    slong sum = 0;

    for (x = start; x < stop; x++)
    {
        r = nmod_add(nmod_mul(nmod_add(nmod_mul(x, x, mod), args->a, mod), x, mod), args->b, mod);
        if (r != 0)
            sum += args->is_square[r] ? 1 : -1;
    }

    args->sums[i] = sum;
}

ulong count_points_parallel(ulong a, ulong b, ulong p)
{
    count_args_t args;
    unsigned char * is_square;
    slong i, sum = 0;
    ulong y;

    nmod_init(&args.mod, p);
    args.a = a;
    args.b = b;
    args.blocks = 64;
    args.sums = flint_malloc(args.blocks * sizeof(slong));

    /* The nonzero squares are y^2 for 1 <= y <= (p - 1)/2. */
    is_square = flint_calloc(p, 1);
    for (y = 1; y <= (p - 1) / 2; y++)
        is_square[nmod_mul(y, y, args.mod)] = 1;
    args.is_square = is_square;

    flint_parallel_do(count_block, &args, args.blocks, 0, FLINT_PARALLEL_UNIFORM);

    for (i = 0; i < args.blocks; i++)
        sum += args.sums[i];

    flint_free(is_square);
    flint_free(args.sums);
    return p + 1 + sum;
}

int main(void)
{
    ulong p = n_nextprime(10000000, 1), N;

    flint_set_num_threads(4);

    N = count_points_parallel(2, 3, 101);
    flint_printf("p = 101: %wu points\n", N);

    N = count_points_parallel(2, 3, p);
    flint_printf("p = %wu: %wu points\n", p, N);

    flint_cleanup_master();
    return 0;
}
```
```output
p = 101: 96 points
p = 10000019: 9997980 points
```

`flint_set_num_threads(4)` allows FLINT to use four threads;
`flint_parallel_do` with a thread limit of 0 then uses as many of them as
it can get. The result does not depend on the number of threads.

## Going further {#further}

The methods above take time linear in $p$. Point counting for
cryptographic sizes, $p \approx 2^{256}$, uses Schoof's algorithm and its
improvements by Elkies and Atkin, which take time polynomial in
$\log p$ [@Sch1985; @Sch1995]. For moderate $p$, baby-step giant-step
methods based on Mestre's theorem need about $p^{1/4}$ operations
[@Coh1996].
