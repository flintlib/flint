---
title: Multiplying matrices
functions: nmod_mat_mul nmod_mat_addmul
---

`nmod_mat_mul` computes $C = AB$ for an $m \times k$ matrix $A$ and a
$k \times n$ matrix $B$ over $\Z/N\Z$ (here $N$ is the modulus,
since $n$ is a dimension). It is the function to call: it picks one of
about ten methods from the dimensions, the size of $N$, the number of
threads and the instruction set, and it accepts any modulus, prime or
not. `nmod_mat_addmul` and `nmod_mat_submul` add the product to, or
subtract it from, another matrix. Almost all of the linear algebra in
this module, from LU decomposition to triangular solving, reduces to
these three functions.

::: summary algorithms/linear-algebra/matrix-multiplication :::

## Usage {#usage}

- $C$ must be initialised as an $m \times n$ matrix with the same modulus
  as $A$ and $B$. The dimensions are checked only by `FLINT_ASSERT`.
- `nmod_mat_mul` allows $C$ to be the same matrix as $A$ or $B$, for
  example `nmod_mat_mul(A, A, A)` to square in place; the product is then
  computed into a temporary matrix. Windows are allowed, but $C$ must not
  partially overlap $A$ or $B$.
- `nmod_mat_addmul(D, C, A, B)` and `nmod_mat_submul(D, C, A, B)` allow
  $D$ to be $C$, which is the usual way to update a matrix in place, but
  not $A$ or $B$.
- With `flint_set_num_threads`, the large products and most kernels use
  several threads.

All methods compute the exact product modulo $N$, so the result never
depends on the method chosen.

## How the method is chosen {#dispatch}

On 64-bit machines `nmod_mat_mul` tests the following cases in order,
where $d = \min(m, k, n)$.

| Case | Method |
|---|---|
| $N \le 255$ and $d \ge 8$ | `nmod_mat_mul_u8`: bit-packed for $N = 2, 3$, byte table lookups for $N \le 16$, `flint_sgemm` above |
| $n = 1$, $m \ge 8$, one thread | `nmod_mat_mul_nmod_vec` |
| at most 4 rows in $A$ (8 with AVX512-IFMA), $k \ge 2$, a vector code path for $N$, $C$ not aliased | the rows of $B$ are streamed once for all rows of $A$, as in `nmod_mat_nmod_vec_mul` |
| 2 to 4 columns in $B$ (8 with AVX512-IFMA), one thread, a vector code path for $N$, $C$ not aliased | each row of $A$ is multiplied by all columns of $B$ at once |
| $N \le 2^{52}$, $n \ge 8$ and $m \ge 8$ (or $m \ge 4$ and $n \ge 32$) | a SIMD kernel, see below, with Strassen–Winograd on top for $d \ge$ `FLINT_NMOD_MAT_MUL_SIMD_STRASSEN_CUTOFF` and one thread |
| $d > 100$ and $d$ above a cutoff | `nmod_mat_mul_blas` |
| several threads | `nmod_mat_mul_classical_threaded` |
| $d < 400$ ($N < 2048$) or $d < 200$ | `nmod_mat_mul_classical` |
| otherwise | `nmod_mat_mul_strassen` |

The SIMD kernel depends on the bit size $b$ of $N$:

- $b \le 32$: `nmod_mat_mul_u32`, or `nmod_mat_mul_fp50` for a tiny inner
  dimension (by default $k \le 3$ with AVX2 and $k = 1$ otherwise), or
  `nmod_mat_mul_blas` when a single floating-point product is exact
  ($k \lfloor N/2 \rfloor^2 < 2^{53}$) and $d$ is at least
  `FLINT_NMOD_MAT_MUL_BLAS_1PASS_CUTOFF` (256 by default on x86-64,
  between 128 and 800 on the tuned x86-64 targets, and 0, meaning never,
  on Zen 4, Tiger Lake, Sapphire Rapids and AArch64).
- $33 \le b \le 50$: `nmod_mat_mul_fp50` on x86-64,
  `nmod_mat_mul_k52` on AArch64.
- $51 \le b \le 52$: `nmod_mat_mul_k52`.
- With AVX512-IFMA, `nmod_mat_mul_u52` replaces these for $b \ge 31$ (and
  on some processors for $b \le 26$) when $k \ge 12$, and for $b \ge 51$
  whatever $k$.
- With an external BLAS, `nmod_mat_mul_blas` replaces `nmod_mat_mul_k52`
  and `nmod_mat_mul_fp50` from dimension
  `FLINT_NMOD_MAT_MUL_K52_BLAS_CUTOFF` on, where this parameter is nonzero
  (384 on Apple M processors, 512 on Broadwell, 3000 on Skylake).

For $N > 2^{52}$ and large matrices, `nmod_mat_mul_blas` reduces the
product to several floating-point products and a Chinese remainder
reconstruction. It is used from $d > 450$ with one thread, where it
competes with Strassen–Winograd, and from
$d > 100 + 5 T \max(b, 32)/2$ with $T$ threads. Below, the classical
algorithm forms each entry as a dot product with a single reduction
([](algorithms/linear-algebra/dot-products)), and Strassen–Winograd is
used from $d = 200$.

The cutoffs are the `FLINT_NMOD_MAT_MUL_*` parameters of the
`flint-mparam.h` file selected for the build machine, measured with
`src/nmod_mat/profile/p-mul_tune.c`. On 32-bit machines only
`nmod_mat_mul_u8` (without bit packing), the classical algorithm and
Strassen–Winograd are available.

## Adding a product {#addmul}

`nmod_mat_addmul` and `nmod_mat_submul` compute $D = C \pm AB$. They
compute $AB$ with `nmod_mat_mul` into a temporary matrix and add it,
except with one thread, a dimension below 400 ($N < 2048$) or 200, and a
modulus for which `nmod_mat_mul` has no faster kernel than the classical
algorithm. In that case they use the classical algorithm and add each dot
product to the entry of $C$ before the final reduction.

On 64-bit machines the kernels of `nmod_mat_mul` are preferred for
$N \le 2^{32}$, and for $N \le 2^{52}$ when the AVX512-IFMA or the
floating-point kernels are available (on 32-bit machines only for
$N \le 255$). The fused classical branch is kept for small products
(fewer than 4 columns of $B$, or $mk \min(n, 128) < 128$), for which the
temporary matrix and the extra pass cost as much as they save, and for
$N \le 255$ when $\min(m, k)$ and $n$ are at least 8 but one of them is
below 16, where `nmod_mat_mul_u8` does not pay for its conversions.

## Examples {#examples}

Products of $120 \times 90$ and $90 \times 110$ matrices for moduli of
several sizes, which go to different methods, checked against the
classical algorithm (64-bit only, since some moduli exceed $2^{32}$):

```c example name=nmod_mat_mul_moduli
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>

/* a simple deterministic pseudo-random generator (64-bit) */
static ulong next(ulong * x)
{
    *x = *x * UWORD(6364136223846793005) + UWORD(1442695040888963407);
    return *x >> 1;
}

int main(void)
{
    /* moduli of 7, 31, 45, 52 and 64 bits */
    ulong ns[5] = { 101, 2147483647, UWORD(35184372088777),
                    UWORD(4503599627370449), UWORD(18446744073709551557) };
    slong m = 120, k = 90, n = 110, i, j, t;
    ulong x = 1;

    for (t = 0; t < 5; t++)
    {
        nmod_mat_t A, B, C, D;

        nmod_mat_init(A, m, k, ns[t]);
        nmod_mat_init(B, k, n, ns[t]);
        nmod_mat_init(C, m, n, ns[t]);
        nmod_mat_init(D, m, n, ns[t]);

        for (i = 0; i < m; i++)
            for (j = 0; j < k; j++)
                nmod_mat_entry(A, i, j) = next(&x) % ns[t];
        for (i = 0; i < k; i++)
            for (j = 0; j < n; j++)
                nmod_mat_entry(B, i, j) = next(&x) % ns[t];

        nmod_mat_mul(C, A, B);              /* chooses the method */
        nmod_mat_mul_classical(D, A, B);    /* reference */

        flint_printf("%2wu bits: C[5][7] = %wu, agrees with classical: %d\n",
            NMOD_BITS(nmod_mat_mod(C)), nmod_mat_entry(C, 5, 7),
            nmod_mat_equal(C, D));

        nmod_mat_clear(A);
        nmod_mat_clear(B);
        nmod_mat_clear(C);
        nmod_mat_clear(D);
    }
    return 0;
}
```
```output
 7 bits: C[5][7] = 90, agrees with classical: 1
31 bits: C[5][7] = 696258119, agrees with classical: 1
45 bits: C[5][7] = 21009516881923, agrees with classical: 1
52 bits: C[5][7] = 3687877138422982, agrees with classical: 1
64 bits: C[5][7] = 12025067871927492558, agrees with classical: 1
```

Updating a matrix in place with `nmod_mat_submul` and
`nmod_mat_addmul`, and squaring in place with `nmod_mat_mul`:

```c example name=nmod_mat_addmul_inplace
#include <flint/flint.h>
#include <flint/nmod_mat.h>

int main(void)
{
    nmod_mat_t A, B, S;
    slong i, j;

    nmod_mat_init(A, 3, 3, 97);
    nmod_mat_init(B, 3, 3, 97);
    nmod_mat_init(S, 3, 3, 97);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            nmod_mat_entry(A, i, j) = (i + 2 * j + 1) % 97;
            nmod_mat_entry(B, i, j) = (i == j) ? 2 : 1;
        }

    /* S = A*B, then S = S - A*B in place: the result is zero */
    nmod_mat_mul(S, A, B);
    nmod_mat_submul(S, S, A, B);
    flint_printf("A*B - A*B is zero: %d\n", nmod_mat_is_zero(S));

    /* S = B + A*A, then square A in place */
    nmod_mat_addmul(S, B, A, A);
    nmod_mat_mul(A, A, A);
    nmod_mat_add(A, A, B);
    flint_printf("addmul agrees: %d\n", nmod_mat_equal(S, A));
    nmod_mat_print_pretty(S);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(S);
    return 0;
}
```
```output
A*B - A*B is zero: 1
addmul agrees: 1
<3 x 3 integer matrix mod 97>
[24 41 59]
[29 54 77]
[35 65 96]
```
