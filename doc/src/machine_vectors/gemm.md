---
title: Floating-point matrix multiplication
functions: flint_sgemm flint_sgemm_fallback flint_sgemm_blas flint_gemm_use_blas
---

`flint_sgemm` and `flint_dgemm` multiply matrices of floats and doubles.
They exist so that FLINT modules can use a fast matrix product whether or
not FLINT was built with BLAS. FLINT uses them for exact computations,
in which every entry of the result is an integer below $2^{24}$ or
$2^{53}$ in absolute value, so the floating-point product has no rounding
error.

## Interface {#interface}

The functions compute $C = AB$ where $A$ is $m \times k$, $B$ is
$k \times n$ and $C$ is $m \times n$, all stored row by row. Entry
$(i, j)$ of $A$ is `A[i * lda + j]`, so a submatrix is passed by a
pointer to its first entry and the row stride of the matrix that contains
it. This is the subset of the BLAS routine `xGEMM` that FLINT needs:

- the dimensions come in the order $m, k, n$, whereas `cblas_dgemm`
  takes $m, n, k$;
- there is no transposition and no scaling: `flint_dgemm(m, k, n, A, lda,
  B, ldb, C, ldc)` is `cblas_dgemm(CblasRowMajor, CblasNoTrans,
  CblasNoTrans, m, n, k, 1.0, A, lda, B, ldb, 0.0, C, ldc)`;
- $C$ is overwritten, never read, and must not overlap $A$ or $B$;
- if $m \le 0$ or $n \le 0$ nothing is done, and if $k = 0$ the
  $m \times n$ block of $C$ is set to zero. A negative $k$ is treated as
  0 by the FLINT kernels but is an invalid argument for the BLAS.

## Choosing an implementation {#implementation}

In a BLAS build `flint_sgemm` and `flint_dgemm` call the BLAS by
default. Setting `flint_gemm_use_blas` to 0 makes them use FLINT's own
kernels instead, which is useful for benchmarks and for avoiding
interference between the thread pools of FLINT and of the BLAS. The
variable is not protected by a lock, so it should only be changed when
no other thread is multiplying.

The FLINT kernels use up to `flint_get_num_threads()` threads from
FLINT's thread pool. A caller that holds worker threads from the pool
should give them back before calling, or the kernels find no free
threads and run serially; `nmod_mat_mul_blas` does this around its calls.

## Accuracy {#accuracy}

The FLINT kernels evaluate each entry as a sum of products accumulated
with fused multiply-add instructions where the backend has them. The
order of summation depends on the shape and on the blocking, so results
for non-integer data may differ in the last bits from those of a BLAS or
from a naive loop. When all partial sums are integers of magnitude at
most $2^{53}$ (doubles) or $2^{24}$ (floats), every operation is exact and
all implementations give the same result. For integer matrices this
holds when

$$
k \cdot \max_{i,l} |A_{il}| \cdot \max_{l,j} |B_{lj}| \le 2^{53}
$$

for `flint_dgemm`, and with $2^{24}$ for `flint_sgemm`.

::: summary algorithms/linear-algebra/gemm :::

## Examples {#examples}

A product of small matrices. The second call multiplies the leading
$2 \times 2$ block of $A$ by the leading $2 \times 2$ block of $B$; only
the leading dimension 3 tells the function where the rows of $A$ start.

```c example name=mv_gemm_basic
#include <flint/flint.h>
#include <flint/machine_vectors.h>

int main(void)
{
    /* A is 2 x 3 and B is 3 x 2, both stored row by row */
    double A[6] = { 1, 2, 3,
                    4, 5, 6 };
    double B[6] = {  7,  8,
                     9, 10,
                    11, 12 };
    double C[4];
    slong i, j;

    flint_dgemm(2, 3, 2, A, 3, B, 2, C, 2);

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
            flint_printf("%5wd", (slong) C[i * 2 + j]);
        flint_printf("\n");
    }

    /* the top left 2 x 2 blocks: only the leading dimensions change */
    flint_dgemm(2, 2, 2, A, 3, B, 2, C, 2);
    flint_printf("%5wd%5wd\n%5wd%5wd\n", (slong) C[0], (slong) C[1],
                 (slong) C[2], (slong) C[3]);
    return 0;
}
```
```output
   58   64
  139  154
   25   28
   73   82
```

A product modulo $p = 1000003$ with one `flint_dgemm`. The entries are
reduced, so every entry of the integer product is at most
$k (p - 1)^2 = 200 \cdot 1000002^2 < 2^{48}$, and the floating-point
result is exact. The program checks it against `nmod_mat_mul`.

```c example name=mv_gemm_mod
#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/machine_vectors.h>

int main(void)
{
    slong m = 50, k = 200, n = 60, i, j;
    ulong p = 1000003;
    nmod_mat_t A, B, C;
    double * dA = flint_malloc(m * k * sizeof(double));
    double * dB = flint_malloc(k * n * sizeof(double));
    double * dC = flint_malloc(m * n * sizeof(double));
    int ok = 1;

    nmod_mat_init(A, m, k, p);
    nmod_mat_init(B, k, n, p);
    nmod_mat_init(C, m, n, p);

    for (i = 0; i < m; i++)
        for (j = 0; j < k; j++)
            nmod_mat_entry(A, i, j) = (UWORD(7919) * (i * k + j) + 13) % p;
    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
            nmod_mat_entry(B, i, j) = (UWORD(104729) * (i * n + j) + 5) % p;

    /* convert to doubles, multiply, reduce */
    for (i = 0; i < m * k; i++)
        dA[i] = (double) nmod_mat_entry(A, i / k, i % k);
    for (i = 0; i < k * n; i++)
        dB[i] = (double) nmod_mat_entry(B, i / n, i % n);

    flint_dgemm(m, k, n, dA, k, dB, n, dC, n);

    nmod_mat_mul(C, A, B);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            if ((ulong) dC[i * n + j] % p != nmod_mat_entry(C, i, j))
                ok = 0;

    flint_printf("C[0][0] = %wu, unreduced %wu\n", nmod_mat_entry(C, 0, 0),
                 (ulong) dC[0]);
    flint_printf("all entries agree: %s\n", ok ? "yes" : "no");

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);
    flint_free(dA);
    flint_free(dB);
    flint_free(dC);
    return 0;
}
```
```output
C[0][0] = 683309, unreduced 42685113738264
all entries agree: yes
```

`nmod_mat_mul_blas` refines this in two ways. It lifts entries to the
symmetric range $(-p/2, p/2]$, which allows four times larger $k$ for the
same modulus, and for larger moduli it multiplies modulo several smaller
coprime numbers and combines the results by Chinese remaindering.
