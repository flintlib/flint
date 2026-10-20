---
title: Combining modules: a multimodular determinant
summary: Integer matrices, matrices modulo word-size primes and the Chinese remainder theorem, combined to compute a determinant over the integers from determinants modulo primes.
order: 2
---

Many algorithms in FLINT reduce a problem over the integers to the same
problem modulo several primes, where all numbers fit in a machine word,
and reconstruct the integer answer with the Chinese remainder theorem.
This tutorial writes such an algorithm for the determinant of an integer
matrix, to show how the modules [](fmpz_mat), [](nmod_mat),
[](ulong_extras) and [](fmpz) work together. For large matrices,
FLINT's own `fmpz_mat_det` does the same, with refinements, as described
in [](algorithms/linear-algebra/multimodular).

## The method {#method}

Let $A$ be an $n \times n$ integer matrix. Hadamard's inequality bounds
the determinant by the product of the Euclidean norms of the rows,
$|\det A| \le \prod_i \lVert a_i \rVert_2$, and `fmpz_mat_det_bound`
computes this bound $H$. If $p_1, \ldots, p_k$ are primes with
$m = p_1 \cdots p_k > 2H$, then $\det A$ is the unique integer in
$(-m/2, m/2]$ that is congruent to $\det A \bmod p_i$ for every $i$. The
determinants modulo $p_i$ are computed with `nmod_mat_det` by Gaussian
elimination over $\Z/p_i\Z$, with word-size arithmetic. The Chinese
remainder step `fmpz_CRT_ui(det, det, m, r, p, 1)` combines the
current value modulo $m$ with the residue $r$ modulo $p$ into the value
modulo $mp$, taken in the symmetric range when the last argument is 1.

Choosing primes close to $2^{\bits - 2}$ (with `n_nextprime`) makes each
prime contribute about 62 bits on a 64-bit machine, while keeping the
modular arithmetic fast.

## The program {#program}

```c example name=tutorial_multimodular_det
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>

#define N 8

int main(void)
{
    fmpz_mat_t A;
    nmod_mat_t Amod;
    fmpz_t det, bound, modulus, d;
    ulong p, r;
    slong i, j, nprimes = 0;

    fmpz_mat_init(A, N, N);
    fmpz_init(det);
    fmpz_init(bound);
    fmpz_init(modulus);
    fmpz_init(d);

    /* A[i][j] = (i + 1)^j, plus 100 on the diagonal */
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
        {
            fmpz_set_ui(fmpz_mat_entry(A, i, j), n_pow(i + 1, j));
            if (i == j)
                fmpz_add_ui(fmpz_mat_entry(A, i, j), fmpz_mat_entry(A, i, j), 100);
        }

    /* |det A| <= bound (Hadamard's inequality); recover det from residues
       modulo primes until their product exceeds 2 * bound */
    fmpz_mat_det_bound(bound, A);
    fmpz_mul_2exp(bound, bound, 1);
    fmpz_one(modulus);
    fmpz_zero(det);
    p = UWORD(1) << (FLINT_BITS - 2);

    while (fmpz_cmp(modulus, bound) <= 0)
    {
        p = n_nextprime(p, 1);
        nmod_mat_init(Amod, N, N, p);
        fmpz_mat_get_nmod_mat(Amod, A);
        r = nmod_mat_det(Amod);
        fmpz_CRT_ui(det, det, modulus, r, p, 1);   /* 1: symmetric range */
        fmpz_mul_ui(modulus, modulus, p);
        nmod_mat_clear(Amod);
        nprimes++;
    }

    fmpz_mat_det(d, A);
    flint_printf("bound: %wu bits, primes used: %wd\n", fmpz_bits(bound), nprimes);
    flint_printf("det (multimodular) = %{fmpz}\n", det);
    flint_printf("det (fmpz_mat_det) = %{fmpz}\n", d);

    fmpz_clear(det);
    fmpz_clear(bound);
    fmpz_clear(modulus);
    fmpz_clear(d);
    fmpz_mat_clear(A);
    return 0;
}
```

```output
bound: 98 bits, primes used: 2
det (multimodular) = 230633716420959157248000
det (fmpz_mat_det) = 230633716420959157248000
```

The entries of an `fmpz_mat_t` are accessed with `fmpz_mat_entry(A, i,
j)`, which returns a pointer to the `fmpz` in row $i$ and column $j$.
`fmpz_mat_get_nmod_mat` reduces all entries modulo the modulus of the
`nmod_mat_t`, which is given when the latter is initialised. The second
argument of `n_nextprime` asks for a proven prime. The determinant has 78
bits, but the bound has 98, so two primes are needed.

## Refinements {#refinements}

The program uses enough primes for the worst case. FLINT's
`fmpz_mat_det_modular` does the same when its argument `proved` is set;
otherwise it also stops once the reconstructed value has not changed
while more than 100 bits of primes were added, which gives the correct
result with high probability. `fmpz_mat_det_modular_accelerated` first
computes a large divisor $d$ of the determinant by solving a linear
system, so that only $\det(A)/d$ has to be reconstructed, which needs
fewer primes. For many primes, reducing all entries modulo
all primes at once and reconstructing with a subproduct tree
(`fmpz_mat_multi_mod_ui`, `fmpz_multi_CRT_ui`) is faster than the one
prime at a time loop above, and these steps run in parallel when FLINT
is allowed to use several threads.
