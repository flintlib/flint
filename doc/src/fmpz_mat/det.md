---
title: Computing determinants
functions: fmpz_mat_det fmpz_mat_det_modular
---

`fmpz_mat_det` computes the determinant of a square integer matrix
exactly. For small dimensions it uses cofactor expansion or fraction-free
elimination, whose intermediate entries never exceed the size of minors
of the matrix. For larger dimensions it computes the determinant modulo
word-size primes, where elimination costs $O(n^\omega)$ word operations
per prime, and reconstructs it by the Chinese remainder theorem. The
number of primes is governed by Hadamard's bound, unless a large divisor
of the determinant is found first, in which case only the cofactor needs
to be reconstructed. All results are proven unless `proved` is set to 0
in the modular functions.

::: summary algorithms/linear-algebra/multimodular :::

## Choice of algorithm {#dispatch}

For an $n \times n$ matrix whose largest entry has $b$ bits,
`fmpz_mat_det` chooses as follows.

| Dimension | Method |
|---|---|
| $n \le 4$ | `fmpz_mat_det_cofactor`: explicit cofactor expansion |
| $5 \le n \le 24$ | `fmpz_mat_det_bareiss`: fraction-free elimination |
| $25 \le n \le 59$ | `fmpz_mat_det_modular` with `proved` = 1 |
| $n \ge 60$, $n < b$ | `fmpz_mat_det_modular` with `proved` = 1 |
| $n \ge 60$, $n \ge b$ | `fmpz_mat_det_modular_accelerated` with `proved` = 1 |

The empty matrix has determinant 1. A nonsquare matrix raises an
exception in `fmpz_mat_det`; the other functions assume squareness.

## Fraction-free elimination {#bareiss}

`fmpz_mat_det_bareiss` runs `fmpz_mat_fflu` with rank check on a copy of
$A$. By Sylvester's identity the $k$-th pivot of Bareiss' algorithm is the
leading $k \times k$ minor of the row-permuted matrix, so the last pivot
is $\pm \det(A)$, with the sign of the permutation, and all
intermediate entries are minors of $A$, bounded by Hadamard's inequality
[@Bar1968]. The cost is $O(n^3)$ multiplications and exact divisions of
integers of up to $n(b + \log_2 n)$ bits, about $O(n^5 b^2)$ word
operations with classical integer multiplication, which is why it is used
only for small $n$; see [](algorithms/linear-algebra/fraction-free).

## Modular computation {#modular}

`fmpz_mat_det_modular_given_divisor(det, A, d, proved)` computes
$x = \det(A)/d$ for a known divisor $d$ of the determinant (1 in
`fmpz_mat_det_modular`):

1. $H$ = `fmpz_mat_det_bound(A)` (Hadamard's bound), and the target is
   $2H/d$, rounded up.
2. For primes $p$ above $2^{59}$ that do not divide $d$: reduce $A$
   modulo $p$, compute its LU decomposition with
   `nmod_mat_lu_with_pivots`, and set $x_p = \pm \prod_i U_{ii} \cdot
   d^{-1} \bmod p$.
3. Combine the $x_p$ with `fmpz_CRT_ui` in the symmetric range until the
   product of the primes exceeds $2H/d$. Then $|x| \le H/d$ is determined.

If the rank modulo some $p$ is below $n$, the function tries to prove
$\det(A) = 0$ at once with an exact kernel vector
(`_fmpz_mat_certify_singular_lu_mod_p`), which succeeds whenever the rank
modulo $p$ equals the rank of $A$. With `proved` = 0 the loop also stops
when the reconstructed value has stayed the same while more than 100
bits of primes were added: the probability that a wrong value survives
that many primes is negligible for nonadversarial input, but the result
is not proven.

The cost is about $\log_2 (2H/d) / 59$ LU decompositions modulo primes,
each $O(n^\omega)$, plus the reduction of $A$, $O(n^2)$ reductions of
$b$-bit entries per prime. With $\log_2 H \approx n (b + \frac12 \log_2
n)$, this is $O(n^{\omega + 1} (b + \log n))$ word operations for
`fmpz_mat_det_modular`.

## Large divisors {#divisor}

`fmpz_mat_det_modular_accelerated` first calls `fmpz_mat_det_divisor`,
which solves $Ax = v$ for a vector $v$ with entries $\pm 1$ by $p$-adic
lifting and returns the least common multiple $d$ of the denominators of
$x$. This $d$ divides the largest invariant factor $s_n$ of $A$, and for
random matrices it equals $s_n$, which is $\det(A)$ up to a small factor
[@AbbottBronsteinMulders1999; @Pan1988]. Only $\det(A)/d$ remains to be
reconstructed, which typically takes one or two primes. Since the
$p$-adic solution costs about the same as $O(n)$ primes in the purely
modular method but saves $O(n(b + \log n))/59$ primes, this pays off when
the determinant is large compared with $n$, which `fmpz_mat_det`
estimates by $n \ge b$. For matrices with a large Smith form
cofactor $\det(A)/s_n$, for example diagonal-dominant structured
matrices, the gain disappears; see
[](algorithms/linear-algebra/multimodular#det).

## Related functions {#related}

- `fmpz_mat_det_bound` computes Hadamard's bound, and
  `fmpz_mat_det_bound_submatrix` a bound valid for all minors.
- `fmpz_mat_is_singular` decides $\det(A) = 0$ without computing the
  determinant, with a prime of 17 bits and an exact certificate; it is
  much faster than `fmpz_mat_det` for large nonsingular matrices.
- `fmpz_mat_det_divisor` returns the divisor described above.
- `nmod_mat_det` computes determinants modulo $n$, and
  `fmpq_mat_det` determinants of rational matrices by clearing
  denominators row by row.

## Example {#example}

A dense $40 \times 40$ matrix with prescribed determinant
$2^{40} \cdot 3^{20} \cdot 1000003$, obtained from a sparse matrix by
random unimodular row and column operations. All methods agree, and the
Hadamard bound is far larger than the determinant. Since `fmpz_mat_randdet`
spreads the prime factors of the determinant over many diagonal entries,
the largest invariant factor, which is what `fmpz_mat_det_divisor` finds,
is only $36 \cdot 1000003$ here; for such matrices the accelerated method
gains nothing.

```c example name=fmpz_mat_det_methods
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

int main(void)
{
    fmpz_mat_t A;
    fmpz_t D, d, t;
    flint_rand_t state;
    slong n = 40;

    flint_rand_init(state);
    fmpz_init(D);
    fmpz_init(d);
    fmpz_init(t);
    fmpz_mat_init(A, n, n);

    /* D = 2^40 * 3^20 * 1000003 */
    fmpz_ui_pow_ui(D, 3, 20);
    fmpz_mul_2exp(D, D, 40);
    fmpz_mul_ui(D, D, 1000003);

    fmpz_mat_randdet(A, state, D);
    fmpz_mat_randops(A, state, 2000);
    flint_printf("entries up to %wd bits\n", FLINT_ABS(fmpz_mat_max_bits(A)));

    fmpz_mat_det(d, A);
    flint_printf("det:                %{fmpz}\n", d);
    flint_printf("correct:            %d\n", fmpz_equal(d, D));

    fmpz_mat_det_bareiss(t, A);
    flint_printf("bareiss agrees:     %d\n", fmpz_equal(t, D));
    fmpz_mat_det_modular(t, A, 1);
    flint_printf("modular agrees:     %d\n", fmpz_equal(t, D));
    fmpz_mat_det_modular_accelerated(t, A, 1);
    flint_printf("accelerated agrees: %d\n", fmpz_equal(t, D));

    fmpz_mat_det_divisor(t, A);
    flint_printf("divisor:            %{fmpz}\n", t);
    flint_printf("divides det:        %d\n", fmpz_divisible(D, t));

    fmpz_mat_det_bound(t, A);
    flint_printf("bits of det %wd, of Hadamard bound %wd\n",
                 fmpz_bits(D), fmpz_bits(t));

    fmpz_mat_clear(A);
    fmpz_clear(D);
    fmpz_clear(d);
    fmpz_clear(t);
    flint_rand_clear(state);
    return 0;
}
```
```output
entries up to 47 bits
det:                3833771493727452464601366528
correct:            1
bareiss agrees:     1
modular agrees:     1
accelerated agrees: 1
divisor:            36000108
divides det:        1
bits of det 92, of Hadamard bound 1865
```
