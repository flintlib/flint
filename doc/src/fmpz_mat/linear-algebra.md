---
title: Exact linear algebra over the integers
summary: How integer matrices and rational results are represented, which algorithms the main functions use and what they cost, which results are proven, and how to choose between rational, modular and lattice methods.
order: 1
---

This guide explains how to use `fmpz_mat` for exact linear algebra with
integer matrices, and how to combine it with `nmod_mat`, `fmpq_mat` and
`fmpz_lll`.

## Storage {#storage}

An `fmpz_mat_t` holds a pointer to its first entry, its dimensions and a
row stride. Entry $(i, j)$ is the `fmpz` at offset
$i \cdot \mathit{stride} + j$, so each row is a contiguous `fmpz` vector
that can be passed to the functions of [](fmpz_vec), for example
`_fmpz_vec_content(g, fmpz_mat_row(A, i), A->c)`. Entries are read and
written through `fmpz_mat_entry`, which returns a pointer usable with every
function of [](fmpz).

A *window* (`fmpz_mat_window_init`) refers to a block of another matrix
without copying it and can be used as input or output of most functions.
Windows must be released with `fmpz_mat_window_clear` and never passed to
`fmpz_mat_clear`, `fmpz_mat_swap` or other functions that free or replace
the storage of their argument.

Small entries (up to 62 bits on a 64-bit machine) are stored inline in
the `fmpz` word, and many kernels of this module read them directly. A
matrix with small entries is therefore several times faster to work with
than one with entries just above that size.

## Rational results {#denominators}

Solutions of linear systems, inverses and reduced echelon forms have
rational entries. They are returned as an integer matrix $X$ and an
integer `den`, the result being $X/\mathit{den}$. The denominator depends
on the algorithm:

| Function | `den` |
|---|---|
| `fmpz_mat_solve_cramer`, `fmpz_mat_solve_fflu`, `fmpz_mat_inv` for $n \le 15$ | $\det(A)$ |
| `fmpz_mat_solve_dixon_den`, `fmpz_mat_solve_multi_mod_den` | least common denominator, positive |
| `fmpz_mat_solve`, `fmpz_mat_inv` | one of the two, depending on the size |
| `fmpz_mat_can_solve` and variants | as the solver used for the pivot block |
| `fmpz_mat_rref_fflu`, `fmpz_mat_fflu` | a nonzero $r \times r$ minor of $A$ |
| `fmpz_mat_rref_mul` | 1 for full column rank, otherwise as `fmpz_mat_solve` |

To compare results or to get canonical output, divide `den` and the
entries by their gcd (`fmpz_mat_content`), or convert to an
`fmpq_mat_t` with `fmpq_mat_set_fmpz_mat_div_fmpz`, which reduces every
entry. `den` may be negative for the fraction-free functions.

For the inverse, the least common denominator returned by the modular
solvers is the largest invariant factor $s_n$ of $A$ (the last entry of the Smith form),
which is often much smaller than $\det(A)$.

## Sizes {#sizes}

The cost of exact linear algebra is governed by the sizes of the
results, which follow from Hadamard's inequality: for an $n \times n$
matrix with entries bounded by $2^b$,

$$
|\det(A)| \le \prod_i \lVert a_i \rVert_2 \le n^{n/2} 2^{nb},
$$

so determinants, the entries of $\det(A) A^{-1}$, the coefficients of the
characteristic polynomial and the entries of fraction-free echelon forms
have up to about $n(b + \frac12 \log_2 n)$ bits. `fmpz_mat_det_bound`,
`fmpz_mat_solve_bound` and `fmpz_mat_charpoly_bound` compute such bounds.
The input has $n^2 b$ bits and a solution vector $n \cdot n(b + \log n)$
bits, so solving one system is a computation with output as large as the
input; the inverse has $n^3 (b + \log n)$ bits, $n$ times the input.

Actual results are often much smaller than the bounds, for example when
the system has a small solution by construction. The modular solvers are
output-sensitive: their running time depends on the size of the actual
solution rather than on the bound. Determinant and characteristic
polynomial computations are not output-sensitive, except for the early
detection of zero determinants and the optional unproven early
termination of `fmpz_mat_det_modular`.

## Choosing functions {#choosing}

The functions without algorithm suffix choose a method from the
dimensions and the entry sizes, and should be the default:

- `fmpz_mat_mul` covers sizes from a few words to millions of bits; see
  [](fmpz_mat/mul).
- `fmpz_mat_solve` for nonsingular square systems, `fmpz_mat_can_solve`
  for arbitrary ones; see [](fmpz_mat/solve).
- `fmpz_mat_det`, and `fmpz_mat_is_singular` when only singularity
  matters: the latter needs one or two primes for nonsingular input,
  against $O(n(b + \log n))/59$ primes for the determinant.
- `fmpz_mat_rank`, `fmpz_mat_rref` and `fmpz_mat_nullspace`.
- `fmpz_mat_charpoly` and `fmpz_mat_minpoly`.
- `fmpz_mat_hnf`, `fmpz_mat_snf` and `fmpz_mat_elementary_divisors` for
  lattices and abelian groups.

Some patterns to avoid:

- Do not compute an inverse to solve a system: `fmpz_mat_inv` costs as
  much as solving with $n$ right-hand sides, and its output is $n$ times
  larger than a solution.
- Do not use the fraction-free functions (`fmpz_mat_det_bareiss`,
  `fmpz_mat_solve_fflu`, `fmpz_mat_rref_fflu`) for large matrices: their
  entries grow to the size of minors at every step, which makes them
  $O(n^5 b^2)$ algorithms against $O(n^4 b)$ or better for the modular
  ones.
- Solve with many right-hand sides at once: a $p$-adic or multimodular
  solve with $k$ columns costs much less than $k$ solves.
- Work with `fmpq_mat_t` only at the interface. Rational matrices store a
  numerator and denominator per entry and canonicalise after every
  operation; the `fmpq_mat` solvers clear denominators and call the
  functions of this module.

## Proofs and randomness {#proofs}

All functions of this module return proven results, with two exceptions:
`fmpz_mat_det_modular` and its variants with `proved` = 0 (never called
that way by `fmpz_mat_det`), and the reducedness tests of [](fmpz_lll)
whose names end in `_d` or `_mpfr`, which are one-sided. Several
algorithms use randomness, for example random primes, random projections
of solutions and random linear combinations, but only to decide what to
compute; the results are checked exactly or follow from proven bounds.

Behind this are three kinds of certificates:

- *Size bounds.* A matrix $X$ known modulo $m$ is determined if $m$
  exceeds twice a bound on its entries; a rational matrix is determined if
  $m > 2ND$ for bounds $N$ and $D$ on numerators and denominator. Bounds
  come from Hadamard's inequality.
- *Residual checks.* A candidate solution $X/d$ with
  $AX \equiv dB \pmod m$ is exact if
  $n \lVert A \rVert \lVert X \rVert + d \lVert B \rVert < m$, or if
  $AX = dB$ is verified by a multiplication.
- *Witnesses.* A nonzero kernel vector proves singularity, a vector $y$
  with $y^T A = 0$ and $y^T B \ne 0$ proves inconsistency, and a reduced
  echelon form whose nullspace annihilates all rows proves the rank.

[](algorithms/linear-algebra/multimodular) gives the details.

## Modular computations {#modular}

Many computations over $\Z$ are best done modulo several primes by the
user, with the tools that the module provides:

- `fmpz_mat_get_nmod_mat` and `fmpz_mat_multi_mod_ui` reduce a matrix
  modulo one or many moduli (the latter with a subproduct tree, and in
  parallel);
- [](nmod_mat) does linear algebra modulo each prime at the speed of
  floating-point arithmetic ([](nmod_mat/linear-algebra)): products
  ([](nmod_mat/mul)), LU decompositions ([](nmod_mat/lu)), solving
  ([](nmod_mat/solve)), determinants ([](nmod_mat/det)) and
  characteristic polynomials ([](nmod_mat/charpoly));
- `fmpz_mat_CRT_ui` combines one more residue into a partial result, and
  `fmpz_mat_multi_CRT_ui` combines many residues at once, in the
  symmetric range when the result may be negative.

Rational reconstruction of the result, with `fmpq_reconstruct_fmpz` per
entry or `fmpq_mat_set_fmpz_mat_mod_fmpz` for a whole matrix, then gives
rational results. A modulus is unlucky when it divides a determinant that
the computation relies on; the functions of this module handle that by
certification, and user code should either do the same or use enough
primes.

## Lattices {#lattices}

The rows of an integer matrix also define a lattice, an additive subgroup
of $\Z^n$, and several functions work with lattices rather than vector
spaces:

- `fmpz_mat_hnf` gives the canonical basis (Hermite normal form), from
  which equality and inclusion of lattices can be decided
  ([](fmpz_mat/hnf));
- `fmpz_mat_snf` and `fmpz_mat_elementary_divisors` give the structure of
  the quotient group $\Z^n / L$ ([](fmpz_mat/snf));
- `fmpz_lll` computes a reduced basis consisting of short, nearly
  orthogonal vectors ([](fmpz_lll), [](algorithms/linear-algebra/lll)).

`fmpz_mat_nullspace` returns a basis of the nullspace over $\Q$, scaled
to integers. This is in general not a basis of the lattice of integer
solutions. That lattice is the left kernel of $A^T$ and can be computed
from a Hermite normal form with transformation: if $UA^T = H$ and the
last $m - r$ rows of $H$ are zero, then the last $m - r$ rows of $U$ are
a basis of $\{x \in \Z^m : x^T A^T = 0\}$. LLL reduction of that basis
then gives short kernel vectors. The example shows both.

## Example {#example}

The integer solutions of $Ax = 0$ for a $2 \times 4$ matrix $A$, first as
the rational nullspace basis, then as a basis of the integer kernel
lattice obtained from the Hermite normal form of $A^T$ with
transformation, and finally LLL-reduced:

```c example name=fmpz_mat_integer_kernel
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>

int main(void)
{
    fmpz_mat_t A, At, N, H, U, K;
    fmpz_lll_t fl;
    slong i, j, nullity;
    slong a[2][4] = { { 6, 10, 15, 7 }, { 4, 0, 6, -2 } };

    fmpz_mat_init(A, 2, 4);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 4; j++)
            fmpz_set_si(fmpz_mat_entry(A, i, j), a[i][j]);

    /* rational nullspace: basis vectors are the columns of N */
    fmpz_mat_init(N, 4, 4);
    nullity = fmpz_mat_nullspace(N, A);
    flint_printf("nullity %wd, nullspace basis (columns):\n", nullity);
    fmpz_mat_print_pretty(N);
    flint_printf("\n");

    /* integer kernel lattice from the HNF of A^T with transformation */
    fmpz_mat_init(At, 4, 2);
    fmpz_mat_init(H, 4, 2);
    fmpz_mat_init(U, 4, 4);
    fmpz_mat_transpose(At, A);
    fmpz_mat_hnf_transform(H, U, At);

    fmpz_mat_init(K, nullity, 4);
    for (i = 0; i < nullity; i++)
        for (j = 0; j < 4; j++)
            fmpz_set(fmpz_mat_entry(K, i, j),
                     fmpz_mat_entry(U, 4 - nullity + i, j));

    fmpz_lll_context_init_default(fl);
    fmpz_lll(K, NULL, fl);
    flint_printf("reduced basis of the integer kernel (rows):\n");
    fmpz_mat_print_pretty(K);
    flint_printf("\n");

    fmpz_mat_clear(A);
    fmpz_mat_clear(At);
    fmpz_mat_clear(N);
    fmpz_mat_clear(H);
    fmpz_mat_clear(U);
    fmpz_mat_clear(K);
    return 0;
}
```
```output
nullity 2, nullspace basis (columns):
[[-60 20 0 0]
[-24 -40 0 0]
[40 0 0 0]
[0 40 0 0]
]
reduced basis of the integer kernel (rows):
[[-1 2 0 -2]
[-7 -4 5 1]
]
```

The nullspace basis consists of the columns $(-60, -24, 40, 0)$ and
$(20, -40, 0, 40)$: each is a multiple of a kernel vector, scaled by the
denominator 40 of the echelon form. The two rows of the second result
generate all integer solutions.
