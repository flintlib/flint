---
title: Getting started with arb_mat
summary: How to create, fill and print matrices of real balls, what the return values of the verified functions mean, and when to use the approximate floating-point functions.
order: 1
---

An `arb_mat_t` is a dense matrix of real balls with dimensions fixed at
initialisation. Each entry is an ordinary `arb_t` accessed by the macro
`arb_mat_entry(A, i, j)`, so all functions of [](arb) can be used on
entries; programs that work with the entries must also include
`arb.h`. This guide explains the conventions that differ from exact
matrix modules such as [](fmpz_mat).

## A first program {#first}

The program builds the symmetric matrix
$$A = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 3 & 1 \\ 0 & 1 & \pi \end{pmatrix},$$
computes its determinant $5\pi - 2$, proves that it is invertible and
positive definite, and computes its exponential:

```c example name=arb_mat_first
#include <flint/arb.h>
#include <flint/arb_mat.h>

int main(void)
{
    arb_mat_t A, B, C, I;
    arb_t d;
    slong prec = 64;

    arb_mat_init(A, 3, 3);
    arb_mat_init(B, 3, 3);
    arb_mat_init(C, 3, 3);
    arb_mat_init(I, 3, 3);
    arb_init(d);

    /* A = [[2, 1, 0], [1, 3, 1], [0, 1, pi]] */
    arb_set_si(arb_mat_entry(A, 0, 0), 2);
    arb_set_si(arb_mat_entry(A, 0, 1), 1);
    arb_set_si(arb_mat_entry(A, 1, 0), 1);
    arb_set_si(arb_mat_entry(A, 1, 1), 3);
    arb_set_si(arb_mat_entry(A, 1, 2), 1);
    arb_set_si(arb_mat_entry(A, 2, 1), 1);
    arb_const_pi(arb_mat_entry(A, 2, 2), prec);

    arb_mat_det(d, A, prec);
    flint_printf("det(A) = ");
    arb_printn(d, 15, 0);
    flint_printf("\n");

    if (arb_mat_inv(B, A, prec))
    {
        arb_mat_mul(C, A, B, prec);
        arb_mat_one(I);
        flint_printf("A * A^-1 contains I: %d\n", arb_mat_contains(C, I));
    }

    if (arb_mat_cho(B, A, prec))
        flint_printf("A is positive definite\n");

    arb_mat_exp(B, A, prec);
    flint_printf("exp(A) =\n");
    arb_mat_printd(B, 10);

    arb_mat_clear(A);
    arb_mat_clear(B);
    arb_mat_clear(C);
    arb_mat_clear(I);
    arb_clear(d);
    return 0;
}
```
```output
det(A) = [13.7079632679490 +/- 3.39e-14]
A * A^-1 contains I: 1
A is positive definite
exp(A) =
[13.83369352 +/- 8.8276e-19, 18.06134102 +/- 1.8257e-18, 9.197078786 +/- 1.2278e-18]
[18.06134102 +/- 1.8257e-18, 41.09211333 +/- 3.9362e-18, 28.56065860 +/- 3.2343e-18]
[9.197078786 +/- 1.2278e-18, 28.56065860 +/- 3.2343e-18, 35.93901398 +/- 6.8848e-18]
```

`arb_mat_init` sets all entries to zero, so only the nonzero entries
need to be set. The entry $\pi$ is a ball of radius about $2^{-64}$, and
every result encloses the corresponding quantity for all matrices
within that ball. The product $A A^{-1}$ is not exactly the identity,
but it contains it.

## Return values and failure {#failure}

Ball arithmetic can prove that a matrix is invertible or positive
definite, but never that it is singular: a ball matrix that contains a
singular matrix also contains invertible ones. The functions that need
such a property therefore return an `int`:

- `arb_mat_solve`, `arb_mat_inv`, `arb_mat_lu` and their variants
  return nonzero if the matrix is proved invertible, and zero
  otherwise; `arb_mat_cho`, `arb_mat_ldl`, `arb_mat_spd_solve` and
  `arb_mat_spd_inv` return nonzero if it is proved positive definite.
  After a zero return value the output is undefined (it is not set to
  indeterminate balls). The usual reaction is to repeat the computation
  at a higher precision, with more accurate input if necessary.
- `arb_mat_det` always succeeds; for a singular or nearly singular
  matrix its result contains zero.
- Predicates such as `arb_mat_eq` and `arb_mat_is_zero` return 1 only
  if the property certainly holds.

For linear systems, the choice of algorithm matters a great deal: plain
Gaussian elimination in ball arithmetic loses accuracy exponentially in
the dimension, and `arb_mat_solve` uses a preconditioned method for all
but small matrices or very high precision. See [](arb_mat/solve).

## Approximate functions {#approx}

The functions `arb_mat_approx_mul`, `arb_mat_approx_lu`,
`arb_mat_approx_solve`, `arb_mat_approx_inv` and the triangular
`arb_mat_approx_solve_triu` and `arb_mat_approx_solve_tril` ignore the
radii and compute floating-point approximations with no error bound.
They are faster than their rigorous counterparts and are the building
blocks of the verified algorithms: an approximate inverse is the
preconditioner of `arb_mat_solve_precond` and `arb_mat_det_precond`,
and `arb_mat_solve_preapprox` certifies an approximate solution
computed by any method. They can also be used for ordinary numerical
linear algebra at arbitrary precision when error bounds are not needed.

## Dimensions, aliasing and windows {#dimensions}

The caller must pass matrices of compatible dimensions. The
multiplication, solving and decomposition functions throw an exception
otherwise, but entrywise functions such as `arb_mat_add` do not check
the dimensions. Outputs may be the same objects as inputs
in the arithmetic and solving functions, which use temporary storage
when needed. `arb_mat_window_init` creates a view of a submatrix that
shares the entries of the original matrix, so that operations on
blocks need no copying; a window must be released with
`arb_mat_window_clear` and must not outlive its parent.

The complex counterpart [](acb_mat) has the same interface and adds
eigenvalues and eigenvectors ([](acb_mat/eigenvalues)). Matrices over
the generic real ball context of [](gr_mat) share the memory layout of
`arb_mat_t`, so that the generic functions can be applied to an
`arb_mat_t` by a cast, as this module does internally for the
characteristic polynomial and for powers.
