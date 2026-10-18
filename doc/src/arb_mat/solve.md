---
title: Solving linear systems with ball matrices
functions: arb_mat_solve arb_mat_solve_preapprox arb_mat_inv arb_mat_det
---

`arb_mat_solve(X, A, B, prec)` solves $AX = B$ for a square matrix $A$
of real balls and returns nonzero if it has *proved* that every matrix
in $A$ is invertible; `X` then contains $A^{-1} B$ for every choice of
$A$ and $B$ in the input balls. A return value of zero means that this
could not be proved: $A$ may be singular, or the input or the working
precision may not be accurate enough. `arb_mat_inv` solves $AX = I$ in
the same way, and `arb_mat_det` encloses the determinant (it cannot
fail, but its enclosure may contain zero). Two families of algorithms
are available, and the difference between them is large.

::: summary algorithms/ball-arithmetic/linear-systems :::

## Gaussian elimination in ball arithmetic {#lu}

`arb_mat_solve_lu` computes an LU decomposition with partial pivoting
directly in ball arithmetic (`arb_mat_lu`) and solves the two
triangular systems. It succeeds if every pivot excludes zero. This is
the fastest method, and at high precision the best one, but in ball
arithmetic the radii grow with every elimination step: the computed
balls must contain the result for all matrices in the input balls *and*
the dependencies between entries are lost at each step, so the radii
typically grow exponentially with $n$ even for well-conditioned
matrices. At 64-bit precision, the enclosures of a random $40 \times 40$
system are already useless and larger systems fail entirely (see the
example below). The method is therefore used by default only when
$\mathit{prec} > 10n$ or $n \le 4$.

## Preconditioning {#precond}

`arb_mat_solve_precond` first computes an approximate inverse $R
\approx A^{-1}$ of the midpoint matrix in floating-point arithmetic
(`arb_mat_approx_solve` with $B = I$). Then $RA$ is close to the
identity matrix, and the system $(RA) X = RB$ has the same solution as
$AX = B$ if $R$ is invertible. Following Hansen and Smith [@HS1967],
the products $RA$ and $RB$ are computed in ball arithmetic and the new
system is solved:

- If $B$ has few columns ($m < 0.1n + 1$), $(RA)X = RB$ is solved by
  ball Gaussian elimination. Since $RA$ is close to $I$, the pivots
  stay close to 1 and the radii grow only mildly.
- Otherwise, with $E = I - RA$, if $\|E\|_\infty \le \delta < 1$ then
  $RA$ is invertible and $X = (I - E)^{-1} RB$ satisfies
  $\|X_j - (RB)_j\|_\infty \le \frac{\delta}{1-\delta} \|(RB)_j\|_\infty$
  for each column $j$, by the Neumann series. The result is $RB$ with
  this bound added to every entry. If $\delta \ge 1$, the first method
  is used instead.

The cost is about that of three matrix multiplications and an
approximate inversion, several times that of `arb_mat_solve_lu`, but
the radii of the result are of the order of $\kappa(A)$ times the
radii of the input and the unit roundoff, as in a numerical backward
error analysis, and do not grow exponentially with $n$. This is the
default for $n > 4$ unless $\mathit{prec} > 10n$. The theory is in
[](algorithms/ball-arithmetic/linear-systems).

## Certifying a given approximation {#preapprox}

`arb_mat_solve_preapprox(X, A, B, R, T, prec)` certifies an approximate
solution $T$ computed by any means, given an approximate inverse $R$.
By Theorem 10.2 of [@Rum2010], if $\|I - RA\|_\infty < 1$ then $A$ is
invertible and each column of the exact solution satisfies
$$\|X_j - T_j\|_\infty \le \frac{\|R(AT_j - B_j)\|_\infty}{1 - \|I - RA\|_\infty}.$$
The function computes $RA$ and the residual $AT - B$ in ball
arithmetic, so the result is rigorous even though $R$ and $T$ are not.
It returns zero if $\|I - RA\|_\infty < 1$ cannot be shown, for example
because $R$ is a poor approximation. This is the cheapest way to obtain
error bounds when good approximations are already available, for
example from a previous solve at lower precision.

## Inverse and determinant {#det}

`arb_mat_inv` sets the output to the identity and calls `arb_mat_solve`,
so it inherits the choice of algorithm.

`arb_mat_det` uses explicit formulas for $n \le 3$ and the product of
the diagonal for triangular matrices. Otherwise, for $n \le 10$ or
$\mathit{prec} > 10n$, it calls `arb_mat_det_lu`, which performs ball
Gaussian elimination as far as pivots can be found, multiplies the
pivots, and bounds the determinant of the remaining submatrix by
Hadamard's inequality $|\det M| \le \prod_i \|m_i\|_2$ over its rows. For
larger matrices it calls `arb_mat_det_precond`, following Rump
[@Rum2010]: from an approximate decomposition $PA \approx LU$ it forms
$M = L^{-1} P A U^{-1}$ in ball arithmetic with approximate inverses of
the triangular factors, which is close to the identity, encloses
$\det M \in [1 \pm (e^{n\varepsilon} - 1)]$ with $\varepsilon$ the
largest Gershgorin radius $\max_i \sum_j |m_{ij} - \delta_{ij}|$ (since
every eigenvalue of $M$ lies within $\varepsilon$ of 1), and
divides by $\det U^{-1}$ (the product of the diagonal of the approximate
$U^{-1}$; $L^{-1}$ has unit diagonal). If the result contains zero,
the LU bound is also computed and the better enclosure is returned.

## Examples {#examples}

The program solves $AX = B$ for random matrices with entries in
$[-1, 1]$ and prints the largest radius of the solution vector. The
radii from ball Gaussian elimination grow by a factor of about $10^4$
each time $n$ doubles, until no pivot can be certified; those of the
preconditioned solver stay near the unit roundoff:

```c example name=arb_mat_solve_compare
#include <flint/ulong_extras.h>
#include <flint/arb.h>
#include <flint/arb_mat.h>

/* largest radius of the entries of X */
static double max_rad(const arb_mat_t X)
{
    slong i, j;
    double r = 0;
    for (i = 0; i < arb_mat_nrows(X); i++)
        for (j = 0; j < arb_mat_ncols(X); j++)
            r = FLINT_MAX(r, mag_get_d(arb_radref(arb_mat_entry(X, i, j))));
    return r;
}

int main(void)
{
    arb_mat_t A, B, X;
    flint_rand_t state;
    slong n, i, j, prec = 64;

    flint_rand_init(state);

    for (n = 10; n <= 80; n *= 2)
    {
        arb_mat_init(A, n, n);
        arb_mat_init(B, n, 1);
        arb_mat_init(X, n, 1);

        /* random entries in [-1, 1] */
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                arb_set_si(arb_mat_entry(A, i, j), (slong) n_randint(state, 2001) - 1000);
                arb_div_ui(arb_mat_entry(A, i, j), arb_mat_entry(A, i, j), 1000, prec);
            }
        arb_mat_ones(B);

        flint_printf("n = %wd:", n);
        if (arb_mat_solve_lu(X, A, B, prec))
            flint_printf("  lu %.2e", max_rad(X));
        if (arb_mat_solve_precond(X, A, B, prec))
            flint_printf("  precond %.2e", max_rad(X));
        flint_printf("\n");

        arb_mat_clear(A);
        arb_mat_clear(B);
        arb_mat_clear(X);
    }
    flint_rand_clear(state);
    return 0;
}
```
```output
n = 10:  lu 3.74e-15  precond 2.10e-18
n = 20:  lu 1.39e-09  precond 2.08e-16
n = 40:  lu 3.34e-01  precond 7.47e-17
n = 80:  precond 2.87e-16
```

The second program computes an approximate solution of a Hilbert
system with the floating-point functions and certifies it with
`arb_mat_solve_preapprox`. The exact solution of $H_6 x = (1, \ldots,
1)^T$ is the vector of column sums of $H_6^{-1}$, which are integers:

```c example name=arb_mat_solve_preapprox
#include <flint/arb.h>
#include <flint/arb_mat.h>

int main(void)
{
    arb_mat_t A, B, R, T, X;
    slong n = 6, prec = 128;

    arb_mat_init(A, n, n);
    arb_mat_init(B, n, 1);
    arb_mat_init(R, n, n);
    arb_mat_init(T, n, 1);
    arb_mat_init(X, n, 1);

    arb_mat_hilbert(A, prec);
    arb_mat_ones(B);

    /* floating-point approximations, no error bounds */
    arb_mat_approx_inv(R, A, prec);
    arb_mat_approx_mul(T, R, B, prec);

    /* certify: X = T + error bound */
    if (arb_mat_solve_preapprox(X, A, B, R, T, prec))
        arb_mat_printd(X, 15);
    else
        flint_printf("failed\n");

    arb_mat_clear(A);
    arb_mat_clear(B);
    arb_mat_clear(R);
    arb_mat_clear(T);
    arb_mat_clear(X);
    return 0;
}
```
```output
[-6.00000000000000 +/- 8.0698e-29]
[210.000000000000 +/- 8.0698e-29]
[-1680.00000000000 +/- 8.0698e-29]
[5040.00000000000 +/- 8.0698e-29]
[-6300.00000000000 +/- 8.0698e-29]
[2772.00000000000 +/- 8.0698e-29]
```

The Hilbert matrix itself is only known to 128 bits (its entries $1/(j
+ k + 1)$ are rounded), and the certified solution is the solution for
every matrix in those balls; the radius $8 \cdot 10^{-29}$ reflects the
condition number $\kappa(H_6) \approx 1.5 \cdot 10^7$ times the rounding
of the entries.
