---
title: Solving linear systems over the rationals
functions: fmpz_mat_solve fmpz_mat_can_solve
---

The solvers of this module find solutions over $\Q$ of $AX = B$ for
integer matrices $A$ and $B$. A solution is returned as an integer matrix
$X$ and a denominator `den` with $AX = \mathit{den} \cdot B$, so that
$X/\mathit{den}$ is the rational solution. `fmpz_mat_solve` and its
variants require a square $A$ and decide whether it is nonsingular;
`fmpz_mat_can_solve` and its variants accept any $A$ and decide whether
the system is consistent. Every answer is proven, including "singular"
and "inconsistent".

::: summary algorithms/linear-algebra/multimodular :::

## Nonsingular systems {#solve}

For a square $A$ of dimension $n$ and an $n \times k$ matrix $B$,
`fmpz_mat_solve(X, den, A, B)` returns 1 and sets $X$ and `den` if $A$ is
nonsingular, and returns 0 if $A$ is singular. It chooses the algorithm
as follows.

| Condition | Function | Denominator returned |
|---|---|---|
| $n \le 3$ | `fmpz_mat_solve_cramer` | $\det(A)$ |
| $4 \le n \le 15$ | `fmpz_mat_solve_fflu` | $\det(A)$ |
| $n \ge 16$, $2k \ge n$ | `fmpz_mat_solve_multi_mod_den` | least common denominator, positive |
| $n \ge 16$, $2k < n$ | `fmpz_mat_solve_dixon_den` | least common denominator, positive |

The fraction-free solvers produce $\det(A) A^{-1} B$, whose entries are
determinants by Cramer's rule; dividing out the gcd of the entries and
`den` would give the reduced form. The modular solvers reconstruct the
reduced solution directly. $X$ may be $A$ or $B$ in all of them.

### Fraction-free elimination {#fflu}

`fmpz_mat_solve_cramer` evaluates Cramer's rule with explicit cofactor
formulas. `fmpz_mat_solve_fflu` computes a fraction-free LU decomposition
of $A$ with `fmpz_mat_fflu` and solves with fraction-free forward and back
substitution (`fmpz_mat_solve_fflu_precomp`), all of whose divisions are
exact; see [](algorithms/linear-algebra/fraction-free#solve). Entries grow
to the size of $n \times n$ minors, about $n(b + \log_2 n)$ bits for
$b$-bit entries, so the cost grows like $n^5 b^2$ and the method is
competitive only for small $n$.

### p-adic lifting {#dixon}

`fmpz_mat_solve_dixon_den` implements Dixon's algorithm [@Dix1982]. It
finds a prime $p$ near $2^{30}$ modulo which $A$ is nonsingular (or
proves $A$ singular), computes the LU decomposition of $A$ modulo $p$
once, and then computes the $p$-adic expansion $X = X_0 + X_1 p + X_2 p^2
+ \cdots$ of $A^{-1}B$, one digit per step: $X_i = A^{-1} R_i \bmod p$
and $R_{i+1} = (R_i - A X_i)/p$, exactly. Each step costs a triangular
solve modulo $p$ and a product of $A$ with an $n \times k$ matrix of
word-size entries. When the solution is long, the digits grow to $p^e$
for increasing $e$, with an inverse of $A$ modulo $p^e$ computed by Newton
iteration and products computed multimodularly from precomputed residues
of $A$.

The rational solution is recovered by rational reconstruction from $X$
modulo $p^i$. Lifting to $p^i > 2ND$, with the bounds $N$ and $D$ of
`fmpz_mat_solve_bound`, always suffices, but the solver stops as soon as
the true solution is determined, which for small solutions happens much
earlier; see [](#termination).

### Chinese remaindering {#multi-mod}

`fmpz_mat_solve_multi_mod_den` solves $AX = B$ modulo primes above
$2^{59}$ with `nmod_mat_solve`, processing batches of primes in parallel,
and combines the solutions modulo the product of the primes by the
Chinese remainder theorem, with the same reconstruction and termination
as the $p$-adic solver. A prime costs an LU decomposition and a triangular
solve, $O(n^\omega + n^2 k)$ word operations, against a product with $A$,
$O(n^2 k (b/\beta + 1))$, per $p$-adic step. Many right-hand sides
make the LU decomposition comparatively cheap, which is why
`fmpz_mat_solve` uses this method when $2k \ge n$.

### Output-sensitive termination {#termination}

Both modular solvers maintain a random linear combination $s$ of the
entries of the solution modulo the current modulus $m$ at negligible
cost, and periodically reconstruct $s$ as a fraction. Once that succeeds,
its denominator is (with high probability) the common denominator of the
solution and its numerator indicates the size of the solution. Only then
is the whole solution reconstructed, and it is accepted when

$$
n \lVert A \rVert \lVert X \rVert + \mathit{den} \lVert B \rVert < m,
$$

with maximum norms: since $AX - \mathit{den} \cdot B \equiv 0 \pmod m$ by
construction and the inequality bounds its entries by $m$ in absolute
value, it is zero. When this inequality would need a much larger modulus,
the solution is verified by computing $AX$ instead, if that is estimated
to be cheaper [@Cab1971; @CheSto2005; @Stef2010]. A wrong guess from the
random combination can therefore never produce a wrong result; it only
costs time. See [](algorithms/linear-algebra/multimodular#termination).

### Singular matrices {#singular}

Both modular solvers search for a prime modulo which $A$ has full rank.
A prime modulo which the rank drops is either a divisor of $\det(A) \ne
0$ (a bad prime) or a witness of singularity. The solvers try to prove
singularity with an exact kernel vector, found by solving the
nonsingular pivot block modulo $p$ against a non-pivot column and checking
the remaining rows exactly. If that fails, they try another prime, and
they also declare $A$ singular once the product of the primes with a rank
drop exceeds Hadamard's bound for $|\det(A)|$. Singular input is thus
detected after one or two primes in practice.

## Systems of any shape {#can-solve}

For an $m \times n$ matrix $A$ and an $m \times k$ matrix $B$,
`fmpz_mat_can_solve(X, den, A, B)` returns 1 if $AX = B$ has a solution
over $\Q$, and sets ($X$, `den`) to the solution whose entries in the
non-pivot rows (the rows of $X$ corresponding to columns of $A$ that are
not pivot columns of its echelon form) are zero. It returns 0 if the
system is inconsistent. For $m \le 15$ it calls `fmpz_mat_can_solve_fflu`,
otherwise the modular method of `fmpz_mat_can_solve_multi_mod_den` with
the pivot block solved by `fmpz_mat_solve`:

1. Reduce $A$ modulo a batch of primes above $2^{59}$ and compute the rank
   $r$ and the pivot rows and columns of $A$ modulo each, in parallel.
2. For the first prime of maximal rank, solve the $r \times r$ pivot
   block of $A$ against the corresponding rows of $B$ exactly.
3. Check the other $m - r$ equations exactly. If they hold, the solution
   is correct whatever the prime.
4. If an equation $i_0$ fails, try to prove inconsistency: solve the
   transposed pivot system for a vector $y$ supported on the pivot rows and
   row $i_0$ with $y^T A = 0$ on the pivot columns, and check exactly that
   $y^T A = 0$ on all columns. Since $y^T B \ne 0$, such a $y$ proves that
   no solution exists (Fredholm alternative) [@GieLobSau1998; @MulSto2004].
5. Otherwise the rank modulo the prime was too small. Continue with other
   primes, redoing exact work only when the rank increases, and declare the
   system inconsistent once the product of the primes exceeds a bound for
   all minors of $A$ (`fmpz_mat_det_bound_submatrix`), since then at least
   one prime gave the true rank.

`fmpz_mat_can_solve_dixon_den` and `fmpz_mat_can_solve_multi_mod_den`
force the solver used in step 2. `fmpz_mat_can_solve_fflu` uses
fraction-free elimination on $A$, solves with the pivot rows, and checks
the other equations. $X$ must be $n \times k$ and must not be aliased
with $A$ or $B$. On a return value of 0, $X$ and `den` are unspecified.

## Related functions {#related}

- `fmpz_mat_inv` solves against the identity matrix.
- `fmpz_mat_solve_dixon` returns the unreconstructed $p$-adic solution.
- `fmpq_mat_solve` and `fmpq_mat_can_solve` solve systems with rational
  entries by clearing denominators row by row and calling these
  functions, returning an `fmpq_mat_t`.
- `nmod_mat_solve` solves modulo a prime.

## Example {#example}

A $20 \times 20$ system with one right-hand side, solved by the default
method (here $p$-adic lifting) and by fraction-free elimination. The
$p$-adic solver returns the least common denominator, the fraction-free
solver the determinant, which is $-4000$ times larger. Then a $3 \times 4$
system of rank 2, consistent for one right-hand side and inconsistent for
another:

```c example name=fmpz_mat_solve_systems
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

int main(void)
{
    fmpz_mat_t A, B, X, Y, AX;
    fmpz_t den, den2;
    slong i, j, n = 20;

    fmpz_init(den);
    fmpz_init(den2);
    fmpz_mat_init(A, n, n);
    fmpz_mat_init(B, n, 1);
    fmpz_mat_init(X, n, 1);
    fmpz_mat_init(Y, n, 1);

    /* A[i][j] = (i + 2j)^2 mod 17 - 8, with 10 added on the diagonal */
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
            fmpz_set_si(fmpz_mat_entry(A, i, j),
                ((i + 2 * j) * (i + 2 * j)) % 17 - 8 + 10 * (i == j));
        fmpz_set_si(fmpz_mat_entry(B, i, 0), i + 1);
    }

    flint_printf("nonsingular: %d\n", fmpz_mat_solve(X, den, A, B));
    flint_printf("den = %{fmpz}\n", den);
    flint_printf("x_0 = %{fmpz} / %{fmpz}\n", fmpz_mat_entry(X, 0, 0), den);

    fmpz_mat_solve_fflu(Y, den2, A, B);
    flint_printf("fflu den = %{fmpz}\n", den2);
    fmpz_divexact(den, den2, den);
    flint_printf("ratio = %{fmpz}\n", den);

    /* check A X = den B */
    fmpz_mat_init(AX, n, 1);
    fmpz_mat_mul(AX, A, Y);
    fmpz_mat_scalar_mul_fmpz(B, B, den2);
    flint_printf("A Y = den B: %d\n", fmpz_mat_equal(AX, B));

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(X);
    fmpz_mat_clear(Y);
    fmpz_mat_clear(AX);

    /* a rank 2 system: the third row is the sum of the first two */
    fmpz_mat_init(A, 3, 4);
    fmpz_mat_init(B, 3, 1);
    fmpz_mat_init(X, 4, 1);
    for (j = 0; j < 4; j++)
    {
        fmpz_set_si(fmpz_mat_entry(A, 0, j), j + 1);
        fmpz_set_si(fmpz_mat_entry(A, 1, j), 3 * j * j - 2);
        fmpz_add(fmpz_mat_entry(A, 2, j),
                 fmpz_mat_entry(A, 0, j), fmpz_mat_entry(A, 1, j));
    }
    fmpz_set_si(fmpz_mat_entry(B, 0, 0), 1);
    fmpz_set_si(fmpz_mat_entry(B, 1, 0), 2);
    fmpz_set_si(fmpz_mat_entry(B, 2, 0), 3);

    flint_printf("consistent: %d\n", fmpz_mat_can_solve(X, den, A, B));
    flint_printf("X = ");
    fmpz_mat_print(X);
    flint_printf(", den = %{fmpz}\n", den);

    fmpz_set_si(fmpz_mat_entry(B, 2, 0), 4);
    flint_printf("consistent: %d\n", fmpz_mat_can_solve(X, den, A, B));

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(X);
    fmpz_clear(den);
    fmpz_clear(den2);
    return 0;
}
```
```output
nonsingular: 1
den = 30710924059369020030
x_0 = -15778958810838794440 / 30710924059369020030
fflu den = -122843696237476080120000
ratio = -4000
A Y = den B: 1
consistent: 1
X = 4 1  -3 4 0 0, den = 5
consistent: 0
```
