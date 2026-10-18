---
title: Computing determinants, ranks and inverses
functions: nmod_mat_det nmod_mat_det_howell nmod_mat_rank nmod_mat_inv
---

The determinant can be computed for any modulus, prime or composite. The
rank and the inverse are defined through Gaussian elimination over a
field and require a prime modulus. All four functions work on a copy of
the input.

::: summary algorithms/linear-algebra/gaussian-elimination :::

## Determinant {#det}

`nmod_mat_det` chooses the method from the dimension $r$:

| $r$ | Method | Divisions |
|---|---|---|
| 0 | returns 1 (0 if $n = 1$) | none |
| 1 to 4 | cofactor expansion | none |
| 5 to 8 | characteristic polynomial by Berkowitz's algorithm, $\det A = (-1)^r \chi_A(0)$ | none |
| 9 and more, $n$ prime | `nmod_mat_lu` with `rank_check` set: the product of the pivots, times the sign of the row permutation | pivots |
| 9 and more, $n$ composite | `nmod_mat_det_howell` | none needed |

The primality of $n$ is decided by `n_is_prime` on every call. The cost
for large $r$ is that of an LU decomposition, $O(r^\omega)$ with the
exponent $\omega$ of matrix multiplication, for a prime modulus, and
$O(r^3)$ operations and $O(r^2)$ extended gcd computations for a
composite one.

`nmod_mat_det_howell` reduces $A$ to upper triangular form with row
operations that never divide. To clear an entry $b$ below a pivot $a$ it
subtracts a multiple of the pivot row if $a$ divides $b$ modulo $n$, and
otherwise first replaces the pivot row by $s$ times itself minus $t$
times the other row, where $sa - tb = \gcd(a, b)$ and $s$ is chosen
invertible modulo $n$. That operation multiplies the determinant by the
unit $s$, which is divided out at the end. Row swaps change the sign.
Despite its name, the function does not compute the Howell form.

## Rank {#rank}

`nmod_mat_rank` returns the number of pivots found by `nmod_mat_lu`
(without rank check) on a copy of $A$, which may have any shape. For a
composite modulus the rank is not well defined; the function then
raises a `FLINT_IMPINV` exception or returns a meaningless value.

## Inverse {#inv}

`nmod_mat_inv(B, A)` solves $AB = I$ with `nmod_mat_solve`, so its cost is
that of an LU decomposition and two triangular solves with $r$
right-hand sides, $O(r^\omega)$. A $1 \times 1$ matrix is inverted
directly. It returns 0 if $A$ is singular. $B$ may be the same matrix as
$A$.

For a composite modulus, a matrix can be invertible (its determinant is
a unit) while elimination meets a pivot that is a zero divisor, for
example $\begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$ modulo 4. The
function then raises a `FLINT_IMPINV` exception. A $1 \times 1$ matrix
with a nonzero, non-invertible entry also raises an exception.

## Example {#example}

The determinant of a $12 \times 12$ matrix modulo a prime and modulo
$10^6$, which agree with the exact integer determinant
696989974644262065485477877973 reduced modulo these numbers, and an
inverse and a rank modulo 11:

```c example name=nmod_mat_det_inv_rank
#include <flint/flint.h>
#include <flint/nmod_mat.h>

int main(void)
{
    nmod_mat_t A, B, I;
    slong i, j, n = 12;

    /* the same integer matrix modulo a prime and a composite number */
    ulong moduli[2] = { 1000003, 1000000 };

    for (i = 0; i < 2; i++)
    {
        nmod_mat_init(A, n, n, moduli[i]);
        /* entry j (in row-major order) is (j^2 + 3j + 1) mod 1000,
           plus 1 on the diagonal */
        for (j = 0; j < n * n; j++)
            nmod_mat_entry(A, j / n, j % n) =
                (j * j + 3 * j + 1) % 1000 + (j % (n + 1) == 0);

        flint_printf("n = %wu: det = %wu, det_howell = %wu\n", moduli[i],
                     nmod_mat_det(A), nmod_mat_det_howell(A));
        nmod_mat_clear(A);
    }

    /* inverse and rank modulo a prime */
    nmod_mat_init(A, 3, 3, 11);
    nmod_mat_init(B, 3, 3, 11);
    nmod_mat_init(I, 3, 3, 11);
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            nmod_mat_entry(A, i, j) = (i <= j) ? i + j + 1 : 0;

    flint_printf("invertible: %d\n", nmod_mat_inv(B, A));
    nmod_mat_print_pretty(B);
    nmod_mat_mul(I, A, B);
    flint_printf("A * A^-1 is one: %d\n", nmod_mat_is_one(I));

    nmod_mat_entry(A, 2, 2) = 0;
    flint_printf("rank %wd, det %wu, invertible: %d\n",
                 nmod_mat_rank(A), nmod_mat_det(A), nmod_mat_inv(B, A));

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(I);
    return 0;
}
```
```output
n = 1000003: det = 920755, det_howell = 920755
n = 1000000: det = 877973, det_howell = 877973
invertible: 1
<3 x 3 integer matrix mod 11>
[ 1  3  8]
[ 0  4 10]
[ 0  0  9]
A * A^-1 is one: 1
rank 2, det 0, invertible: 0
```
