---
title: Smith normal form and elementary divisors
functions: fmpz_mat_snf fmpz_mat_snf_transform fmpz_mat_elementary_divisors
---

The Smith normal form describes an integer matrix up to unimodular
transformations on both sides. If the rows of the $m \times n$ matrix $A$
generate a subgroup $L$ of $\Z^n$, the diagonal entries
$d_1 \mid d_2 \mid \cdots \mid d_r$ of its Smith form give the structure
of the quotient,

$$
\Z^n / L \cong \Z/d_1\Z \times \cdots \times \Z/d_r\Z \times \Z^{n - r},
$$

which is how finitely presented abelian groups, class groups and
homology groups are computed. FLINT computes the Smith form by
alternating Hermite normal forms, optionally with the transformation
matrices, and the elementary divisors alone by a $p$-adic method.

::: summary algorithms/linear-algebra/hermite-smith :::

## Definition {#definition}

For every $m \times n$ integer matrix $A$ of rank $r$ there are
unimodular matrices $U$ ($m \times m$) and $V$ ($n \times n$) such that
$S = UAV$ is diagonal with diagonal entries
$d_1, \ldots, d_r, 0, \ldots, 0$, where $d_i > 0$ and $d_i \mid d_{i+1}$.
The $d_i$, the *elementary divisors* or *invariant factors* of $A$, are
unique: $d_1 d_2 \cdots d_k$ is the gcd of the $k \times k$ minors of
$A$. $S$ is the Smith normal form of $A$; $U$ and $V$ are not unique.
`fmpz_mat_is_in_snf` checks this definition.

## Usage {#usage}

- `fmpz_mat_snf(S, A)` sets $S$, of the dimensions of $A$, to its Smith
  form. $S$ may be $A$.
- `fmpz_mat_snf_transform(S, U, V, A)` also computes $U$ and $V$; either
  may be `NULL`. $U$ and $V$ must be distinct from each other and from
  $A$ and $S$.
- `fmpz_mat_elementary_divisors(ed, A)` writes $d_1, \ldots, d_r$ to the
  vector `ed`, which needs room for $\min(m, n)$ entries, and returns $r$.
- `fmpz_mat_snf_kannan_bachem` and `fmpz_mat_snf_iliopoulos` are
  alternative algorithms; the latter needs a nonsingular square matrix and
  a multiple of its largest elementary divisor, for example
  $|\det A|$, and computes modulo that number.
- `fmpz_mat_snf_diagonal` computes the Smith form of a diagonal matrix.

## Algorithms {#algorithms}

`fmpz_mat_snf` returns at once with `fmpz_mat_snf_diagonal` if $A$ is
diagonal, which repeatedly replaces pairs of diagonal entries by their
gcd and lcm. Otherwise it alternates row-style Hermite normal forms of
the matrix and of its transpose (`fmpz_mat_hnf`, hence the Pernet–Stein
algorithm for large matrices) until the matrix is diagonal, then makes the
diagonal a divisibility chain by gcd and lcm steps and makes it
nonnegative. Each Hermite form is upper triangular and the next one of
the transpose clears the entries above the diagonal while possibly
creating new ones below; the number of rounds is small in practice, and
Kannan and Bachem showed that it is polynomial [@KanBac1979]. A
defensive bound (`_fmpz_mat_snf_iter_bound`) raises an exception if it is
ever exceeded.

`fmpz_mat_snf_transform` runs the same iteration with
`fmpz_mat_hnf_transform`, accumulating the transforms in $U$ and $V$. The
divisibility chain is then fixed with extended gcds: for a pair of
diagonal entries $(a, b)$ with $g = \gcd(a, b) = sa + tb$, a few row and
column operations, applied to $U$ and $V$ as well, turn $\operatorname{diag}
(a, b)$ into $\operatorname{diag}(g, ab/g)$. The construction follows the
implementation by Jean Kieffer in `acb_theta`.

`fmpz_mat_snf_kannan_bachem` eliminates the first row and column
alternately with extended gcd operations until the leading entry divides
all others [@KanBac1979]. `fmpz_mat_snf_iliopoulos` does the same
eliminations modulo a multiple of the largest elementary divisor, which
bounds the entries [@Iliopoulos1989].

`fmpz_mat_elementary_divisors` computes the Hermite normal form $H$ of $A$
and factors its pivots: the product of the pivots is a multiple of
$d_1 \cdots d_r$, so the elementary divisors are products of powers of
the primes dividing the pivots. For each such prime $p$, the valuations
$v_p(d_i)$ follow from the ranks modulo $p$ of a sequence of lattices
obtained by Lübeck's algorithm [@Lüb2002]: as long as the $r$ nonzero
rows $R$ of $H$ have rank less than $r$ modulo $p$, each vector $y$ of a
basis of the left nullspace of $R$ modulo $p$ gives the integral vector
$yR/p$, which replaces a row of $R$; the nullities at the successive
steps count the $d_i$ with $v_p(d_i) \ge 1, 2, \ldots$. If a pivot has
more than $2\beta$ bits or a prime factor that does not fit in a word,
the function computes the full Smith form instead.

## Example {#example}

The quotient of $\Z^3$ by the lattice of [](fmpz_mat/hnf#example), and
the same computation with transformation matrices:

```c example name=fmpz_mat_snf_group
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_mat.h>

int main(void)
{
    fmpz_mat_t A, S, U, V, T;
    fmpz * ed;
    slong i, j, r;
    slong rows[3][3] = { { 4, 6, 10 }, { 2, -2, 8 }, { 6, 14, -6 } };

    fmpz_mat_init(A, 3, 3);
    fmpz_mat_init(S, 3, 3);
    fmpz_mat_init(U, 3, 3);
    fmpz_mat_init(V, 3, 3);
    fmpz_mat_init(T, 3, 3);
    ed = _fmpz_vec_init(3);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            fmpz_set_si(fmpz_mat_entry(A, i, j), rows[i][j]);

    fmpz_mat_snf(S, A);
    fmpz_mat_print_pretty(S);
    flint_printf("\nin SNF: %d\n", fmpz_mat_is_in_snf(S));

    r = fmpz_mat_elementary_divisors(ed, A);
    flint_printf("Z^3 / L = ");
    for (i = 0; i < r; i++)
        flint_printf("%sZ/%{fmpz}", i ? " x " : "", ed + i);
    flint_printf("\n");

    fmpz_mat_snf_transform(T, U, V, A);
    fmpz_mat_mul(T, U, A);
    fmpz_mat_mul(T, T, V);
    flint_printf("U A V = S: %d\n", fmpz_mat_equal(T, S));

    fmpz_mat_clear(A);
    fmpz_mat_clear(S);
    fmpz_mat_clear(U);
    fmpz_mat_clear(V);
    fmpz_mat_clear(T);
    _fmpz_vec_clear(ed, 3);
    return 0;
}
```
```output
[[2 0 0]
[0 2 0]
[0 0 90]
]
in SNF: 1
Z^3 / L = Z/2 x Z/2 x Z/90
U A V = S: 1
```
