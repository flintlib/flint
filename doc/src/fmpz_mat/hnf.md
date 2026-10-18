---
title: Hermite normal form
functions: fmpz_mat_hnf fmpz_mat_hnf_modular fmpz_mat_hnf_transform
---

The Hermite normal form is the canonical basis of a lattice in $\Z^n$
given by generators, and the integer analogue of the reduced row echelon
form. Two integer matrices have the same Hermite normal form exactly when
their rows generate the same lattice. FLINT computes it by several
algorithms, from classical elimination with gcds to a randomised
algorithm of Pernet and Stein that reduces most of the work to linear
system solving; `fmpz_mat_hnf` chooses among them.

::: summary algorithms/linear-algebra/hermite-smith :::

## Definition {#definition}

An $m \times n$ integer matrix $H$ is in (row-style) Hermite normal form
if

- it is in row echelon form, with its $r$ nonzero rows first: the first
  nonzero entry (pivot) of row $i$ is in column $j_i$, with
  $j_0 < j_1 < \cdots < j_{r-1}$;
- the pivots are positive;
- the entries above each pivot are reduced: $0 \le H_{k j_i} < H_{i j_i}$
  for $k < i$.

For every integer matrix $A$ there is a unique $H$ in Hermite normal form
with $H = UA$ for some unimodular $U$ ($\det U = \pm 1$). Its nonzero
rows are a basis of the lattice spanned by the rows of $A$, the product
of its pivots is the determinant (covolume) of that lattice when $A$ has
full column rank, and the pivot columns are the column rank profile of
$A$. `fmpz_mat_is_in_hnf` checks the definition.

Column-style Hermite normal forms, $H = AV$ with the columns spanning the
lattice, are obtained by transposing: the column form of $A$ is the
transpose of `fmpz_mat_hnf` of $A^T$.

## Usage {#usage}

- $H$ must have the dimensions of $A$ and may be $A$.
- `fmpz_mat_hnf_minors` and the modular functions require $A$ to have
  full column rank (in particular $m \ge n$), which is not checked; the
  other functions accept any matrix.
- `fmpz_mat_hnf_modular(H, A, D)` requires $D$ to be a positive multiple
  of the determinant of the lattice, for example the absolute value of a
  nonzero $n \times n$ minor of $A$ when $A$ has full column rank, or
  $|\det A|$ for a nonsingular square $A$.
  `fmpz_mat_hnf_modular_eldiv(A, D)` transforms $A$ in place and requires
  $D$ to be a positive multiple of the largest elementary divisor of $A$.
- `fmpz_mat_hnf_transform(H, U, A)` also sets the $m \times m$
  unimodular $U$ with $UA = H$; $U$ is not unique when the rank is less
  than $m$.

## Choice of algorithm {#dispatch}

`fmpz_mat_hnf` uses `fmpz_mat_hnf_classical` when the number of rows $m$
is below a cutoff that decreases with the bit size $b$ of the largest
entry:

| $b$ | $\le 2$ | $\le 4$ | $\le 8$ | $\le 32$ | $\le 64$ | $\le 128$ | $\le 512$ | larger |
|---|---|---|---|---|---|---|---|---|
| cutoff | 52 | 38 | 30 | 11 | 5 | 4 | 3 | 2 |

Otherwise it calls `fmpz_mat_hnf_pernet_stein`. That function in turn
uses `fmpz_mat_hnf_minors` for matrices of full column rank with fewer
than about 35 to 52 columns (depending on $b$), and falls back to
`fmpz_mat_hnf_minors` or `fmpz_mat_hnf_classical` when the rank is too
small for its decomposition.

## Algorithms {#algorithms}

- `fmpz_mat_hnf_classical` eliminates column by column with the entry of
  smallest absolute value as pivot, by integer division with remainder,
  until the pivot divides the other entries of its column [@Coh1996,
  Algorithm 2.4.4]. Entries can grow exponentially in the worst case.
- `fmpz_mat_hnf_xgcd` combines pairs of rows with the cofactors of an
  extended gcd, which clears an entry in one step [@Coh1996, Algorithm
  2.4.5].
- `fmpz_mat_hnf_minors` is the algorithm of Kannan and Bachem
  [@KanBac1979]: it brings the leading $k \times k$ principal submatrix
  into Hermite normal form for $k = 1, 2, \ldots$, reducing the entries
  above the diagonal at each step, which keeps all entries polynomially
  bounded. It requires full column rank.
- `fmpz_mat_hnf_modular` is the algorithm of Domich, Kannan and Trotter
  [@DomKanTro1987; @Coh1996, Algorithm 2.4.8]: since $D\Z^n$ is contained
  in the lattice, all computations can be done modulo $D$, which bounds
  the entries by $D$.
- `fmpz_mat_hnf_modular_eldiv` computes the strong echelon form of $A$
  modulo $D$ and puts $D$ on zero diagonal entries; the lattice contains
  $D\Z^n$ when $D$ is a multiple of the largest elementary divisor
  [@FieHof2014].
- `fmpz_mat_hnf_pernet_stein` implements [@PernetStein2010]. It chooses
  $r - 1$ rows forming, together with one of two further rows, a
  nonsingular submatrix modulo a random prime, computes the two
  determinants obtained with each of the two further rows by solving one
  linear system over $\Q$, and combines the two rows with the cofactors of
  the gcd $g$ of these determinants. The Hermite normal form of the
  resulting $(r - 1) \times (r - 1)$ block is computed modulo $g$, which is
  small for random matrices; then the remaining columns are added by
  solving a linear system, and the remaining rows are added one at a time
  by an incremental update. A random prime giving wrong pivots is detected
  and the computation restarted.

See [](algorithms/linear-algebra/hermite-smith) for the mathematics and
the complexities.

## Transformation matrix {#transform}

`fmpz_mat_hnf_transform` uses `fmpz_mat_hnf_minors_transform` when
$m \ge n$ and $A$ has full column rank modulo a random prime; that
function performs the Kannan–Bachem elimination on $A$ and applies the
same operations to the identity matrix. Otherwise it computes the Hermite
normal form of the $m \times (n + m)$ matrix $(A \mid I_m)$, which has
full row rank, and reads off $H$ from the first $n$ columns and $U$ from
the last $m$.

## Example {#example}

The lattice spanned by three vectors in $\Z^3$, its Hermite normal form
computed by three algorithms, the determinant of the lattice as the
product of the pivots, and a transformation matrix:

```c example name=fmpz_mat_hnf_lattice
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

int main(void)
{
    fmpz_mat_t A, H, H2, U, UA;
    fmpz_t d;
    slong i, j;
    slong rows[3][3] = { { 4, 6, 10 }, { 2, -2, 8 }, { 6, 14, -6 } };

    fmpz_init(d);
    fmpz_mat_init(A, 3, 3);
    fmpz_mat_init(H, 3, 3);
    fmpz_mat_init(H2, 3, 3);
    fmpz_mat_init(U, 3, 3);
    fmpz_mat_init(UA, 3, 3);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            fmpz_set_si(fmpz_mat_entry(A, i, j), rows[i][j]);

    fmpz_mat_hnf(H, A);
    fmpz_mat_print_pretty(H);
    flint_printf("\nin HNF: %d\n", fmpz_mat_is_in_hnf(H));

    fmpz_mat_hnf_classical(H2, A);
    flint_printf("classical agrees: %d\n", fmpz_mat_equal(H, H2));

    /* the lattice contains |det A| Z^3 */
    fmpz_mat_det(d, A);
    fmpz_abs(d, d);
    flint_printf("|det A| = %{fmpz}\n", d);
    fmpz_mat_hnf_modular(H2, A, d);
    flint_printf("modular agrees: %d\n", fmpz_mat_equal(H, H2));

    fmpz_mat_hnf_transform(H2, U, A);
    fmpz_mat_mul(UA, U, A);
    flint_printf("U A = H: %d\n", fmpz_mat_equal(UA, H));
    fmpz_mat_det(d, U);
    flint_printf("det U = %{fmpz}\n", d);

    fmpz_mat_clear(A);
    fmpz_mat_clear(H);
    fmpz_mat_clear(H2);
    fmpz_mat_clear(U);
    fmpz_mat_clear(UA);
    fmpz_clear(d);
    return 0;
}
```
```output
[[2 8 2]
[0 10 12]
[0 0 18]
]
in HNF: 1
classical agrees: 1
|det A| = 360
modular agrees: 1
U A = H: 1
det U = 1
```

The product of the pivots, $2 \cdot 10 \cdot 18 = 360$, is the index of
the lattice in $\Z^3$, which equals $|\det A|$.
