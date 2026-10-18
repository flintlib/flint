---
title: Roots of polynomials and eigenvalues
functions: qqbar_roots_fmpz_poly qqbar_real_roots_fmpz_poly qqbar_eigenvalues_fmpz_mat QQBAR_ROOTS_IRREDUCIBLE
---

The roots of a polynomial with integer or rational coefficients are
algebraic numbers, and `qqbar_roots_fmpz_poly` returns all of them, with
multiplicity, as exact `qqbar_t` values. This is the usual way to create
algebraic numbers of higher degree, and the basis of exact eigenvalue
computations for integer and rational matrices.

::: summary algorithms/exact/qqbar-representation :::

## Method {#method}

Unless `QQBAR_ROOTS_IRREDUCIBLE` is given, the polynomial is factored
over $\Z$ with `fmpz_poly_factor` (see [](fmpz_poly_factor)). For each
irreducible factor $P$ of degree $d \ge 2$, the complex roots are
computed numerically with `arb_fmpz_poly_complex_roots` at 128 bits
([](arb_fmpz_poly/roots), [](algorithms/ball-arithmetic/polynomial-roots)),
doubling the precision until every ball is certified by an interval
Newton step to contain exactly one root of $P$; then each root is
stored as $P$, made primitive with positive leading coefficient,
together with its ball. A factor of degree 1 gives a rational root
directly. Roots of a factor with multiplicity $e$ are repeated $e$ times.
Finally the roots are sorted with `qqbar_cmp_root_order` unless
`QQBAR_ROOTS_UNSORTED` is given: real roots in decreasing order, then
nonreal roots by decreasing real part, increasing absolute imaginary
part and positive imaginary part first.

The cost is dominated by the factorisation over $\Z$ for polynomials
with many factors or large coefficients, and otherwise by numerical root
isolation, which is fast for degrees up to a few hundred.
`qqbar_real_roots_fmpz_poly` saves work by computing the squarefree
factorisation first, skipping squarefree parts without real roots
(counted with `fmpz_poly_num_real_roots`) before factoring the others,
and by isolating only real roots (`arb_fmpz_poly_real_roots`). `qqbar_eigenvalues_fmpz_mat` computes the characteristic
polynomial with `fmpz_mat_charpoly` and calls `qqbar_roots_fmpz_poly`.

For polynomials with algebraic coefficients, use `gr_poly_roots` with a
context created by `gr_ctx_init_complex_qqbar` (see [](gr_domains)),
which performs a squarefree factorisation and then calls
`_qqbar_roots_poly_squarefree`.

## Example {#example}

The polynomial $(x^3 - 2)(x^2 + x + 1)^2$ has the real root
$\sqrt[3]{2}$, the two primitive cube roots of unity with multiplicity
2, and the two other cube roots of 2. The eigenvalues of
$\left(\begin{smallmatrix} 1 & 2 \\ 3 & 4 \end{smallmatrix}\right)$ are
$(5 \pm \sqrt{33})/2$, printed with their minimal polynomial
$x^2 - 5x - 2$ and isolating ball.

```c example name=qqbar_roots
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/qqbar.h>

int main(void)
{
    fmpz_poly_t f;
    fmpz_mat_t A;
    qqbar_ptr r;
    slong i, n;

    /* f = (x^3 - 2) (x^2 + x + 1)^2 */
    fmpz_poly_init(f);
    fmpz_poly_set_str(f, "8  -2 -4 -6 -3 0 3 2 1");
    flint_printf("f = "); fmpz_poly_print_pretty(f, "x"); flint_printf("\n");

    r = _qqbar_vec_init(7);
    qqbar_roots_fmpz_poly(r, f, 0);
    for (i = 0; i < 7; i++)
    {
        qqbar_printnd(r + i, 10);
        flint_printf("\n");
    }

    n = qqbar_real_roots_fmpz_poly(r, f, 0);
    flint_printf("%wd real root(s)\n", n);

    /* eigenvalues of [[1, 2], [3, 4]] */
    fmpz_mat_init(A, 2, 2);
    fmpz_set_si(fmpz_mat_entry(A, 0, 0), 1);
    fmpz_set_si(fmpz_mat_entry(A, 0, 1), 2);
    fmpz_set_si(fmpz_mat_entry(A, 1, 0), 3);
    fmpz_set_si(fmpz_mat_entry(A, 1, 1), 4);
    qqbar_eigenvalues_fmpz_mat(r, A, 0);
    for (i = 0; i < 2; i++)
    {
        qqbar_print(r + i);
        flint_printf("\n");
    }

    _qqbar_vec_clear(r, 7);
    fmpz_mat_clear(A);
    fmpz_poly_clear(f);
    return 0;
}
```
```output
f = x^7+2*x^6+3*x^5-3*x^3-6*x^2-4*x-2
1.259921050 (deg 3)
-0.5000000000 + 0.8660254038*I (deg 2)
-0.5000000000 + 0.8660254038*I (deg 2)
-0.5000000000 - 0.8660254038*I (deg 2)
-0.5000000000 - 0.8660254038*I (deg 2)
-0.6299605249 + 1.091123636*I (deg 3)
-0.6299605249 - 1.091123636*I (deg 3)
1 real root(s)
deg 2 [-2, -5, 1] [5.3722813232690143299253057341094646591 +/- 4.84e-38]
deg 2 [-2, -5, 1] [-0.37228132326901432992530573410946465911 +/- 4.21e-39]
```

The array passed to `qqbar_roots_fmpz_poly` must have room for $\deg f$
entries; `qqbar_real_roots_fmpz_poly` writes only as many entries as it
returns.
