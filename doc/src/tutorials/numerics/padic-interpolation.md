---
title: Interpolation over the p-adic numbers: precision and Mahler's basis
summary: Recovering a polynomial from its values at 0, 1, ..., n-1 when the values are only known modulo p^N: how many p-adic digits of each coefficient are determined (a closed formula with factorials and Stirling numbers), how FLINT's p-adic linear algebra attains exactly this precision, and how much the LU solver and the Newton and Lagrange interpolation formulas lose.
order: 10
---

The previous tutorials compute with real balls $[m \pm r]$. The $p$-adic
numbers have their own balls: in [](padic_radix) an element is
$p^v u + O(p^N)$, a residue known modulo $p^N$, and every operation
propagates these error terms ([](algorithms/p-adic/precision)). As with
real balls, the guaranteed precision of a result depends on the
algorithm, and the interesting question is how much precision is really
determined by the data. This tutorial answers it for a classical
problem, where the answer has a pleasant closed form.

## The problem {#problem}

Let $f \in \mathbf{Z}_p[x]$ have degree less than $n$, and suppose that
we know its values $f(0), f(1), \ldots, f(n-1)$ modulo $p^N$, for
example because they were computed $p$-adically. The coefficients
$c_0, \ldots, c_{n-1}$ of $f$ are the solution of the Vandermonde system
$$
\sum_{k=0}^{n-1} c_k\, i^k = f(i), \qquad i = 0, \ldots, n - 1.
$$
The Vandermonde determinant $\prod_{i < j}(j - i) = \prod_{m < n} m!$ is
divisible by a large power of $p$ as soon as $n > p$: the nodes
$0, \ldots, n-1$ collide modulo $p$. So some precision is lost. How much?

## Mahler's basis {#mahler}

Write $f$ in the binomial basis,
$$
f(x) = \sum_{m=0}^{n-1} \Delta^m f(0) \binom{x}{m}, \qquad
\Delta^m f(0) = \sum_{j=0}^{m} (-1)^{m-j} \binom{m}{j} f(j),
$$
which is the finite form of Mahler's expansion of continuous functions
on $\mathbf{Z}_p$ [@Mah1958]. The map from the values $f(0), \ldots,
f(n-1)$ to the differences $\Delta^m f(0)$ is given by a triangular
integer matrix with ones on the diagonal, and so is its inverse
($f(j) = \sum_m \binom{j}{m} \Delta^m f(0)$). It maps the set of value
vectors known modulo $p^N$ onto the set of difference vectors known
modulo $p^N$: no precision is lost or gained.

All the loss is in the change of basis to the monomials. With the
Stirling numbers of the first kind, $\binom{x}{m} = \frac{1}{m!}
\sum_{k} s(m, k)\, x^k$, so
$$
c_k = \sum_{m=k}^{n-1} \frac{s(m, k)}{m!}\, \Delta^m f(0).
$$
Each difference is an independent unknown modulo $p^N$, so $c_k$ is
determined exactly modulo $p^{N - \ell_k}$ with
$$
\ell_k = \max_{k \le m < n,\ s(m, k) \ne 0} \big(\operatorname{ord}_p(m!) -
\operatorname{ord}_p s(m, k)\big),
$$
and not to any higher precision. By Legendre's formula
$\operatorname{ord}_p(m!) = (m - s_p(m))/(p - 1)$, where $s_p(m)$ is the
sum of the base-$p$ digits of $m$, the loss grows like $n/(p - 1)$. For
example $c_1$ loses $\max_{m < n} \operatorname{ord}_p(m)$ digits, since
$s(m, 1) = \pm (m - 1)!$, and $c_0 = f(0)$ loses nothing.

## The program {#program}

The program interpolates $f(x) = 1 + 2x + 3x^2 + \cdots + n x^{n-1}$
over $\mathbf{Q}_p$ from values known modulo $p^N$ in four ways:

- with `gr_mat_nonsingular_solve`, which for `padic_radix` matrices uses
  Gauss–Jordan elimination with full pivoting on valuations
  ([](algorithms/p-adic/linear-algebra#solve));
- with `gr_mat_nonsingular_solve_lu`, an LU decomposition with pivots of
  minimal valuation in each column followed by back substitution;
- with `_gr_poly_interpolate_newton`, which computes divided differences
  and converts from the Newton basis to the monomial basis;
- with `_gr_poly_interpolate_fast`, the Lagrange formula
  $f = \sum_i y_i w_i \prod_{j \ne i} (x - x_j)$ with weights
  $w_i = 1 / \prod_{j \ne i} (x_i - x_j)$, evaluated with a subproduct
  tree.

It prints the absolute precision of each coefficient for all four, the
precision $N - \ell_k$ predicted above, and whether all four results
agree with $k + 1$ to their claimed precisions.

A few points about the code:

- The values are made inexact by adding the element $0 + O(p^N)$, which
  is parsed from a string by `gr_set_str`.
- The context has absolute precision $3N$. This is the precision to
  which exact intermediate quantities (here, everything computed from the
  Vandermonde matrix alone) are rounded; it only needs to be well beyond
  the uncertainty of the data. With both context precisions infinite the
  arithmetic would be exact, and the division of an exact element by
  $3$ would fail.
- `PADIC_RADIX_N` reads the precision $N$ of an element and
  `arith_stirling_number_1` gives $s(m, k)$.

```c example name=tutorial_padic_interpolation requires=64
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/arith.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/gr_poly.h>
#include <flint/padic_radix.h>

/* The p-adic valuation of a nonzero integer. */
static slong val(const fmpz_t x, ulong p)
{
    fmpz_t t, q;
    slong v;
    fmpz_init(t);
    fmpz_init_set_ui(q, p);
    v = fmpz_remove(t, x, q);
    fmpz_clear(t);
    fmpz_clear(q);
    return v;
}

/* Whether x agrees with the integer c to the precision of x. */
static int agrees(gr_srcptr x, slong c, gr_ctx_t Qp)
{
    gr_ptr d = gr_heap_init(Qp);
    int ok = (gr_sub_si(d, x, c, Qp) == GR_SUCCESS && gr_is_zero(d, Qp) != T_FALSE);
    gr_heap_clear(d, Qp);
    return ok;
}

/* Interpolate f(x) = 1 + 2x + ... + n x^(n-1) from f(0), ..., f(n-1),
   known modulo p^N, over Q_p. */
static void interpolate(ulong p, slong n, slong N)
{
    gr_ctx_t Qp;
    gr_mat_t V, y, c, cLU;
    padic_radix_struct * xs, * ys, * cN, * cL;
    gr_ptr err;
    fmpz_t t, s, fac;
    char str[64];
    slong i, k, m;
    int status = GR_SUCCESS;

    /* exact intermediate quantities are rounded to 3N digits, far beyond
       the errors of the data */
    gr_ctx_init_padic_radix(Qp, p, PADIC_RADIX_PREC_INF, 3 * N, 0);
    gr_mat_init(V, n, n, Qp);
    gr_mat_init(y, n, 1, Qp);
    gr_mat_init(c, n, 1, Qp);
    gr_mat_init(cLU, n, 1, Qp);
    err = gr_heap_init(Qp);
    xs = gr_heap_init_vec(n, Qp);
    ys = gr_heap_init_vec(n, Qp);
    cN = gr_heap_init_vec(n, Qp);
    cL = gr_heap_init_vec(n, Qp);
    fmpz_init(t);
    fmpz_init(s);
    fmpz_init(fac);

    /* the error term 0 + O(p^N) */
    flint_sprintf(str, "O(%wu^%wd)", p, N);
    status |= gr_set_str(err, str, Qp);

    for (i = 0; i < n; i++)
    {
        /* the node x_i = i, row i of the Vandermonde matrix, and
           y_i = f(i) + O(p^N) */
        fmpz_zero(s);
        for (k = n - 1; k >= 0; k--)
        {
            fmpz_mul_ui(s, s, i);
            fmpz_add_ui(s, s, k + 1);
        }
        status |= gr_add_fmpz(gr_mat_entry_ptr(y, i, 0, Qp), err, s, Qp);
        status |= gr_set(ys + i, gr_mat_entry_srcptr(y, i, 0, Qp), Qp);
        status |= gr_set_ui(xs + i, i, Qp);

        for (k = 0; k < n; k++)
        {
            fmpz_ui_pow_ui(t, i, k);
            status |= gr_set_fmpz(gr_mat_entry_ptr(V, i, k, Qp), t, Qp);
        }
    }

    status |= gr_mat_nonsingular_solve(c, V, y, Qp);
    status |= gr_mat_nonsingular_solve_lu(cLU, V, y, Qp);
    status |= _gr_poly_interpolate_newton(cN, xs, ys, n, Qp);
    status |= _gr_poly_interpolate_fast(cL, xs, ys, n, Qp);

    flint_printf("p = %wu, n = %wd, data known modulo p^%wd\n", p, n, N);
    flint_printf("  k  predicted  solve  LU  Newton  Lagrange  correct\n");

    for (k = 0; k < n; k++)
    {
        const padic_radix_struct * ck = gr_mat_entry_srcptr(c, k, 0, Qp);
        const padic_radix_struct * dk = gr_mat_entry_srcptr(cLU, k, 0, Qp);
        slong loss = 0;

        /* c_k = sum_m (Delta^m f)(0) s(m, k) / m!, and the map from the
           values to the Mahler coefficients (Delta^m f)(0) is unimodular */
        for (m = k; m < n; m++)
        {
            arith_stirling_number_1(s, m, k);
            fmpz_fac_ui(fac, m);
            if (!fmpz_is_zero(s))
                loss = FLINT_MAX(loss, val(fac, p) - val(s, p));
        }

        /* every result must agree with k + 1 to its claimed precision */
        flint_printf("%3wd  %9wd  %5wd  %2wd  %6wd  %8wd  %s\n", k, N - loss,
            PADIC_RADIX_N(ck), PADIC_RADIX_N(dk),
            PADIC_RADIX_N(cN + k), PADIC_RADIX_N(cL + k),
            (agrees(ck, k + 1, Qp) && agrees(dk, k + 1, Qp)
                && agrees(cN + k, k + 1, Qp) && agrees(cL + k, k + 1, Qp)) ? "yes" : "no");
    }

    if (status != GR_SUCCESS)
        flint_printf("unexpected status %d\n", status);

    gr_heap_clear(err, Qp);
    gr_heap_clear_vec(xs, n, Qp);
    gr_heap_clear_vec(ys, n, Qp);
    gr_heap_clear_vec(cN, n, Qp);
    gr_heap_clear_vec(cL, n, Qp);
    gr_mat_clear(V, Qp);
    gr_mat_clear(y, Qp);
    gr_mat_clear(c, Qp);
    gr_mat_clear(cLU, Qp);
    fmpz_clear(t);
    fmpz_clear(s);
    fmpz_clear(fac);
    gr_ctx_clear(Qp);
}

int main(void)
{
    interpolate(3, 12, 30);
    flint_printf("\n");
    interpolate(2, 10, 30);
    flint_cleanup_master();
    return 0;
}
```
```output
p = 3, n = 12, data known modulo p^30
  k  predicted  solve  LU  Newton  Lagrange  correct
  0         30     30  30      30        30  yes
  1         28     28  26      28        28  yes
  2         28     28  26      27        27  yes
  3         26     26  26      26        26  yes
  4         26     26  26      26        26  yes
  5         26     26  26      26        26  yes
  6         27     27  27      26        26  yes
  7         27     27  27      26        26  yes
  8         27     27  27      26        26  yes
  9         26     26  26      26        26  yes
 10         26     26  26      26        26  yes
 11         26     26  26      26        26  yes

p = 2, n = 10, data known modulo p^30
  k  predicted  solve  LU  Newton  Lagrange  correct
  0         30     30  30      30        30  yes
  1         27     27  23      27        27  yes
  2         25     25  23      25        25  yes
  3         25     25  23      24        24  yes
  4         23     23  23      23        23  yes
  5         23     23  23      23        23  yes
  6         24     24  24      23        23  yes
  7         24     24  24      23        23  yes
  8         23     23  23      23        23  yes
  9         23     23  23      23        23  yes
```

The full pivoting solver attains the predicted precision for every
coefficient, and all four methods are correct (their claims are always
valid). The LU solver claims up to two digits less for $p = 3$ and up to
four bits less for $p = 2$, for the low coefficients: in back
substitution, the coefficients of high degree are computed first, rounded
to their own precision, and then divided by the small pivots that belong
to the low coefficients, which amplifies their error terms. Gauss–Jordan
elimination with full pivoting performs the divisions before the rounded
quantities are combined, so that the error terms stay as small as the
data allow; see [](algorithms/p-adic/linear-algebra) for the theory and
for the same phenomenon with determinants and characteristic polynomials.

The interpolation formulas do better than LU on the low coefficients but
lose one digit (one bit for $p = 2$) on a few others, the same ones for
both. In exact arithmetic each coefficient is a linear combination
$c_k = \sum_i (V^{-1})_{ki}\, y_i$ of the values, and since the errors of
the values are independent, its optimal precision is
$N + \min_i \operatorname{ord}_p (V^{-1})_{ki}$, which is the formula
above. But neither method forms these combinations term by term. Newton's
method goes through the divided differences, here
$\Delta^m f(0) / m!$, and expands the Newton form
$d_0 + x(d_1 + (x - 1)(d_2 + \cdots))$ from the inside out. The fast
Lagrange formula sums $y_i w_i \prod_{j \ne i}(x - j)$ along a subproduct
tree. Either way, intermediate quantities carry the errors of several
values at once, rounded to the precision of the least precise one. The
extra divisibility by $p$ of the final coefficients (the Stirling numbers
$s(m, k)$ in the formula for $\ell_k$) comes from cancellations between
these intermediate quantities, which a bound computed one operation at a
time cannot see.

## Variations {#variations}

- With more nodes the loss grows: for $p = 3$ and $n = 30$ the
  coefficient of $x^{27}$ loses $\operatorname{ord}_3(27!) = 13$ digits.
  Interpolating at nodes that are distinct modulo $p$ (possible only for
  $n \le p$) loses nothing at all.
- If the nodes are only known modulo $p^N$ as well, the Vandermonde
  matrix becomes inexact and its own error terms contribute; the solver
  still returns the optimal precision $\min(N_B, N_A + \ldots)$ of
  [](algorithms/p-adic/linear-algebra#solve).
- The program computes in $\mathbf{Q}_p$ with a word-size $p$. For the
  fixed-precision modules [](padic) and [](qadic), and for the choice
  between the two models, see [](algorithms/p-adic/overview).
