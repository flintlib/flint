---
title: Finding the complex roots of a polynomial
functions: acb_poly_find_roots _acb_poly_validate_roots _acb_poly_root_inclusion
---

`acb_poly_find_roots(roots, poly, initial, maxiter, prec)` computes
enclosures for all $n$ complex roots of a polynomial of degree $n$ with
complex ball coefficients and returns the number of them that it has
*isolated*. The enclosures are certified: if the return value is $n$,
each of the $n$ boxes contains exactly one root of every polynomial in
the ball `poly`, so all roots have been found. If it is smaller, the
first return-value boxes still contain exactly one root each, and the
others contain at least one root each but overlap, so that the roots
inside them are not separated. The numerical part of the computation
uses floating-point arithmetic; only the final validation uses ball
arithmetic, which is why the method is fast and why the user has to
retry with a higher precision when the validation fails.

::: summary algorithms/ball-arithmetic/polynomial-roots :::

## Algorithm {#algorithm}

1. Degenerate cases. For degree 0 nothing is done (the return value is
   0). If the leading coefficient contains zero, the roots can be
   anywhere, so all outputs are set to $[\pm \infty] + [\pm \infty]i$
   and 0 is returned. For degree 1 the root $-a_0/a_1$ is computed in
   ball arithmetic.
2. Double precision. If $\mathit{prec} \le 53$, the polynomial is
   rescaled ($x \to \lambda x$ and a constant factor, with $\lambda$
   from the line that best fits the points $(i, \log_2 |a_i|)$ in the
   minimax sense) so that its coefficients and roots are representable
   as doubles, and `cd_poly_find_roots` runs up to $\max(150,
   \mathit{maxiter})$ damped Weierstrass–Durand–Kerner iterations in
   double precision. If the last relative correction is below $1/8$,
   the result is passed to the validation.
3. Multiprecision. Otherwise the starting values are `initial`, or, if
   `initial` is `NULL`, points on circles whose radii come from the
   Newton polygon of the coefficients [@Bin1996]: if $k_1 < \cdots <
   k_q$ are the vertices of the upper convex hull of the points $(k,
   \log_2 |a_k|)$, then $m_i = k_i - k_{i-1}$ points are placed on the
   circle of radius $|a_{k_{i-1}} / a_{k_i}|^{1/m_i}$, at angles that
   avoid symmetric configurations. The iteration
   $$z_k \leftarrow z_k - \frac{f(z_k)}{a_n \prod_{j \ne k} (z_k - z_j)}$$
   (`_gr_poly_refine_roots_wdk`) then runs on floating-point numbers of
   `prec` bits (the `nfloat_complex` context, or `acf` if `nfloat` does
   not support the precision), for at most `maxiter` iterations, the
   default being $2n + \lfloor\sqrt{\mathit{prec}}\rfloor$. The iteration
   stops two, three or four iterations after the largest correction,
   relative to the largest root, falls below $2^{-\mathit{prec}/2}$,
   $2^{-\mathit{prec}/3}$ or $2^{-\mathit{prec}/4}$, since it converges
   quadratically from then on.
4. Validation. `_acb_poly_validate_roots` replaces each approximation
   $m$ by the box around the disk of radius $n |f(m) / f'(m)|$, which
   contains a root ([@Kob2010], see `_acb_poly_root_inclusion`), sorts
   the boxes that overlap no other box to the front and returns their
   number. If $n$ disjoint disks each contain at least one root, each
   contains exactly one.

The iteration does not depend on the coefficient radii, and the
validation uses upper bounds of $|f(m)|$ and $1/|f'(m)|$ over all
polynomials in the coefficient balls, so the enclosures are valid for
every polynomial in `poly`. The method converges for squarefree
polynomials and generic starting values; each iteration costs $O(n^2)$
operations at `prec` bits. Since the validation computes $f(m)$ with
cancellation, the attainable accuracy for a root $\zeta$ is limited by
the condition number of the root: the precision has to exceed the
number of bits lost in evaluating $f$ near $\zeta$ by a margin.

## Using the function {#usage}

The polynomial must be squarefree: at a root of multiplicity $m > 1$,
the iteration converges only linearly to about $\mathit{prec}/m$ bits,
and the inclusion disks of the $m$ approximations overlap for every
precision. Remove repeated roots first, for example with
`fmpz_poly_factor_squarefree` for integer polynomials. Clusters of
close roots and roots of very different magnitudes need a higher
precision before they can be separated.

The recommended use is a loop that doubles the precision until all
roots are isolated, reusing the previous approximations as `initial`
values; `arb_fmpz_poly_complex_roots` implements this, together with
real root isolation and deflation, for integer polynomials and is the
function to use in that case. For a real polynomial,
`acb_poly_validate_real_roots` proves which of the isolated roots are
real, after which their imaginary parts can be set to zero.

The output array must have room for $n$ entries. The function throws an
exception if `poly` is the zero polynomial. The underscore version
takes the coefficient array and the length $n + 1$.

## Examples {#examples}

The program finds the five roots of $x^5 - x - 1$, proves that exactly
one of them is real, and shows what happens at a double root: the
simple root $-2$ of $x^3 - 3x + 2 = (x-1)^2(x+2)$ is isolated, while
the double root 1 is found only to about half the precision, in two
overlapping boxes.

```c example name=acb_poly_find_roots_basic
#include <flint/acb_poly.h>

static void show(acb_srcptr roots, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
    {
        acb_printn(roots + i, 15, 0);
        flint_printf("\n");
    }
}

int main(void)
{
    acb_poly_t f;
    acb_ptr roots;
    slong i, isolated;

    acb_poly_init(f);
    roots = _acb_vec_init(5);

    /* x^5 - x - 1 */
    acb_poly_set_coeff_si(f, 5, 1);
    acb_poly_set_coeff_si(f, 1, -1);
    acb_poly_set_coeff_si(f, 0, -1);

    isolated = acb_poly_find_roots(roots, f, NULL, 0, 64);
    flint_printf("x^5 - x - 1: %wd isolated\n", isolated);
    if (isolated == 5 && acb_poly_validate_real_roots(roots, f, 64))
    {
        for (i = 0; i < 5; i++)
            if (arb_contains_zero(acb_imagref(roots + i)))
                arb_zero(acb_imagref(roots + i));
    }
    _acb_vec_sort_pretty(roots, 5);
    show(roots, 5);

    /* (x - 1)^2 (x + 2) = x^3 - 3x + 2 has a double root */
    acb_poly_zero(f);
    acb_poly_set_coeff_si(f, 3, 1);
    acb_poly_set_coeff_si(f, 1, -3);
    acb_poly_set_coeff_si(f, 0, 2);
    isolated = acb_poly_find_roots(roots, f, NULL, 0, 64);
    flint_printf("x^3 - 3x + 2: %wd isolated\n", isolated);
    show(roots, 3);

    _acb_vec_clear(roots, 5);
    acb_poly_clear(f);
    return 0;
}
```
```output
x^5 - x - 1: 5 isolated
[1.16730397826142 +/- 1.32e-15]
[-0.764884433600585 +/- 2.75e-16] + [0.352471546031726 +/- 2.51e-16]*I
[-0.764884433600585 +/- 2.75e-16] + [-0.352471546031726 +/- 2.51e-16]*I
[0.181232444469875 +/- 3.85e-16] + [-1.08395410131771 +/- 6.70e-16]*I
[0.181232444469875 +/- 3.86e-16] + [1.08395410131771 +/- 6.70e-16]*I
x^3 - 3x + 2: 1 isolated
[-2.00000000000000 +/- 1.30e-18] + [+/- 1.08e-18]*I
[1.00000 +/- 9.87e-6] + [+/- 1.17e-5]*I
[1.0000 +/- 6.10e-6] + [+/- 7.18e-6]*I
```

The radii printed for the first polynomial are dominated by the
conversion to 15 decimal digits; the enclosures computed at 64 bits are
accurate to about 59 bits.

A polynomial whose roots are hard to separate needs a precision loop.
Wilkinson's polynomial $\prod_{k=1}^{20} (x - k)$ has integer
coefficients of up to 62 bits and roots that are extremely sensitive to
perturbations of the coefficients:

```c example name=acb_poly_find_roots_wilkinson
#include <flint/acb_poly.h>

int main(void)
{
    acb_poly_t f;
    acb_ptr roots, xs;
    slong i, isolated, prec, n = 20;

    acb_poly_init(f);
    roots = _acb_vec_init(n);
    xs = _acb_vec_init(n);

    for (i = 0; i < n; i++)
        acb_set_si(xs + i, i + 1);
    acb_poly_product_roots(f, xs, n, 128);     /* exact */

    for (prec = 32; ; prec *= 2)
    {
        isolated = acb_poly_find_roots(roots, f, (prec == 32) ? NULL : roots,
            0, prec);
        flint_printf("prec %wd: %wd isolated\n", prec, isolated);
        if (isolated == n)
            break;
    }

    _acb_vec_sort_pretty(roots, n);
    acb_printn(roots + 0, 10, 0); flint_printf("\n");
    acb_printn(roots + 19, 10, 0); flint_printf("\n");

    _acb_vec_clear(roots, n);
    _acb_vec_clear(xs, n);
    acb_poly_clear(f);
    return 0;
}
```
```output radius-tol
prec 32: 0 isolated
prec 64: 20 isolated
[1.000000000 +/- 7.32e-16] + [+/- 7.24e-16]*I
[20.00000 +/- 4.55e-6] + [+/- 4.49e-6]*I
```

No root can be isolated at 32 bits. At 64 bits all of them are, but the
enclosure of the root 20 has radius $5 \cdot 10^{-6}$ while that of the
root 1 has radius $10^{-15}$: a relative perturbation $\varepsilon$ of
the coefficients can move the root 20 by about $3 \cdot 10^{12} \varepsilon$
(and the root 1 by about $420 \varepsilon$), and
the validation sees the corresponding rounding errors in $f(m)$.

The radii, and the last digits of the midpoints, can differ slightly
between platforms, since the starting values come from the
double-precision iteration at 32 bits (step 2 above), whose rounding
depends on the compiler. Every enclosure is still rigorous.
