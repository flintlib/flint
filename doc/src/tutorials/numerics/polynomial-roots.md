---
title: Certified roots of polynomials
summary: All real and complex roots of integer polynomials with arb_fmpz_poly, multiplicities via squarefree factorisation, and roots of polynomials with inexact coefficients with acb_poly_find_roots and a precision loop; the example program poly_roots.
order: 7
---

Every polynomial of degree $n$ has $n$ complex roots, so a root finder
for polynomials can do something a general root finder cannot: certify
that it has found all of them. FLINT provides two levels of interface.
For polynomials with integer coefficients, [](arb_fmpz_poly) computes all
real or complex roots, isolated and to any requested accuracy, with no
choices left to the user. For polynomials whose coefficients are balls,
`acb_poly_find_roots` in [](acb_poly) computes approximations and
certifies them at a given precision, and the precision loop is left to
the caller.

## Integer polynomials {#integer}

`arb_fmpz_poly_complex_roots(roots, f, flags, prec)` writes the $n$
roots of a squarefree $f \in \Z[x]$ to an array of `acb_t`, each with at
least `prec` bits of relative accuracy, in disjoint balls. The real roots
come first, in increasing order and with imaginary parts exactly zero,
followed by the nonreal roots in conjugate pairs.
`arb_fmpz_poly_real_roots(roots, f, flags, prec)` computes only the real
roots and returns their number. Both require $f$ to be squarefree;
multiple roots are handled by computing a squarefree factorisation
$f = \prod_i f_i^{e_i}$ first, with `fmpz_poly_factor_squarefree`, and
finding the roots of each factor $f_i$, which are roots of $f$ of
multiplicity $e_i$.

```c example name=numerics_poly_roots_fmpz
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>

int main(void)
{
    fmpz_poly_t f, g;
    fmpz_poly_factor_t fac;
    acb_ptr roots;
    arb_ptr real;
    slong i, j, deg, nreal;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_factor_init(fac);

    /* f = (x^5 - x - 1) (x^2 - 2)^2 */
    fmpz_poly_set_str(f, "6  -1 -1 0 0 0 1");
    fmpz_poly_set_str(g, "3  -2 0 1");
    fmpz_poly_mul(f, f, g);
    fmpz_poly_mul(f, f, g);

    fmpz_poly_factor_squarefree(fac, f);

    for (i = 0; i < fac->num; i++)
    {
        deg = fmpz_poly_degree(fac->p + i);
        flint_printf("multiplicity %wd: ", fac->exp[i]);
        fmpz_poly_print_pretty(fac->p + i, "x");
        flint_printf("\n");

        /* the real roots only */
        real = _arb_vec_init(deg);
        nreal = arb_fmpz_poly_real_roots(real, fac->p + i, 0, 64);
        flint_printf("  %wd real roots\n", nreal);
        _arb_vec_clear(real, deg);

        /* all complex roots */
        roots = _acb_vec_init(deg);
        arb_fmpz_poly_complex_roots(roots, fac->p + i, 0, 64);
        for (j = 0; j < deg; j++)
        {
            flint_printf("  ");
            acb_printn(roots + j, 15, 0);
            flint_printf("\n");
        }
        _acb_vec_clear(roots, deg);
    }

    fmpz_poly_factor_clear(fac);
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    flint_cleanup();
    return 0;
}
```

```output
multiplicity 1: x^5-x-1
  1 real roots
  [1.16730397826142 +/- 1.32e-15]
  [-0.764884433600585 +/- 2.74e-16] + [0.352471546031726 +/- 2.50e-16]*I
  [-0.764884433600585 +/- 2.74e-16] + [-0.352471546031726 +/- 2.50e-16]*I
  [0.181232444469875 +/- 3.84e-16] + [1.08395410131771 +/- 6.69e-16]*I
  [0.181232444469875 +/- 3.84e-16] + [-1.08395410131771 +/- 6.69e-16]*I
multiplicity 2: x^2-2
  2 real roots
  [-1.41421356237310 +/- 4.96e-15]
  [1.41421356237310 +/- 4.96e-15]
```

The guarantees are strong: the balls are pairwise disjoint and each
contains exactly one root, so the output proves for example that
$x^5 - x - 1$ has exactly one real root, and the zero imaginary parts
are exact (the real roots are proved real, not merely close to the real
axis). Internally, the two functions work as follows.

- *Real roots.* `arb_fmpz_poly_real_roots` first isolates the real roots
  exactly, with Descartes' rule of signs on the integer polynomial
  (`fmpz_poly_isolate_real_roots`), which yields rational numbers that
  are roots and disjoint intervals with dyadic endpoints each containing
  one root. It then refines each interval in ball arithmetic with
  `arb_fmpz_poly_refine_root_arb`, which combines interval Newton steps
  (see [](tutorials/numerics/interval-newton)) with bisection. A
  polynomial in $x^d$ is handled by computing the roots of the deflated
  polynomial and taking $d$-th roots.
- *Complex roots.* `arb_fmpz_poly_complex_roots` isolates the real roots
  in the same way. If all roots are real, it is done. Otherwise the
  nonreal roots are computed numerically with the Durand–Kerner
  (Weierstrass) iteration of `acb_poly_find_roots`, with the precision
  doubled from 53 bits until all roots are isolated and accurate enough;
  when there are many real roots or clusters of them, the real roots are
  first divided out, which keeps the iteration fast.

## Polynomials with inexact coefficients {#inexact}

`acb_poly_find_roots(roots, f, initial, maxiter, prec)` takes a
polynomial with ball coefficients. It runs the Durand–Kerner iteration
in floating point from initial values (given, or placed on circles
chosen from the Newton polygon of the coefficients [@Bin1996] if
`initial` is `NULL`) and then *validates* the approximations $z_k$ with
the following inclusion theorem, applied with ball arithmetic so that it
holds for every polynomial with coefficients in the coefficient balls.

::: proposition Root inclusion {#prop-inclusion}
Let $f$ be a polynomial of degree $n$ and $z \in \mathbb{C}$ with
$f'(z) \ne 0$. Then the closed disk with centre $z$ and radius
$n |f(z) / f'(z)|$ contains a root of $f$.
:::

::: proof
Let $\zeta_1, \ldots, \zeta_n$ be the roots. If $|z - \zeta_i| > r =
n |f(z)/f'(z)|$ for all $i$, then
$|f'(z)/f(z)| = |\sum_i 1/(z - \zeta_i)| < n / r = |f'(z)/f(z)|$, a
contradiction.
:::

If the $n$ disks obtained from $n$ approximations are pairwise disjoint,
each contains at least one root, and since there are only $n$ roots,
each contains exactly one. `acb_poly_find_roots` returns the number of
disks that are disjoint from all others. If it equals $n$, all roots are
isolated; otherwise the precision (or the number of iterations) was
insufficient, or the polynomial has multiple roots, or roots so close
that the working precision cannot separate them. The caller decides
what to do, typically double the precision and try again, recomputing
the coefficients at the new precision.

For a polynomial with real coefficients, the nonreal roots come in
conjugate pairs, but a numerical root very close to the real axis has a
disk that overlaps it. `acb_poly_validate_real_roots` decides which
disks contain real roots: if it returns nonzero, every disk that touches
the real axis contains a real root, and its imaginary part can be set to
zero.

```c example name=numerics_poly_roots_acb
#include <flint/acb.h>
#include <flint/acb_poly.h>

int main(void)
{
    acb_poly_t f;
    acb_ptr roots;
    acb_t t;
    slong i, prec, deg = 6, isolated, acc, goal = 60;

    acb_poly_init(f);
    acb_init(t);
    roots = _acb_vec_init(deg);

    for (prec = 32; ; prec *= 2)
    {
        /* f = x^6 + pi x^3 + e x - 1/3, rounded to prec bits */
        acb_poly_zero(f);
        acb_poly_set_coeff_si(f, 6, 1);
        arb_const_pi(acb_realref(t), prec);
        acb_poly_set_coeff_acb(f, 3, t);
        arb_const_e(acb_realref(t), prec);
        acb_poly_set_coeff_acb(f, 1, t);
        acb_set_si(t, -1);
        acb_div_ui(t, t, 3, prec);
        acb_poly_set_coeff_acb(f, 0, t);

        isolated = acb_poly_find_roots(roots, f, NULL, 0, prec);
        acc = WORD_MAX;
        for (i = 0; i < deg; i++)
            acc = FLINT_MIN(acc, acb_rel_accuracy_bits(roots + i));
        flint_printf("prec = %wd: %wd isolated roots, accuracy %wd bits\n",
            prec, isolated, acc);

        if (isolated == deg && acc >= goal &&
                acb_poly_validate_real_roots(roots, f, prec))
            break;
    }

    for (i = 0; i < deg; i++)
    {
        if (arb_contains_zero(acb_imagref(roots + i)))
            arb_zero(acb_imagref(roots + i));       /* proved real */
    }

    _acb_vec_sort_pretty(roots, deg);
    for (i = 0; i < deg; i++)
    {
        acb_printn(roots + i, 18, 0);
        flint_printf("\n");
    }

    _acb_vec_clear(roots, deg);
    acb_poly_clear(f);
    acb_clear(t);
    flint_cleanup();
    return 0;
}
```

```output
prec = 32: 6 isolated roots, accuracy 25 bits
prec = 64: 6 isolated roots, accuracy 57 bits
prec = 128: 6 isolated roots, accuracy 121 bits
[-1.62045603282535763 +/- 2.78e-18]
[0.120598234331127720 +/- 9.67e-22]
[-0.139021211652347645 +/- 3.33e-19] + [0.881122870929321444 +/- 4.89e-19]*I
[-0.139021211652347645 +/- 3.33e-19] + [-0.881122870929321444 +/- 4.89e-19]*I
[0.888950110899462602 +/- 2.82e-19] + [1.16335399922150313 +/- 2.00e-19]*I
[0.888950110899462602 +/- 2.82e-19] + [-1.16335399922150313 +/- 2.00e-19]*I
```

The roots are isolated already at 32 bits, but only to 25 bits; the
accuracy goal of 60 bits requires 128. `_acb_vec_sort_pretty` sorts the
roots by their real and imaginary parts for display; it is not meant as
a rigorous ordering.

Durand–Kerner converges for simple roots, but slowly in the presence of
clusters and not at all to multiple roots; the polynomial should be made
squarefree first when possible. Since the iteration runs on all roots
simultaneously with $O(n^2)$ operations per step, the cost grows at
least quadratically with the degree. Well-separated roots are found
quickly even for fairly large degrees (`poly_roots a 400`, degree 400
with all roots nonreal, takes less than 0.1 seconds), while clustered
roots need more iterations and a higher precision.

## The example program {#programs}

`examples/poly_roots.c` computes the roots of an integer polynomial with
`arb_fmpz_poly_complex_roots`, after a squarefree factorisation, and
prints the progress (the flag `ARB_FMPZ_POLY_ROOTS_VERBOSE`):

```
poly_roots [-refine d] [-print d] [-real] <poly>
```

`-refine d` computes the roots to a relative accuracy of $10^{-d}$ (by
default only to the accuracy needed for isolation), `-print d` prints
them with $d$ digits, and `-real` computes only the real roots with
`arb_fmpz_poly_real_roots`. The polynomial is a product of factors given
as `a n` ($1 + 2x + \cdots + (n+1)x^n$), `t n` and `u n` (Chebyshev
$T_n$, $U_n$), `p n` (Legendre $P_n$), `c n` (cyclotomic $\Phi_n$),
`s n` (Swinnerton-Dyer $S_n$), `b n` (Bernoulli $B_n$), `w n`
(Wilkinson, $\prod_{k=1}^n (x - k)$), `e n` (the Taylor polynomial of
$e^x$ of degree $n$), `m n m` ($x^n + (100x + 1)^m$) and
`coeffs c0 c1 ... cn`; for example `p 5 t 6 coeffs 1 2 3` is
$P_5(x) T_6(x) (1 + 2x + 3x^2)$.

```
> build/examples/poly_roots -refine 100 -print 20 b 16
computing squarefree factorization...
cpu/wall(s): 0 0
searching for 16 roots, 16 deflated
prec=53: 16 isolated roots (8 real) | cpu/wall(s): 0.001 0.003
prec=128: 16 isolated roots (8 real) | cpu/wall(s): 0.001 0
prec=256: 16 isolated roots (8 real) | cpu/wall(s): 0 0
prec=512: 16 isolated roots (8 real) | cpu/wall(s): 0.001 0.001
done!
16 roots with multiplicity 1
[-0.94308706466055783383 +/- 2.02e-21]
[-0.75534059252067985752 +/- 2.70e-21]
[-0.24999757119077421009 +/- 4.27e-21]
[0.24999757152512726002 +/- 4.43e-21]
[0.75000242847487273998 +/- 4.43e-21]
[1.2499975711907742101 +/- 1.43e-20]
[1.7553405925206798575 +/- 1.74e-20]
[1.9430870646605578338 +/- 3.21e-20]
[-0.99509334829256233279 +/- 9.41e-22] + [0.44547958157103608805 +/- 3.59e-21]*I
[-0.99509334829256233279 +/- 9.41e-22] + [-0.44547958157103608805 +/- 3.59e-21]*I
[1.9950933482925623328 +/- 1.10e-20] + [0.44547958157103608805 +/- 3.59e-21]*I
[1.9950933482925623328 +/- 1.10e-20] + [-0.44547958157103608805 +/- 3.59e-21]*I
[-0.92177327714429290564 +/- 4.68e-21] + [1.0954360955079385542 +/- 1.71e-21]*I
[-0.92177327714429290564 +/- 4.68e-21] + [-1.0954360955079385542 +/- 1.71e-21]*I
[1.9217732771442929056 +/- 3.54e-20] + [1.0954360955079385542 +/- 1.71e-21]*I
[1.9217732771442929056 +/- 3.54e-20] + [-1.0954360955079385542 +/- 1.71e-21]*I
cpu/wall(s): 0.004 0.005
```

Roots of different multiplicities are separated by the squarefree
factorisation; for example
`poly_roots -print 5 p 5 p 5 t 7 coeffs 1 5 10 10 5 1`, that is $P_5^2 T_7 (1 + x)^5$, reports the six roots of
$T_7 / x$ with multiplicity 1, the four nonzero roots of $P_5$ with
multiplicity 2, the root 0 with multiplicity 3 (from $P_5^2 T_7$) and
the root $-1$ with multiplicity 5. The Wilkinson polynomial
$\prod_{k=1}^{100}(x - k)$, notorious for the ill-conditioning of its
roots in floating point, is solved with `poly_roots -print 15 w 100` in a
few milliseconds, since its roots are all real and the exact real root
isolation does not suffer from the conditioning.
