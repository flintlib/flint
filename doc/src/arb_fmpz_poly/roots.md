---
title: Isolating the roots of integer polynomials
functions: arb_fmpz_poly_real_roots arb_fmpz_poly_refine_root_arb
---

`arb_fmpz_poly_real_roots(roots, poly, flags, prec)` and
`arb_fmpz_poly_complex_roots(roots, poly, flags, prec)` compute all real,
respectively all complex, roots of a squarefree integer polynomial. The
results are balls that are pairwise disjoint and each contain exactly
one root, so the roots are isolated and counted correctly, and each has
a relative accuracy of at least `prec` bits. Unlike
`acb_poly_find_roots`, which may fail to isolate the roots at a given
precision, these functions always succeed: they increase the precision
internally until the roots are separated, so the running time depends
on the root separation of the polynomial, and `prec` only determines the
accuracy of the output.

::: summary algorithms/ball-arithmetic/polynomial-roots :::

## Real roots {#real}

The real roots are computed in four steps.

1. Deflation. A root at 0 is removed (the polynomial being squarefree,
   it is simple), and $f(x) = g(x^d)$ is written with the largest
   possible $d$ (`arb_fmpz_poly_deflation`). If $d$ is even, only the
   positive roots of $g$ are needed.
2. Exact isolation. `fmpz_poly_isolate_real_roots` (or
   `fmpz_poly_isolate_positive_roots`) computes the dyadic roots of $g$
   exactly and open intervals $(c \cdot 2^k, (c + 1) 2^k)$ that contain
   exactly one root each, by Descartes' rule of signs and the
   Vincent–Collins–Akritas bisection algorithm
   ([](algorithms/polynomials/real-roots-integers)).
3. Refinement. After dividing out the exact roots (so that the interval
   endpoints are not roots), each interval is refined to the target
   accuracy by `arb_fmpz_poly_refine_root_arb`, which combines interval
   Newton steps, ordinary Newton steps confirmed by sign changes, and
   bisection ([](algorithms/ball-arithmetic/real-root-isolation)).
4. Inflation. The roots of $f$ are the $d$-th roots of the positive
   roots of $g$ (and, for odd $d$, of the negative ones), plus their
   negatives for even $d$. If rounding makes two enclosures overlap
   (which is unlikely), the computation is repeated at a higher
   precision.

With the flag `ARB_FMPZ_POLY_ROOTS_RATIONAL`, the rational roots are
first computed exactly by `fmpz_poly_roots_fmpq` and divided out; they
are returned as exact balls when they are dyadic and to `prec` bits
otherwise. This is much faster when there are many rational roots, and
adds the cost of the rational root search otherwise.

## Complex roots {#complex}

`arb_fmpz_poly_complex_roots` first deflates the polynomial and isolates
its real roots as above. If all roots are real it is done. Otherwise it
runs `acb_poly_find_roots` on the deflated polynomial at 53, 128, 256,
... bits, reusing the previous approximations as starting values, until
all roots are isolated and accurate to `prec` bits, and takes $d$-th
roots of the results. Two variants are used:

- If the real roots are numerous (at least a quarter of all roots) or
  clustered (two of them closer than $1/(16n)$ times their magnitude),
  they are divided out in ball arithmetic, with enough precision to keep
  the quotient accurate, and only the nonreal roots are computed
  numerically. Clusters of real roots would otherwise slow down the
  Weierstrass–Durand–Kerner iteration considerably.
- Otherwise all roots are computed numerically and
  `acb_poly_validate_real_roots` proves which of the isolated roots are
  real.

The real roots are written first, in increasing order and with
imaginary part exactly zero, followed by the nonreal roots in conjugate
pairs, each pair starting with the root in the upper half-plane; the
second member of a pair is computed as the conjugate of the first.

## Examples {#examples}

The polynomial $f = (2x - 1)(3x + 1)(x^2 - 2)(x^2 + x + 1)$ has the
rational roots $1/2$ and $-1/3$, the real roots $\pm\sqrt{2}$ and a pair
of complex conjugate roots. With `ARB_FMPZ_POLY_ROOTS_RATIONAL`, the
dyadic root $1/2$ is returned as an exact ball. The output arrays must
have room for $\deg f$ entries even for the real roots.

```c example name=arb_fmpz_poly_roots_basic
#include <flint/fmpz_poly.h>
#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>

int main(void)
{
    fmpz_poly_t f;
    arb_ptr r;
    acb_ptr z;
    slong i, n;

    fmpz_poly_init(f);

    /* f = (2x - 1)(3x + 1)(x^2 - 2)(x^2 + x + 1) */
    fmpz_poly_set_str(f, "7  2 4 -9 -12 -8 5 6");
    r = _arb_vec_init(6);
    z = _acb_vec_init(6);

    n = arb_fmpz_poly_real_roots(r, f, ARB_FMPZ_POLY_ROOTS_RATIONAL, 53);
    flint_printf("%wd real roots:\n", n);
    for (i = 0; i < n; i++)
    {
        arb_printn(r + i, 15, 0);
        flint_printf("\n");
    }

    arb_fmpz_poly_complex_roots(z, f, 0, 53);
    flint_printf("complex roots:\n");
    for (i = 0; i < 6; i++)
    {
        acb_printn(z + i, 15, 0);
        flint_printf("\n");
    }

    _arb_vec_clear(r, 6);
    _acb_vec_clear(z, 6);
    fmpz_poly_clear(f);
    return 0;
}
```
```output
4 real roots:
[-1.41421356237310 +/- 4.96e-15]
[-0.333333333333333 +/- 3.34e-16]
0.500000000000000
[1.41421356237310 +/- 4.96e-15]
complex roots:
[-1.41421356237310 +/- 4.96e-15]
[-0.333333333333333 +/- 3.34e-16]
0.500000000000000
[1.41421356237310 +/- 4.96e-15]
[-0.500000000000000 +/- 2.01e-38] + [0.866025403784439 +/- 3.54e-16]*I
[-0.500000000000000 +/- 2.01e-38] + [-0.866025403784439 +/- 3.54e-16]*I
```

The roots of the Chebyshev polynomial $T_{100}$ are $\cos((2k+1)\pi /
200)$, $0 \le k < 100$; they are all real and cluster near $\pm 1$.
Since $T_{100}(x) = g(x^2)$, only the 50 positive roots of $g$ are
isolated and refined:

```c example name=arb_fmpz_poly_roots_chebyshev
#include <flint/fmpz_poly.h>
#include <flint/arb.h>
#include <flint/arb_fmpz_poly.h>

int main(void)
{
    fmpz_poly_t f;
    arb_ptr r;
    arb_t c;
    slong i, num, ok = 0, n = 100;

    fmpz_poly_init(f);
    arb_init(c);
    r = _arb_vec_init(n);

    /* the Chebyshev polynomial T_n has the roots cos((2k+1) pi / (2n)) */
    fmpz_poly_chebyshev_t(f, n);
    num = arb_fmpz_poly_real_roots(r, f, 0, 100);

    for (i = 0; i < num; i++)
    {
        /* the i-th smallest root is cos((2(n-1-i)+1) pi / (2n)) */
        arb_set_si(c, 2 * (n - 1 - i) + 1);
        arb_div_si(c, c, 2 * n, 128);
        arb_cos_pi(c, c, 128);
        ok += arb_overlaps(r + i, c);
    }

    flint_printf("%wd real roots, %wd as expected\n", num, ok);
    arb_printn(r + 0, 30, 0);
    flint_printf("\n");
    flint_printf("accuracy of r[0]: %wd bits\n", arb_rel_accuracy_bits(r + 0));

    _arb_vec_clear(r, n);
    arb_clear(c);
    fmpz_poly_clear(f);
    return 0;
}
```
```output
100 real roots, 100 as expected
[-0.999876632481660598638907127731 +/- 2.53e-31]
accuracy of r[0]: 141 bits
```

The accuracy is higher than the requested 100 bits because the
refinement converges quadratically and stops only once the target is
exceeded.
