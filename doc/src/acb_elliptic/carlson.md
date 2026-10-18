---
title: Carlson symmetric integrals
functions: acb_elliptic_rf acb_elliptic_rg acb_elliptic_rj acb_elliptic_rc1
---

Carlson's symmetric integrals are

$$
R_F(x,y,z) = \frac{1}{2} \int_0^{\infty} \frac{dt}{\sqrt{(t+x)(t+y)(t+z)}},
\qquad
R_J(x,y,z,p) = \frac{3}{2} \int_0^{\infty} \frac{dt}{(t+p)\sqrt{(t+x)(t+y)(t+z)}},
$$

$$
R_G(x,y,z) = \frac{1}{4} \int_0^{\infty}
\frac{t}{\sqrt{(t+x)(t+y)(t+z)}}
\left( \frac{x}{t+x} + \frac{y}{t+y} + \frac{z}{t+z}\right) dt,
$$

with the special cases $R_C(x, y) = R_F(x, y, y)$ and
$R_D(x, y, z) = R_J(x, y, z, z)$. The square root is the continuous
branch that is positive for large $t$ [@NIST2012, §19.16]. Every
elliptic integral, and in particular every Legendre integral of
[](acb_elliptic#legendre), can be written in terms of $R_F$, $R_D$, $R_J$
and $R_G$; the symmetric forms have simpler transformation laws and,
unlike the Legendre forms, a uniform algorithm. The definitions and
algorithms follow [@Car1995] and [@NIST2012, Ch. 19], with the
high-order Taylor expansion of [@Joh2019].

## Domains {#domains}

The integral for $R_F$ is well defined for $x, y, z$ in
$\mathbb{C} \setminus (-\infty, 0)$ with at most one of them zero; it is
extended to the negative real axis by the same formula, which then
depends on the side of the cut chosen by the principal square root. The
integral for $R_J$ also requires $p \notin (-\infty, 0]$. The functions
return indeterminate results when more than one of $x, y, z$, or when
$p$, contains zero, and when an input is not finite. $R_G$ has no
restrictions on its arguments: `acb_elliptic_rg` uses the formula
[@NIST2012, 19.21.10]

$$
2 R_G(x, y, z) = z R_F(x, y, z) - \tfrac13 (x - z)(y - z) R_D(x, y, z) + \sqrt{x y / z},
$$

after permuting the arguments so that $z$ does not contain zero (the
cases where two arguments are exactly zero give $\sqrt{w}/2$ directly).

## Method {#method}

`acb_elliptic_rf` and `acb_elliptic_rj_carlson` use Carlson's duplication
algorithm. The duplication theorem

$$
R_F(x, y, z) = R_F\!\left(\tfrac{x + \lambda}{4}, \tfrac{y + \lambda}{4}, \tfrac{z + \lambda}{4}\right),
\qquad \lambda = \sqrt{x}\sqrt{y} + \sqrt{y}\sqrt{z} + \sqrt{z}\sqrt{x},
$$

(and its analogue for $R_J$, which adds a sum of $R_C$ terms, computed
with `acb_elliptic_rc1`) brings the arguments closer together by a factor
about 4 per step. Once the relative deviations from their mean are small,
a multivariate Taylor series in the elementary symmetric polynomials of
the deviations is summed. Instead of Carlson's fixed series of degree 5
or 7 and many duplication steps, FLINT sums about $2.1 p^{0.4}$ terms
(real arguments) or $2.5 p^{0.4}$ terms (complex arguments) at precision
$p$, capped at 500, which balances the two costs at high precision, with
exact rational coefficients and a rigorous bound for the tail. The
working precision is lowered to the relative accuracy of the inputs plus
20 bits when the inputs are inexact. The details and proofs are in
[](algorithms/special-functions/elliptic-modular).

## Correctness of the duplication algorithm for $R_J$ {#rj}

For $R_F$ the duplication algorithm with principal square roots is valid
everywhere. For $R_J$ it is not: the square roots and the $R_C$ terms can
cross branch cuts, and then `acb_elliptic_rj_carlson` returns a wrong
value with a small radius. According to [@Car1995], a sufficient
condition for correctness is that $x, y, z$ have nonnegative real part
and $p$ has positive real part. `acb_elliptic_rj_carlson` does not check
anything.

`acb_elliptic_rj` uses the duplication algorithm when one of the
following sufficient conditions holds, and `acb_elliptic_rj_integration`
otherwise:

- one of $x, y, z$ is equal to $p$, that is, the function is $R_D$ (the
  duplication algorithm for $R_D$ is valid everywhere);
- $x, y, z$ have nonnegative real part and $p$ has positive real part;
- all four arguments are real;
- $\operatorname{Re}(p) \ge 0$ or $\operatorname{Im}(p) \ne 0$, and one
  of $x, y, z$ is real and nonnegative while the other two are complex
  conjugates.

Before integrating, it returns an indeterminate result if some argument
straddles the negative real axis. `acb_elliptic_rj_integration` chooses
$N$ with $x + N, y + N, z + N, p + N$ in the right half-plane and writes

$$
R_J(x, y, z, p) = \frac32 \int_0^N \frac{dt}{(t+p)\sqrt{(t+x)(t+y)(t+z)}} + R_J(x+N, y+N, z+N, p+N),
$$

computing the integral with `acb_calc_integrate` along a path deformed
into the upper or lower half-plane (and substituting $t \to t^2$ when an
argument is zero) and the second term by duplication. This is reliable
but much slower than the duplication algorithm alone.

## Example {#example}

The program computes $R_F(0, 1, 2) = \Gamma(1/4)^2 / (4 \sqrt{2\pi})$,
checks the Legendre relation
$E(m) = R_F(0, 1-m, 1) - \frac{m}{3} R_D(0, 1-m, 1)$ against
`acb_elliptic_e`, and evaluates an $R_J$ at complex arguments where the
duplication algorithm alone gives a wrong value.

```c example name=acb_elliptic_carlson
#include <flint/acb.h>
#include <flint/acb_elliptic.h>

int main(void)
{
    acb_t x, y, z, p, r, s, m;
    slong prec = 128;

    acb_init(x); acb_init(y); acb_init(z); acb_init(p);
    acb_init(r); acb_init(s); acb_init(m);

    acb_set_ui(y, 1);
    acb_set_ui(z, 2);
    acb_elliptic_rf(r, x, y, z, 0, prec);
    flint_printf("R_F(0,1,2)    = "); acb_printn(r, 30, 0); flint_printf("\n");

    /* E(m) = R_F(0, 1-m, 1) - (m/3) R_D(0, 1-m, 1) with m = 3/4 */
    acb_set_d(m, 0.75);
    acb_sub_ui(y, m, 1, prec);
    acb_neg(y, y);
    acb_one(z);
    acb_elliptic_rf(r, x, y, z, 0, prec);
    acb_elliptic_rj(s, x, y, z, z, 0, prec);
    acb_mul(s, s, m, prec);
    acb_div_ui(s, s, 3, prec);
    acb_sub(r, r, s, prec);
    flint_printf("via R_F, R_D  = "); acb_printn(r, 30, 0); flint_printf("\n");
    acb_elliptic_e(r, m, prec);
    flint_printf("E(3/4)        = "); acb_printn(r, 30, 0); flint_printf("\n");

    acb_set_d_d(x, 2.0, 3.0);
    acb_set_d_d(y, -1.0, -1.0);
    acb_set_d_d(z, -2.0, -1.0);
    acb_set_d_d(p, -4.0, 1.0);
    acb_elliptic_rj(r, x, y, z, p, 0, 64);
    flint_printf("R_J           = "); acb_printn(r, 10, 0); flint_printf("\n");
    acb_elliptic_rj_carlson(r, x, y, z, p, 0, 64);
    flint_printf("R_J (carlson) = "); acb_printn(r, 10, 0); flint_printf("\n");

    acb_clear(x); acb_clear(y); acb_clear(z); acb_clear(p);
    acb_clear(r); acb_clear(s); acb_clear(m);
    return 0;
}
```
```output
R_F(0,1,2)    = [1.31102877714605990523241979495 +/- 4.45e-30]
via R_F, R_D  = [1.21105602756845952480356289955 +/- 1.03e-30]
E(3/4)        = [1.21105602756845952480356289955 +/- 1.03e-30]
R_J           = [0.1948611350 +/- 1.75e-11] + [-1.020158039 +/- 3.81e-10]*I
R_J (carlson) = [-0.9832009839 +/- 1.42e-11] + [1.009854664 +/- 2.79e-10]*I
```

The value from `acb_elliptic_rj` agrees with a direct numerical
evaluation of the defining integral, $0.1948611\ldots - 1.0201580\ldots i$;
the duplication algorithm, applied outside its domain of validity,
returns a different number with a tight but wrong enclosure.
