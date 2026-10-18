---
title: Precision and accuracy
summary: What the precision argument controls, how accuracy is measured, how to choose the precision with guard bits and doubling, and the polynomial time guarantee for huge inputs.
order: 3
---

Every ball function that can round takes a precision `prec` in bits. The
precision controls the *working precision*, not the accuracy of the
result: how many bits of the output are correct depends on the function,
on the input and on the conditioning of the problem. This page explains
the relation between the two and how to choose the precision.

## Measuring accuracy {#accuracy}

For a ball $[m \pm r]$ approximating a number $x$, define

- the absolute error $e_{\mathrm{abs}} = r$, an upper bound for $|x -
  m|$;
- the relative error $e_{\mathrm{rel}} = r / \max(0, |m| - r)$ (and 0
  if $r = m = 0$), an upper bound for $|x - m|/|x|$;
- the absolute and relative accuracy $1/e_{\mathrm{abs}}$ and
  $1/e_{\mathrm{rel}}$.

Measured in bits, one takes $\log_2$ of these quantities.
`arb_rel_accuracy_bits` returns a cheap estimate of the relative accuracy
in bits, computed from the exponents of $m$ and $r$, and `arb_get_mag` and
`arb_get_rad_arb` give the absolute error.

## Working precision and guard bits {#guard-bits}

When a fixed expression is evaluated, so that the sequence of operations
does not depend on the precision, the absolute or relative error of the
result is generically bounded by

$$
2^{O(1) - \mathit{prec}},
$$

where the $O(1)$ term depends on the expression and on the
implementation of the functions used. To get $p$ correct bits one
therefore needs a working precision of $p + O(1)$. For a numerically
well-behaved expression the $O(1)$ term is small, which leads to the
usual advice of adding a few *guard bits*: for most basic calculations,
10 or 20 extra bits suffice. If the $O(1)$ term is not known, increasing
the precision in exponential steps until the result is accurate enough is
a good strategy (it is due to Ziv for correctly rounded floating-point
functions). A partially accurate result can be used to estimate the
$O(1)$ term: if 120 bits of precision give 80 bits of accuracy, then a
little over 140 bits will probably give 100 bits.

The next program evaluates $(1 + 1/n)^n - e$ for $n = 10^{15}$, which
loses about 51 bits to cancellation, and doubles the precision until 53
bits are correct:

```c example name=ball_precision_loop
#include <flint/flint.h>
#include <flint/arb.h>

/* f = (1 + 1/n)^n - e for n = 10^15, which suffers from cancellation */
void f(arb_t res, slong prec)
{
    arb_t t, e;
    arb_init(t); arb_init(e);
    arb_set_ui(t, 1000000000000000);   /* n = 10^15 */
    arb_inv(t, t, prec);
    arb_log1p(t, t, prec);             /* log(1 + 1/n) */
    arb_mul_ui(t, t, 1000000000000000, prec);
    arb_exp(t, t, prec);
    arb_const_e(e, prec);
    arb_sub(res, t, e, prec);
    arb_clear(t); arb_clear(e);
}

int main(void)
{
    arb_t y;
    slong prec, goal = 53;

    arb_init(y);
    for (prec = goal + 10; ; prec *= 2)
    {
        f(y, prec);
        flint_printf("prec = %4wd: accuracy %3wd bits, f = ", prec,
                     arb_rel_accuracy_bits(y));
        arb_printn(y, 16, 0); flint_printf("\n");
        if (arb_rel_accuracy_bits(y) >= goal)
            break;
    }
    arb_clear(y);
    flint_cleanup();
    return 0;
}
```
```output
prec =   63: accuracy   8 bits, f = [-1.36e-15 +/- 2.24e-18]
prec =  126: accuracy  71 bits, f = [-1.359140914229521e-15 +/- 3.72e-31]
```

Note the use of `arb_log1p` to compute $\log(1 + 1/n)$: computing $1 +
1/n$ first and then its logarithm would lose another 50 bits.

## Simple and complex functions {#simple-complex}

The functions of FLINT lie on a spectrum between two extremes.

- *Simple operations*, such as the arithmetic operations and most
  elementary functions. For an input $x = [m \pm r]$, $f(x)$ is computed
  by evaluating $f(m)$ and bounding the *propagated error* $|f(m) - f(m +
  \varepsilon)|$, $|\varepsilon| \le r$, separately. The working precision
  is increased internally so that $f(m)$ has `prec` bits of relative
  accuracy up to a few units in the last place, and the propagated error
  is bounded tightly (see
  [](algorithms/ball-arithmetic/error-propagation)). The result is close
  to the best possible enclosure at the given precision, and the precision
  needed for a given accuracy is easy to estimate.
- *Complex operations*, such as many higher transcendental functions
  (for example the Riemann zeta function). They are evaluated by a
  sequence of simpler operations in ball arithmetic at a working precision
  of about `prec` bits, and the sequence may depend on `prec`; for
  example, an infinite series is truncated where the remainder is below
  $2^{-\mathit{prec}}$. The result can be far from the best enclosure,
  and convergence as $\mathit{prec} \to \infty$ is not guaranteed in
  general, although it holds in most cases in practice.

In short, the inclusion principle is the contract; how tight the
enclosures are is a matter of quality of implementation. Tighter bounds
let the user work at lower precision but may be more expensive to
compute, so each function makes a trade-off.

## The polynomial time guarantee {#polynomial-time}

FLINT provides a soft guarantee that the time used to evaluate a ball
function depends polynomially on `prec` and on the bit size of the input,
uniformly in the numerical value of the input. This allows ball
arithmetic to be used as a black box for numerical evaluation without the
risk that "bad" input, such as nested exponentials, makes a computation
hang or crash. By controlling the precision, the user can bound the
resources a computation uses; a result that is feasible but expensive can
still be forced by a high enough precision.

As motivation, consider $\sin x$ or $\exp x$ for the exact number $x =
2^{2^n}$, whose size in memory grows only linearly with $n$. Computing an
accurate approximation of $\sin x$ requires subtracting a multiple of
$2\pi$ that is accurate to $2^n$ bits, and $\exp x$ has an exponent of
$2^n$ bits; already $n = 33$ needs more than a gigabyte of memory, and
$n = 100$ is physically impossible. For functions computed directly from
their power series, without fast argument reduction, the cost is
exponential already for $x = 2^n$.

FLINT therefore caps the internal work parameters (the working precision,
the number of terms of a series, ...) by polynomial, usually linear,
functions of `prec`, and returns a crude bound when the cap is exceeded.
For example, `arb_sin` returns $[\pm 1]$ for arguments larger than
$2^{\max(65536, 4\,\mathit{prec})}$, and `arb_exp` returns $[\pm \infty]$
for large positive arguments and a tiny ball $[\pm 2^{-N}]$ for large
negative ones:

```c example name=ball_huge_input
#include <flint/flint.h>
#include <flint/arb.h>

int main(void)
{
    arb_t x, y;
    arb_init(x); arb_init(y);

    /* x = 2^(2^40): an exact number with a 41-bit exponent */
    arb_one(x);
    arb_mul_2exp_si(x, x, WORD(1) << 40);

    arb_sin(y, x, 64);
    flint_printf("sin(x)  = "); arb_printn(y, 10, 0); flint_printf("\n");
    arb_exp(y, x, 64);
    flint_printf("exp(x)  = "); arb_printn(y, 10, 0); flint_printf("\n");
    arb_neg(x, x);
    arb_exp(y, x, 64);
    flint_printf("exp(-x) = "); arb_printn(y, 10, 0); flint_printf("\n");

    arb_clear(x); arb_clear(y);
    flint_cleanup();
    return 0;
}
```
```output
sin(x)  = [+/- 1.00]
exp(x)  = [+/- inf]
exp(-x) = [+/- 1.69e-102435199438739363750012109250103232700]
```

The last bound is $2^{-2^{128}}$. This is not only a safeguard but
sometimes a useful optimisation. In an expression such as

$$
\log x + \sin(x) \exp(-x),
$$

which has the same structure as the reflection formula of the digamma
function, the second term is negligible for large $x$, and computing it
to high relative accuracy would be expensive and pointless; the crude
bound $[\pm 1] \cdot [\pm 2^{-N}]$ is enough to evaluate the whole
expression accurately.

The guarantee is soft: there are exceptions where the input is most
likely chosen deliberately to evaluate the function itself. For example,
the cost of the Riemann zeta function $\zeta(\sigma + it)$ grows linearly
with $|t|$ in the current implementation (the best known algorithms cost
$|t|^{\alpha}$ with $\alpha$ about $0.3$), and large $|t|$ are not cut
off.
