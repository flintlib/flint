---
title: Correctly rounded elementary functions
summary: How arf rounds the elementary functions correctly - certification of a fixed-point enclosure, Ziv's strategy of retrying with more guard bits, why it terminates, the exact and asymptotic cases that must be recognised separately - and the relaxed and fixed-precision variants of arf, nfloat and acf.
complexity: That of one evaluation at the target precision plus 20 to 24 guard bits; a retry, with probability about $2^{-20}$, doubles the guard bits.
order: 40
---

The elementary functions of `arf` (`arf_exp`, `arf_log`, `arf_sin`,
`arf_atan`, `arf_pow` and the others of [](arf#elementary)) return the
function value rounded to `prec` bits in a given direction. In the
strict modes (`ARF_RND_DOWN`, `ARF_RND_UP`, `ARF_RND_FLOOR`,
`ARF_RND_CEIL` and `ARF_RND_NEAR`) the result is *correctly rounded*: it
is the representable number that rounding the exact value would give,
as for the arithmetic operations of IEEE 754 [@IEEE754]. This page
describes how that is achieved on top of the fixed-point kernels of
[](algorithms/special-functions/elementary), and the variants with a
weaker guarantee: the relaxed modes of `arf`, the fixed-precision
functions of `nfloat` and the complex functions of `acf`.

A kernel returns a fixed-point enclosure, a *frame* $Y 2^e$ with error
counts $\varepsilon_-$, $\varepsilon_+$ such that the value $V$ lies in
$[(Y - \varepsilon_-) 2^e, (Y + \varepsilon_+) 2^e]$ (one of the counts
is zero for a one-sided bound). Correct rounding needs more than a
small error: it needs to know on which side of every rounding boundary
$V$ lies, which no fixed precision guarantees, since $V$ may lie
arbitrarily close to a boundary. This is the table maker's dilemma
[@MBD2018].

## Certifying a frame {#certify}

Let $k$ be the position (in units of $2^e$) of the last bit kept at
`prec` bits, so that $Y = h 2^{k-1} + D$ with the bits of $Y$ above
$2^{k-1}$ in $h$ and the discarded part $0 \le D < 2^{k-1}$, and let
$2^{eb}$ bound the error counts. The rounding of $V$ is decided when the
enclosure lies strictly between two consecutive multiples of
$2^{k-1}$, that is between two rounding boundaries (for the directed
modes the boundaries are the representable numbers $\pm$ multiples of
$2^k$, for rounding to nearest also the midpoints between them):

- **two-sided bounds:** the bits $[eb, k - 1)$ of $Y$ are neither all
  zero nor all one;
- **one-sided bounds:** when $V$ is known to lie above $Y$, all zero is
  fine as well (and so is a $Y$ exact at `prec` bits); when $V$ lies
  below $Y$, all one is fine;
- **rounding to nearest:** when the window is all equal to the bit at
  $2^{k-1}$ and $eb \le k - 4$, $V$ lies within a quarter ulp of a
  representable number, far from the ties, and rounds to it from either
  side. A value exponentially close to a representable number (such as
  $\tanh x = 1 - 2e^{-2x} + \cdots$ for large $x$) therefore needs no
  more precision than any other in this mode.

In each decided case $V$ rounds like $Y$ with its bits below the half
bit replaced by a single sticky bit, which an exact rounding of the
integer $Y$ then performs. The test reads a window of bits and costs
nothing next to the evaluation.

## Ziv's strategy {#ziv}

If the window does not decide the rounding, the evaluation is repeated
with more guard bits [@Ziv1991]:

::: algorithm Rounding by retries {#alg-ziv}
**Input:** $x$, `prec`, a strict rounding mode.

**Output:** $f(x)$ correctly rounded.

1. Set $\mathit{extra} = 0$.
2. Evaluate the kernel at $\mathit{prec} + 24 + \mathit{extra}$ bits and
   certify the frame. If the rounding is decided, return.
3. Set $\mathit{extra} = 32$ if it was 0, else double it, and go to 2.
:::

With $24$ guard bits the first attempt fails with probability about
$2^{-20}$ for a value whose bits are random, and each retry squares the
chance of failure, so that the expected cost is that of one evaluation.
The doubling makes the total cost at most a constant factor more than
the last attempt.

## Termination {#termination}

The loop terminates whenever $V$ is not a dyadic number, since then a
frame accurate enough excludes every rounding boundary. For the
elementary functions at dyadic arguments this is a consequence of the
Lindemann–Weierstrass theorem: $e^\alpha$ is transcendental for every
algebraic $\alpha \ne 0$, and so are $\log \alpha$ for algebraic
$\alpha \ne 0, 1$, $\sin \alpha$, $\cos \alpha$, $\tan \alpha$,
$\operatorname{atan} \alpha$ and the hyperbolic functions and other
inverse functions at algebraic $\alpha \ne 0$. The remaining cases have
algebraic values and are recognised before any evaluation:

- the values at $0$ (and $\cos 0 = 1$, $\operatorname{acos} 1 = 0$, and
  so on), and the logarithms $\log 1 = 0$, $\log_2 2^k = k$,
  $\log_{10} 10^k = k$ and the powers $2^k$, $10^k$ of `arf_exp2` and
  `arf_exp10` at integers;
- the functions of $\pi x$, whose values at rational $x$ are algebraic:
  by Niven's theorem the only rational values of $\sin(\pi x)$ and
  $\cos(\pi x)$ are $0$, $\pm 1/2$ and $\pm 1$, which at dyadic $x$ means
  the multiples of $1/2$; the only rational values of $\tan(\pi x)$ are
  $0$ and $\pm 1$, at the multiples of $1/4$; and correspondingly
  $\operatorname{atan}(x)/\pi$, $\operatorname{asin}(x)/\pi$ and
  $\operatorname{acos}(x)/\pi$ are dyadic only at $x \in \{0, \pm 1\}$;
- powers $x^y$, which are algebraic for dyadic $x$ and $y$: `arf_pow`
  computes exactly every power whose result is dyadic (integer powers,
  and the $2^f$-th roots of perfect powers for $y = n/2^f$) before
  evaluating, so that the loop only sees irrational values.

Other arguments give irrational values (such as $2^{1/2}$ from
`arf_exp2`), for which the loop terminates as well. No bound on the
number of retries is known at multiple precision (exhaustive searches
have bounded the precision needed for the hardest cases of the
double-precision formats [@MBD2018]). In practice a retry is rare, and
the asymptotic cases below, which are the ones where an argument of
special form puts the value very close to a representable number, are
decided without retries.

## Asymptotic cases {#asymptotic}

Near some points the value is a dyadic number $c$ times $1 + \delta$
with a tiny $\delta$ of known sign: $e^x = 1 + x + \cdots$ and
$\cos x = 1 - x^2/2 + \cdots$ for tiny $x$, $\sin x = x(1 - x^2/6 +
\cdots)$, $\operatorname{atan} x$, $\tanh x$ and $\sinh x$ likewise,
$\log(1 + x) = x(1 - x/2 + \cdots)$, $\cot x = (1/x)(1 - x^2/3 -
\cdots)$, and, for huge arguments, $\tanh x = 1 - 2e^{-2x} + \cdots$,
$\operatorname{atan} x = \pi/2 - 1/x + \cdots$ and so on. When
$|\delta|$ lies below the last bit, a two-sided bound would never decide
a directed rounding of such a value (its window is all zero or all
one), whatever the precision. The functions recognise these arguments
from their exponents and round $c(1 + \delta)$ directly from $c$ and the
sign of $\delta$ with a one-sided bound, without evaluating the
function, at any exponent; near the singular points of the inverse
functions (such as $\operatorname{acos} x$ near $x = 1$, which is
$\sqrt{2(1 - x)}(1 + \delta)$) the same is done with the exact square
root when it is dyadic. A quotient $(1/x)(1 + \delta)$ is decided by one
division with remainder.

## Relaxed modes {#relaxed}

The modes `ARF_RND_FAST` and `ARF_RND_ACCURATE` drop correct rounding
for speed: the first frame is final. `ARF_RND_ACCURATE` rounds the frame
to nearest, which with 20 guard bits gives an error below $0.51$ ulp,
and `ARF_RND_FAST` truncates the frame plus its upper error count,
within one ulp. In these modes the inexact flag is set whenever a kernel
was evaluated; the exact cases above are still detected before any
evaluation and reported as exact, except the exact powers of `arf_pow`
(and `arf_exp10`), which would cost a test on every call: the register
path evaluates them within its error bound like any other power.

## Compositions and quotients {#compositions}

The less common functions are compositions of the kernels: the inverse
trigonometric and hyperbolic functions of $\operatorname{atan2}$,
$\operatorname{log1p}$ and square roots, $\log_{10}$ of $\log$, the
$\operatorname{sinc}$ functions of the sine, and so on. They evaluate
the inner functions in a relaxed mode at a few more bits on temporaries,
bound the error of the composition explicitly (the relative errors of
the factors add up, with the condition numbers of the outer operations),
and certify the result like a frame, retrying at a higher working
precision until the rounding is decided. The reciprocal functions
($\cot$, $\sec$, $\csc$ and their hyperbolic counterparts) and the
tangent from a sine and cosine divide a frame by another frame or by an
exact mantissa at the precision of the current attempt, with the error
of the quotient bounded from those of the operands, and certify the
quotient.

## Large exponents {#exponents}

The arguments and results of `arf` may have exponents beyond a word. The
exponential of a large argument and the powers return their binary
exponent as an exact integer; the frame is rounded at exponent zero and
shifted afterwards, so that the rounding does not depend on the size of
the exponent. The cost of the argument reduction, however, grows with
the exponent of the argument (see
[](algorithms/special-functions/elementary#reduction)), so `arf` uses
the same magnitude limits as `arb`: the exponential and the hyperbolic
functions of $|x| > 2^{\max(128, 2\,\mathit{prec})}$ and the
trigonometric functions of $|x| > 2^{\max(65536, 4\,\mathit{prec})}$
give NaN rather than a result whose cost would not be polynomial in
`prec`. The functions of $\pi x$ reduce exactly and have no limit.

## Fixed precision: nfloat {#nfloat}

The elementary functions of `nfloat` evaluate a kernel once, at the
function precision of the context plus 12 guard bits (16 from 4096
bits on), and round the
frame to the context's precision: there is no retry, and the result is
within about one ulp at the function precision. With the directed
rounding flags `NFLOAT_RND_FLOOR` and `NFLOAT_RND_CEIL` the error count
is added in the direction of rounding, so that the result is a valid
lower or upper bound within a few ulps, and the asymptotic cases above
give the bound $c(1 \pm 2^{-q})$ in the right direction. The exact cases
give exact results, as in `arf`.

## Complex numbers: acf {#acf}

The complex functions of `acf` round each part of the result. In the
strict modes each part is correctly rounded, whatever its size against
the other part; in the relaxed modes each part is within 1 resp. 0.51
ulp of the larger part, which spares the cost of resolving a part that
is exponentially smaller than the other. The parts are compositions of
the real functions and their certification follows the same scheme; see
[](algorithms/special-functions/complex-elementary).

## Use in FLINT {#use}

- The elementary functions of [](arf#elementary), [](nfloat#elementary)
  and [](acf#elem).
- `arf_const_pi` certifies the floor of $\pi/4$ of `mp_real` at
  increasing precision in the same way.

::: functions arf_exp arf_log arf_sin_cos arf_atan arf_atan2 arf_pow arf_const_pi :::
