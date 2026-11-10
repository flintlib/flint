---
title: Elementary functions
summary: How FLINT evaluates exp, log, the trigonometric and hyperbolic functions and their inverses - the fixed-point kernels of mp_real (argument reduction, lookup tables, series, bit-burst evaluation, multi-prime reduction, Newton–Taylor inversion, the AGM) with their error bounds, and how arf, nfloat, arb, acb and acf wrap them.
complexity: $O(n^2)$ word operations up to a few hundred limbs; $O(M(n)\log^2 n)$ at high precision ($O(M(n)\log n)$ for the AGM logarithm), plus the cost of the argument reduction for huge arguments.
order: 10
---

All elementary functions in FLINT -- the exponential, the logarithm,
sine, cosine and tangent, the arctangent, the functions of $\pi x$, and
the hyperbolic, inverse and reciprocal functions derived from them --
are evaluated by one set of fixed-point *kernels* in the `mp_real`
module. The kernels cover precisions from one limb, where an evaluation
takes a few dozen nanoseconds, to hundreds of millions of bits. The
number types wrap them: `arf` rounds their results correctly, `nfloat`
rounds them once to a fixed precision, `arb` turns them into balls,
and `acb`, `acf` and `nfloat_complex` combine the real functions into
complex ones. The design follows [@Joh2014c] (medium precision) and
[@HJ2024] (high precision with precomputation), with the classical
algorithms of Brent [@Bre1976agm] [@Bre1976b]; [@Mul2016] and
[@BZ2011] are general references.

Throughout, $B = 2^{\bits}$, $n$ is a number of limbs and $p$ a number
of bits.

## Overview {#overview}

An evaluation has three stages.

1. **Argument reduction.** The argument $x$, a floating-point number,
   is reduced to a fixed-point number $v \in [0, 1)$: by multiples of
   $\log 2$ for the exponential, by the exponent for the logarithm, by
   multiples of $\pi/2$ for the trigonometric functions (exactly by
   $1/2$ for the functions of $\pi x$), and to $|x| < 1$ for the
   arctangent ([](#reduction)).
2. **The kernel** computes $f(v)$ for $v \in [0, 1)$ at $n$ limbs with
   an error bound in units of $B^{-n}$. The kernel is chosen by size:

   | Limbs | exp, sin, cos | log | atan |
   |---|---|---|---|
   | 1–3 | [table-driven](algorithms/special-functions/elementary-tables#small) | table-driven | table-driven |
   | 4–32 | [dense tables](algorithms/special-functions/elementary-tables#dense) | dense tables | dense tables |
   | 33–600 | [bitwise](algorithms/special-functions/elementary-tables#bitwise) | bitwise | bitwise |
   | 601–65536 | [multi-prime](algorithms/special-functions/diophantine-reduction) | [Newton–Taylor](#newton) | Newton–Taylor |
   | above | [table-free](algorithms/special-functions/elementary-series#notab) | Newton–Taylor | Newton–Taylor |

   The table-based kernels reduce $v$ further, against tables of
   logarithms or arctangents, to a remainder $t < 2^{-r}$; the
   multi-prime kernels reduce by integer combinations of logarithms of
   primes (or arguments of Gaussian primes); the table-free kernels halve
   the argument.
3. **The reduced series.** The function at the remainder $t$ is a Taylor
   series, evaluated by [rectangular
   splitting](algorithms/special-functions/elementary-series#series) at
   moderate precision and by the [bit-burst
   algorithm](algorithms/special-functions/elementary-series#bit-burst)
   at high precision. The series is also called directly when the
   argument itself is small enough, since its cost falls with $r$ while
   that of a table-based reduction does not.

The reconstruction undoes the reductions: a product of the factors
$1 + 2^{-i}$ or a power of two for the exponential, a sum of tabulated
logarithms for the logarithm, a rotation or the half-angle formulas for
the sine and cosine. Each stage has an error bound in units of its last
limb, and the bounds add up to the error count that the kernel returns
([](#frames)).

## Argument reduction {#reduction}

Below, $m$ is the argument, an exact dyadic number (the midpoint, for a
ball).

**Exponential.** For $|m| \ge \log 2$, write $|m| = q\log 2 + t$ with
$t \in [0, \log 2)$. The quotient $q$ is computed from the top limbs of
$m$ and two limbs of $1/\log 2$ (stored as the fraction
$1/\log 2 - 1$), which gives $q$ or $q - 1$; one multiply-subtract by
the cached floor of $\log 2$ and a correction give $t$. For $m < 0$ the
quotient is rounded up and the remainder complemented, so that
$\exp(m) = 2^{-q}\exp(t')$ with $t' \in [0, \log 2)$ and no division.
In the ball functions of `mp_real` the bound $|m| < 2^{\bits - 5}$
keeps the exponent of the result in the safe range of
[](algorithms/floating-point/mp-real-balls); the floating-point entry
points return a larger quotient as a multiprecision integer.

Both reductions turn a short argument into a full-length one, which
costs a full exponential however few bits the argument has, whereas the
cost of the bit-burst evaluation grows with the number of bits. So
beyond the bitwise kernels, $m \in (-1, 0)$ is evaluated as
$1/\exp(|m|)$ when $|m| < 2^{-128}$, or when $|m|$ is short enough for
the next rule (always above 65536 limbs, where the table-free kernel is
a bit-burst anyway); and in the range of the diophantine kernel, an
argument $v \in [0, 1)$ with $z$ leading zero bits and $\ell$
significant bits below them goes to the bit-burst directly when $\ell$
is below a bound that grows with $z$ and the precision: from 64 bits
($z \ge 16$) at 600 to 2000 limbs to 1024 ($z < 16$), 4096 ($z \ge 16$)
and 32768 bits ($z \ge 64$) from 30000 limbs. At 51907 limbs, for
example, $\exp(-x)$ with $x$ a 64-bit number below $2^{-64}$ takes
0.046 s instead of 0.31 s, and $\exp(x)$ 0.039 s instead of 0.30 s (the
profile program `p-exp_window` measures these choices).

**Logarithm.** For $m = 2^E u$ with $u \in [1/2, 1)$,
$\log m = (E - 1)\log 2 + \operatorname{log1p}(2u - 1)$ (up to 600
limbs) or $E\log 2 - (-\log u)$ (beyond). Near $m = 1$ both forms cancel
up to $z + 1$ bits for $|m - 1| < 2^{-z}$, and the working precision is
raised by that much; for larger $z$ the function uses
$\log(1 \pm s) = \pm 2\operatorname{atanh}(s/(2 \pm s))$ with $s = |m - 1|$
and the reduced atanh series.

**Sine and cosine.** For $|m| \ge 1$, the quotient $q$ approximating
$2|m|/\pi$ from below makes $t = |m| - q\pi/2 \ge 0$; a fold
$v = \pi/2 - t$ for $t > \pi/4$ brings the argument into $[0, \pi/4]$,
and the outputs are permuted and negated according to $q \bmod 4$. The
constant $2/\pi$ for the quotient and $\pi/4$ for the remainder come
from the caches. The reduction is absolute: near a nonzero multiple of
$\pi$ the reduced argument is small and its leading bits cancel. The
floating-point entry points used by `arf` and `arb` repeat the
reduction with as many more bits as cancel, so that the result keeps
its relative accuracy; the ball functions of `mp_real` show the
cancellation in the radius instead.

**Tangent.** The tangent is reduced the same way and finished as
$\pm\tan v$ or $\mp 1/\tan v$. Near a zero or a pole, $v$ is small and
an absolute error is a relative error of the result, so the kernel runs
with the leading zero bits of $v$ as extra precision, repeating the
reduction at the higher precision when needed. This terminates because
$\pi$ has finite irrationality measure, at most $7.11$ [@ZZ2020]: a
$b$-bit dyadic $m$ satisfies $|m - k\pi/2| \gtrsim 2^{-7.2b}$ for large
$b$, so a repetition beyond $8b + 1024$ extra bits is never needed (and
is not made: the ball function of `mp_real` then reports failure).

**Functions of $\pi x$.** For $\sin(\pi m)$, $\cos(\pi m)$ and
$\tan(\pi m)$ the reduction is exact: the units bit and the fraction
bits of the dyadic midpoint give $|m| = a/2 + st$ with
$a \in \{0, 1, 2, 3\}$, $s = \pm 1$ and $t \in [0, 1/4]$, and
$v = \pi t$ goes to the kernel with the leading zero bits of $t$ as
extra precision. The results therefore keep their relative accuracy
near every zero and pole, and multiples of $1/2$ (and odd multiples of
$1/4$ for the tangent) give exact results.

**Inverse functions.** $\operatorname{atan}$ reduces $|m| \ge 1$ by
$\operatorname{atan}(m) = \pi/2 - \operatorname{atan}(1/m)$ (or $\pi/4 +
\operatorname{atan}((m - 1)/(m + 1))$ for $1 < m < 2$). The other
inverse and hyperbolic functions are compositions without cancellation,
listed under `mp_real_log1p_bits`; for example
$\operatorname{acosh}(m) = \operatorname{log1p}(d + \sqrt{d(m + 1)})$
with $d = m - 1$ computed exactly.

**Hyperbolic functions and powers.** $\sinh$, $\cosh$ and $\tanh$ are
formed from $E = e^{|x|}$ and one approximate reciprocal, with the
precision raised by the bits that cancel in $E - 1/E$ for $|x| < 1$.
$x^y$ is $\exp(t)$ with $t = y \log|x|$: the logarithm is computed to as
many more bits as the magnitude of $t$ needs, multiplied exactly by the
mantissa of $y$, and the exponential of the result is returned with the
binary exponent of the reduction.

**Huge arguments.** The reductions are exact or carried to the precision
the exponent of the argument requires, so their cost grows with the
exponent: linearly for the exponential, whose quotient by $\log 2$ has
as many bits as the exponent, and for the trigonometric functions, which
need $\pi$ to about $e + p$ bits for an argument $2^e$. The wrappers
therefore bound the exponent in terms of the precision (see [](#arb)).

## The kernels {#kernels}

The algorithms of the kernels have pages of their own:

::: summary algorithms/special-functions/elementary-tables :::

::: summary algorithms/special-functions/elementary-series :::

::: summary algorithms/special-functions/diophantine-reduction :::

::: summary algorithms/special-functions/rational-sines :::

The logarithm and the arctangent at more than 600 limbs, and the AGM
logarithm, are described here.

## Newton–Taylor inversion {#newton}

Above 600 limbs the logarithm and the arctangent are computed from the
forward functions. If $t$ is any real number,
$$
-\log x = t - \log(1 + w), \quad w = x e^t - 1, \qquad
\operatorname{atan} x = t + \operatorname{atan} w, \quad
w = \frac{x\cos t - \sin t}{x\sin t + \cos t},
$$
exactly (the second for $|{\operatorname{atan} x - t}| < \pi/2$). If $t$
is correct to $p_0$ bits, then $|w| < 2^{-p_0}$, and $p/p_0$ terms of the
Taylor series of $\log(1 + w)$ or $\operatorname{atan} w$ give $p$ bits;
since the $k$-th term has $kp_0$ leading zero bits, the series costs a
fraction of one multiplication. The starting value is computed at
$p_0 \approx p/(N + 1)$ bits (for the arctangent, whose series is in
$w^2$, $p/(2N)$) by the bitwise kernel or recursively, so the whole
costs one forward evaluation plus about a tenth of one. Everything
around the forward function is ball arithmetic: the cancellation in $w$,
the forward function's error, the series and its tail are carried by
the radii. This is a single Newton step of high order in the sense of
[@BZ2011].

## The AGM logarithm {#agm}

The logarithm also has an $O(M(n)\log n)$ algorithm, from the
arithmetic-geometric mean [@Sal1976] [@Bre1976agm]. `_mp_real_neglog_agm`
uses the formula of Sasaki and Kanada [@SK1982],
$$
\log\frac{1}{q} = \frac{\pi}{\operatorname{agm}(\theta_2(q)^2,\, \theta_3(q)^2)},
\qquad 0 < q < 1,
$$
with $\theta_2(q) = 2q^{1/4}\sum_{k \ge 0} q^{k(k+1)}$ and $\theta_3(q) =
1 + 2\sum_{k \ge 1} q^{k^2}$. For $x \in [1/2, 1)$ it sets $r = x 2^{-e}$
and $q = r^4$, so that $q^{1/4} = r$ needs no root:
$$
-\log x = \frac{\pi/4}{\operatorname{agm}(4r^2 S_2^2,\, \theta_3^2)} - e\log 2,
\qquad S_2 = \sum_{k \ge 0} q^{k(k+1)}.
$$
The shift $e$ makes $q$ small enough that $N$ terms of the theta series
suffice ($e = \lceil p/(4(N + 1)^2)\rceil$); the merged exponents of the
two series are $\lfloor j^2/4 \rfloor$, consecutive ones differing by
$\lfloor j/2 \rfloor$, so the series cost about one product per term. The
first AGM step needs no square root, since $\theta_2^2\theta_3^2 =
(2rS_2\theta_3)^2$. The AGM itself runs about $\log_2 e$ iterations in
its slow phase and is finished by the series of `mp_real_agm`.

With warm caches the Newton–Taylor logarithm over the diophantine
exponential is faster up to about $7 \cdot 10^5$ limbs; without tables,
the AGM logarithm wins from about 20000 limbs, which is why
`mp_real_exp_agm` uses it.

## Fixed-point frames and error bounds {#frames}

A kernel returns a fixed-point *frame*: an array of limbs $Y$ and an
error count $\varepsilon$ such that the true value lies within
$\varepsilon$ units of $Y 2^{e}$, where $2^e$ is the unit of the frame
(a power of $B$, shifted by the binary exponent of the reduction). Some
bounds are one-sided -- the value lies in $[Y, Y + \varepsilon]$ or
$[Y - \varepsilon, Y]$ -- which the wrappers use to decide roundings
near exact values. The floating-point entry points of `mp_real`
(`_mp_real_exp_fp`, `_mp_real_log_fp`, `_mp_real_log1p_fp`,
`_mp_real_sin_cos_fp`, `_mp_real_tan_fp`, `_mp_real_atan_fp`,
`_mp_real_atan2_fp`, `_mp_real_sinh_cosh_fp`, `_mp_real_tanh_fp`,
`_mp_real_pow_fp` and the functions of $\pi x$; see [](mp_real#fp))
take the mantissa and exponent of a floating-point argument and a target
precision in bits, perform the reduction and choose the number of limbs
of the kernel, and return the frame. Their contract makes the error
bound independent of the path taken: the reduction contributes the
errors of the tabulated constants times the quotient (for example
$|q| (\log 2 - L)$ for the floor $L$ of $\log 2$, which the constant's
extra limbs keep below a unit), the kernel its own count, and the
reconstruction one unit per truncated product. The counts are bounded in
terms of the kernel's parameters (they grow linearly with the reduction
depth of the bitwise kernels, to a few thousand units at most) and are a
few units on the common paths; each wrapper turns the count into its
own error bound or rounding decision, so that a large count costs
accuracy (or a retry) but never correctness. Three situations need
special care:

- **Relative accuracy near zeros.** Near a zero of the function
  ($\sin x$ near a multiple of $\pi$, $\log x$ near 1, $\tan$ and
  $\operatorname{atan}$ near 0) an absolute error is a relative error of
  the result, so the entry points add the leading zero bits of the
  reduced argument to the number of limbs, and repeat the reduction at
  the higher precision when needed. The functions of $\pi x$, whose
  reduction is exact, keep their relative accuracy near every zero and
  pole.
- **Reciprocals.** For the cotangent, cosecant and the reciprocal
  hyperbolic functions the wrappers divide by the frame; the entry points
  choose their limbs so that the frame normalises with a short shift,
  and the division is done at the frame's own precision.
- **Tiny arguments.** For $|x| < 2^{-p/2}$ or so, $e^x$, $\log(1 + x)$,
  $\sin x$, $\operatorname{atan} x$ and so on are $1$ or $x$ times
  $1 + \delta$ with a tiny $\delta$ of known sign; the wrappers return
  that directly, at any exponent, without evaluating a kernel (and a
  Taylor polynomial of degree two or three serves a little further).

The thresholds between kernels are tuned by measurement and do not
affect the error bounds.

## Register fast paths {#fast}

On 64-bit machines, arguments with a mantissa of one or two limbs at a
target precision of at most one or two limbs are evaluated entirely in
registers by inline functions shared by `arf`, `nfloat`, `arb`, `acb`
and `acf` (`arf/elem.h`, `acf/elem.h`): they perform the reduction in
`double` arithmetic with exact integer corrections, evaluate the
table-driven kernels of three or four words, and return a frame with
an error count of a few units, falling back to the general entry points
for arguments outside their range. These paths give the elementary
functions their speed at 64 and 128 bits, where the general code would
spend most of its time in setup.

## From kernels to number types {#wrappers}

**arf** ([](algorithms/floating-point/correct-rounding)) rounds a frame
correctly when the bits between the rounding point and the error bound
decide the rounding, and otherwise repeats the evaluation with more
guard bits (Ziv's strategy); exact results and values of the form
$c(1 + \delta)$ with a tiny $\delta$ are handled without evaluation. The
relaxed modes `ARF_RND_FAST` and `ARF_RND_ACCURATE` round the first
frame.

**nfloat** evaluates once, at the function precision of the context plus
12 guard bits (16 from 4096 bits on, where the error bounds of the
bitwise kernels grow), and rounds the frame to the
precision of the context: the result is within about one ulp, and with
the directed rounding flags it is a valid bound, the error count being
added in the direction of rounding.

**arb** ([](#arb)) converts a frame for the exact midpoint into a ball
and adds the propagated error of the radius.

**acb, acf and nfloat_complex** combine the real functions by formulas
that keep each part accurate; see
[](algorithms/special-functions/complex-elementary).

## Balls in arb {#arb}

A point evaluation such as `arb_exp_arf`, `arb_log_arf` or
`arb_atan_arf` evaluates the entry point at the precision plus the same
guard bits as nfloat and converts the frame into a ball: the midpoint is $Y 2^e$
rounded to `prec` bits and the radius is the error count times $2^e$
plus the rounding error of the midpoint. The radius therefore stays
within a few ulps, and the result has close to `prec` bits of relative
accuracy, also near zeros. A ball function such as `arb_exp` receives
$[m \pm r]$ and proceeds as follows.

1. Special values (zero, infinities, NaN, an infinite radius) are
   handled directly, and so are arguments beyond the magnitude limits,
   which keep the cost polynomial in the precision
   ([](guide/types/ball-arithmetic#polynomial-time)): `arb_exp` returns
   $[0 \pm \infty]$ or a tiny ball when $|m| > 2^{\max(128, 2p)}$, and
   `arb_sin_cos` returns $[0 \pm 1]$ when $|m| > 2^{\max(65536, 4p)}$.
2. If $r > 0$, the precision is lowered to about the relative accuracy of
   the input plus 20 to 30 bits; more bits would be hidden by the
   propagated error. If the input has less than about 20 bits of
   accuracy, the function is evaluated on the interval instead, using
   its monotonicity (for $\sin$ and $\cos$, the endpoints and the extrema
   in between).
3. Otherwise the function of the exact midpoint is computed by the point
   evaluation, and the propagated error is added.

The propagated error uses a bound on the derivative over the ball
([](algorithms/ball-arithmetic/error-propagation#elementary)), which the
point evaluation provides nearly for free: the sine and cosine return
upper bounds of $|\sin m|$ and $|\cos m|$ read from the top bits of the
frames, so that the error of $\sin$ is bounded by $r |\cos m| + r^2/2$
and that of $\cos$ by $r |\sin m| + r^2/2$; the hyperbolic functions use
$r \cosh m \,(1 + 2r)$ for $r \le 1/2$ from the frame of $\cosh m$, and
the exponential $e^m (e^r - 1)$. For the functions of $\pi x$ the
bounds are multiplied by $\pi$. The logarithm bounds the error of
$\log [m \pm r]$, $r < m$, by $\log(1 + r/(m - r))$ and the arctangent by
$r / (1 + (|m| - r)^2)$.

## Derived functions {#derived}

- `arb_expm1` and `arb_log1p` use the entry points with the subtraction
  of 1 built in, at a precision increased by the bits that would be lost
  to cancellation for small arguments.
- `arb_sinh_cosh`, `arb_tanh` and `arb_coth` are computed from the
  hyperbolic entry points; `arb_asinh`, `arb_acosh` and `arb_atanh` from
  `arb_log`, `arb_log1p` and `arb_sqrt`.
- `arb_asin` and `arb_acos` use $\operatorname{atan}(x / \sqrt{1 -
  x^2})$, and `arb_atan2` combines `arb_atan` with $\pi$.
- `arb_sin_pi` and `arb_cos_pi` reduce $x$ modulo $1/2$ exactly before
  multiplying by $\pi$, so they are exact at integers and half-integers
  and accurate near the zeros; a dyadic $x$ with a small denominator,
  like a rational $x$ in `arb_sin_cos_pi_fmpq` and `arb_tan_pi_fmpq`, is
  evaluated as an algebraic number where that is faster
  ([](algorithms/special-functions/rational-sines)).
- `arb_pow` with exact arguments and a positive base (or an integer
  exponent) evaluates $x^y$ in one go; otherwise it is $\exp(y \log x)$
  in ball arithmetic, with the special cases of [](arb#powers).
- In `arf`, `nfloat` and `acf` the less common functions (`arf_asin`,
  `arf_acosh`, `arf_log10`, `arf_sinc_pi`, ...) are compositions of the
  ones above at a few more bits, with an explicit bound for the error of
  the composition; one-sided expansions near the singular points (such
  as $\operatorname{acos} x$ near $x = 1$, or $\operatorname{acoth} x$
  for huge $x$) avoid the stalls of two-sided bounds.
- Logarithms of smooth integers at high precision are assembled from
  cached logarithms of primes (`arb_log_ui_from_prev`,
  `arb_log_primes_vec_bsplit`).

The constants $\pi$ and $\log 2$ used by the reductions, and by
`arb_const_pi` and `arb_const_log2`, are kept as exact floors in static
tables and a shared cache ([](algorithms/special-functions/constants#cache)).

## Use in FLINT {#use}

- The ball functions of `mp_real` (`mp_real_exp_bits`,
  `mp_real_log_bits`, `mp_real_sin_cos_bits`, ...) apply the same
  reductions and kernels to `mp_real_t` balls, with the radius
  propagated like that of `arb` ([](mp_real/elementary)).
- The floating-point entry points serve `arf` ([](arf#elementary)),
  `nfloat` ([](nfloat#elementary)), `arb` ([](arb/elementary)), and
  through them `acb`, `acf` and `nfloat_complex`.

::: functions arb_exp arb_log arb_sin_cos arb_atan arf_exp arf_log arf_sin_cos arf_atan _mp_real_exp_fp _mp_real_log_fp _mp_real_sin_cos_fp _mp_real_atan_fp _mp_real_pow_fp _mp_real_neglog_newton _mp_real_atan_newton _mp_real_neglog_agm :::
