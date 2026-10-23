---
title: Elementary functions in mp_real
summary: How mp_real evaluates exp, log, the trigonometric functions and their inverses from one limb to millions of limbs: argument reduction of balls, table-driven and dense-table kernels, bitwise (BKM-type) reduction with signed digits, rectangular splitting, bit-burst evaluation, Newton–Taylor inversion, the AGM logarithm, and sines of rational multiples of π.
complexity: $O(M(n))$ times a small constant for the table-based kernels at fixed $n$; $O(M(n) \log^2 n)$ for the bit-burst and diophantine kernels; $O(M(n) \log n)$ for the AGM logarithm and $O(M(n) \log q)$ for $\sin(\pi p/q)$ by the Chebyshev iteration.
order: 40
---

The elementary functions of `mp_real` cover precisions from one limb,
where an evaluation takes a few dozen nanoseconds, to hundreds of
millions of bits, where it takes minutes. No single algorithm is good
over that range, and the module combines several: table lookups
with short polynomials at a few limbs, argument reduction against
precomputed tables of logarithms and arctangents up to a few hundred
limbs, multi-prime reduction and binary splitting beyond, and table-free
methods for single evaluations at huge precision. The design follows
[@Joh2014c] (medium precision) and [@HJ2024] (high precision with
precomputation), with the classical algorithms of Brent
[@Bre1976agm] [@Bre1976b]; [@Mul2016] and [@BZ2011] are general
references.

Throughout, $B = 2^{\bits}$, $n$ is a number of limbs and $p$ a number
of bits.

## Layers {#layers}

Every function is evaluated in three layers.

1. **The ball layer** (`mp_real_exp_bits` and so on) evaluates the
   function at the exact midpoint $m$ of the argument and adds the
   effect of the radius $\rho$ through a bound on the derivative over the
   ball. It reduces $m$ to a fixed-point argument $v \in [0, 1)$
   and chooses the number of limbs of the kernel from the target
   precision, the leading zero bits of $v$ (to keep relative accuracy
   near zeros of the function) and the kernel's error bound.
2. **The kernel** computes $f(v)$ for $v \in [0, 1)$ with an error bound
   in ulps (see [](mp_real/fixed-point) for the contract). The kernel is
   chosen by size:

   | Limbs | exp, sin, cos | log | atan |
   |---|---|---|---|
   | 1–3 | [table-driven](#small) | table-driven | table-driven |
   | 4–32 | [dense tables](#dense) | dense tables | dense tables |
   | 33–600 | [bitwise](#bitwise) | bitwise | bitwise |
   | 601–65536 | [diophantine](#diophantine) | [Newton–Taylor](#newton) | Newton–Taylor |
   | above | [table-free](#notab) | Newton–Taylor | Newton–Taylor |

3. **The reduced series** evaluate the function at a small argument
   $t < 2^{-r}$ left by a reduction: [rectangular splitting](#series)
   for moderate $r$ and [bit-burst evaluation](#bit-burst) at high
   precision. They are also called directly when the argument itself is
   small enough, since their cost falls with $r$ while that of a
   table-based reduction does not.

## Argument reduction of the ball functions {#reduction}

**Exponential.** For $|m| \ge \log 2$, write $|m| = q\log 2 + t$ with
$t \in [0, \log 2)$. The quotient $q$ is computed from the top limbs of
$m$ and two limbs of $1/\log 2$ (stored as the fraction
$1/\log 2 - 1$), which gives $q$ or $q - 1$; one multiply-subtract by
the cached floor of $\log 2$ and a correction give $t$. For $m < 0$ the
quotient is rounded up and the remainder complemented, so that
$\exp(m) = 2^{-q}\exp(t')$ with $t' \in [0, \log 2)$ and no division.
The bound $|m| < 2^{\bits - 5}$ keeps the exponent of the result in the
safe range of [](algorithms/floating-point/mp-real-balls).

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
$\pi$ the sine loses relative accuracy, which shows in the radius.

**Tangent.** The tangent is reduced the same way and finished as
$\pm\tan v$ or $\mp 1/\tan v$. Near a zero or a pole, $v$ is small and
an absolute error is a relative error of the result, so the kernel runs
with the leading zero bits of $v$ as extra precision, repeating the
reduction at the higher precision when needed. This terminates because
$\pi$ has finite irrationality measure, at most $7.11$ [@ZZ2020]: a
$b$-bit dyadic $m$ satisfies $|m - k\pi/2| \gtrsim 2^{-7.2b}$ for large
$b$, so a repetition beyond $8b + 1024$ extra bits is never needed (and
is not made; the status is then 0).

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

## Series and rectangular splitting {#series}

A Taylor series $\sum_{k<N} a_k v^k$ with rational coefficients is
evaluated at a fixed-point $v$ by rectangular splitting [@Smi1989]
[@PS1973]: with blocks of $m$ terms,

$$
S_b = \sum_{i<m} a_{mb+i}\, v^i + v^m S_{b+1},
$$

the powers $v, \ldots, v^m$ are computed once ($m$ full products), every
term is one scalar multiply-add of a precomputed power by a one-limb
coefficient, and each block boundary costs one full product. The
coefficients $1/k!$ are generated as ratios of consecutive integers,
chained downwards so that one `mpn_divrem_1` closes each run of factors
that fills a limb; for $1/(2k + 1)$ the integers $D/(2k+1)$ for a group
of odd numbers with product $D$ are used, with one rescaling per group.

**Tapering.** Block $b$ is multiplied by $v^{mb} < 2^{-rmb}$ when it
enters the sum, so only its top $\lceil(\bits n + g - rmb)/\bits\rceil$
limbs matter. Computing each block in its shortened window reduces the
total cost to about $m$ products for the powers and $N/(3m)$ for the
block boundaries; with $m \approx \sqrt{N/1.5}$ this is about
$2\sqrt{N/3}$ full products, against $N/3$ for the tapered Horner scheme
(`_mp_real_series_tapered`), which is faster only below about a dozen
limbs. Errors committed in block $b \ge 1$ are damped by
$2^{-\min(g, rm)}$, and the final multiplication by $v$ (for the odd
functions) by $2^{-r}$, which is why the series is accurate to 2 ulps
except for exp, whose sum is the result.

**The tangent** has coefficients $T_{k+2}/(2k+3)!$ involving the
tangent numbers, which are not ratios of small integers. Their leading
denominators are small, though: in `_mp_real_series_rs_tan`, the terms
are grouped into chunks whose common denominator $Q_j$ fits $j + 1$
limbs, and the sum is carried as $Q_j$ times the partial sum, so that a
term costs $j + 1$ scalar passes and moving to the next chunk one
division by $Q_j/Q_{j-1}$.

For the sizes on the hot paths (arguments below $2^{-32}$ or $2^{-64}$
after a reduction, up to a few dozen limbs), the series are generated
straight-line code; a generic engine making the same calls in the same
order serves the other sizes with identical results.

## Table-driven kernels {#small}

At one to three limbs on 64-bit machines, the kernels are short
table-driven evaluations in registers [@Tan1989]:

- exp: $v = j_1 2^{-7} + j_2 2^{-14} (+ j_3 2^{-21}) + r$ with tables of
  $\exp(j 2^{-7})$, $\exp(j 2^{-14}) - 1$ and $\exp(j 2^{-21}) - 1$ (128
  entries each) and a Taylor polynomial in $r$;
- log1p: a multiplicative reduction by short factors $1 - D_k[j]$ over
  three or four levels of 7 bits, with tables of $-\log(1 - D_k[j])$,
  then the Taylor polynomial of $\log(1 + u)$;
- sin and cos: the pairs $(\sin c, 1 - \cos c)$ for $c = j 2^{-7}$ and
  $j 2^{-14}$ combined with the polynomials of the remainder by the
  addition formulas;
- atan: $\operatorname{atan}(v) = \operatorname{atan}(c) +
  \operatorname{atan}((v - c)/(1 + vc))$ with $c = j 2^{-9}$ (513
  entries) and the odd series of the quotient; for a ratio $v = Y/X$ the
  quotient is $(512Y - jX)/(512X + jY)$ with an exact numerator and
  denominator, so $\operatorname{atan}(Y/X)$ costs one division.

The polynomials taper their precision (one limb for the innermost terms,
then two, then three) and keep only the high parts of products. The
tables hold three limbs per entry, read at their top limbs at smaller
sizes.

## Dense tables {#dense}

From four to 32 limbs, the kernels reduce the argument seven bits per
level against tables, with multipliers that are exact single limbs, so
that applying a level costs one row of one-limb products. For exp,
$t \gets v - \log(2^s(1 + C_1[j_1]))$ with $j_1$ the top seven bits,
then $t \gets t - \log(1 + j_k 2^{-7k})$ with $j_k = \lfloor t 2^{7k}
\rfloor$ for the following levels (the tables hold floors, so $t$ stays
nonnegative), the series of $\exp(t)$, and the product of the factors
$1 + j_k 2^{-7k}$. Logarithms use the dual multiplicative reduction,
arctangents rotate the vector $(1, v)$ by factors $1 - i c_k$, and the
sine and cosine rotate the pair $(\cos t, \sin t)$ by $1 + i c_k$ and
absorb the scale $\prod (1 + c_k^2)$ by the double-angle formulas
$\cos v = (X^2 - Y^2)/N$, $\sin v = 2XY/N$ with $N = X^2 + Y^2$, so that
no square root is needed.

At four to ten limbs the code is generated straight-line code with
static tables, at 11 to 32 limbs a generic loop (unrolled to
add-with-carry chains at 11 to 16 limbs, with an `mulx`/`adcx`/`adox`
row in assembly where available). The deeper table levels at those
sizes are computed at run time when first needed, entry by entry, and
published in a process-wide cache by a compare-and-swap on a state byte,
so that no reader waits. The dense kernels are about 1.5 to 4 times
faster than the bitwise kernels at the same sizes.

## Bitwise reduction {#bitwise}

Up to 600 limbs, the argument is reduced by subtracting tabulated
values $L_i = \log(1 + 2^{-i})$ or $A_i = \operatorname{atan}(2^{-i})$,
in the manner of the BKM and CORDIC algorithms [@BKM1994] [@Mul2016],
but at multiple precision and with a series for the remainder. The
tables are computed once at the precision and index range in use and
cached.

::: lemma {#lem-greedy}
Let $c_0 > c_1 > \cdots > c_r > 0$ with $c_{i-1} < 2c_i$, and
$0 \le x < 2c_0$. The greedy reduction, which for $i = 0, \ldots, r$
subtracts $c_i$ from $x$ if $x \ge c_i$, leaves $0 \le x < c_r$, using
each $c_i$ at most once.
:::

::: proof
Before step $i$ we have $0 \le x < 2c_i$: this holds for $i = 0$, and
if $x \ge c_i$ the step leaves $x - c_i < c_i$, otherwise $x < c_i$; in
both cases $x < c_i < 2c_{i+1}$. After step $r$ the same argument gives
$x < c_r$.
:::

The logarithms satisfy the hypothesis, since
$2L_{i+1} = \log(1 + 2^{-i} + 2^{-2i-2}) > L_i$, and so do the
arctangents, by concavity of $\operatorname{atan}$.

**Exponential.** After the reduction $x = \sum_{i \in U} L_i + t$ with
$t < L_r < 2^{-r}$,
$$
\exp(x) = \exp(t) \prod_{i \in U} (1 + 2^{-i}),
$$
and each factor is one shift-and-add. When the series for $\exp(t)$ is
long, $\sinh t$ is computed instead (half the terms) and
$\exp(t) = \sinh t + \sqrt{1 + \sinh^2 t}$.

**Logarithm.** Dually, a product $P = \prod_{i \in U}(1 + 2^{-i})$ is
built greedily below $X = 1 + x$, each accepted factor one shift-and-add,
and
$$
\log X = \sum_{i \in U} L_i + 2\operatorname{atanh}\frac{X - P}{X + P}.
$$
The numerator is the exact deficit maintained through the reduction,
and the single division both normalises by $P$ and produces the atanh
argument, whose series has half the terms of that of $\log(1 + u)$.

**Arctangent.** The vector $(X, Y) = (1, x)$ is rotated by the factors
$1 - i 2^{-i}$ whenever $Y \ge 2^{-i} X$, that is whenever the angle
$\operatorname{atan}(Y/X) \ge A_i$. This is the greedy reduction of the
angle with steps $A_i$, so $\operatorname{atan}(Y/X) < 2^{-r}$ at the
end, and
$\operatorname{atan} x = \sum_{i \in U} A_i + \operatorname{atan}(Y/X)$
with one division. The decisions compare $Y$ with a shift of $X$ and
never read the table; the growth of $|X + iY|$ needs no compensation
because angles are scale invariant.

### Half-angle reconstruction {#half-angle}

For the sine and cosine, the angle $x/2$ is reduced greedily by the
$A_i$, $i \ge 1$ (by [](#lem-greedy), since $x/2 < 1/2 < 2A_1$), and the
used indices define $W = w_x + i w_y = \prod_{i \in U} (1 + i 2^{-i})$,
whose argument is $\sum_{i \in U} A_i$. Halving is needed because the
decisions are taken on a short model of the remainder that requires
table entries below $2^{-i}$; the doubled angles $2A_i$ would not
qualify.

::: theorem {#thm-half-angle}
Let $x/2 = \arg W + t'$ and $u = \tan t'$. With $T = w_y + w_x u$ and
$D = w_x - w_y u$,
$$
\sin x = \frac{2TD}{D^2 + T^2}, \qquad
\cos x = 1 - \frac{2T^2}{D^2 + T^2}, \qquad
\tan x = \frac{2TD}{D^2 - T^2}.
$$
:::

::: proof
By the addition formula, $\tan(x/2) = \tan(\arg W + t') = (w_y/w_x +
u)/(1 - u w_y/w_x) = T/D$. The formulas follow from
$\sin x = 2\tau/(1 + \tau^2)$, $\cos x = (1 - \tau^2)/(1 + \tau^2)$ and
$\tan x = 2\tau/(1 - \tau^2)$ with $\tau = T/D$, after multiplying
numerator and denominator by $D^2$.
:::

The quotient $T/D$ is never formed: each output costs one reciprocal and
one high product. The modulus $|W|$ cancels, so no scale factor
$\prod(1 + 2^{-2i})^{1/2}$ is needed, and in the limb version the
residual enters through $\sin t'$ and $g = 1 - \cos t'$ as
$T = w_y + w_x \sin t' - w_y g$ and $D = w_x - w_x g - w_y \sin t'$, the
common factor $\cos t'$ cancelling as well.

### Signed digits on vector lanes {#signed-digits}

From a few limbs on 64-bit machines, the bitwise reductions use signed
digits $d_i \in \{-1, 0, 1\}$: for the exponential,
$x = \sum_{d_i = 1} a_i - \sum_{d_i = -1} b_i + t$ with $a_i = L_i$
and $b_i = -\log(1 - 2^{-i})$, the digit $d_i$ chosen by comparing
$x$ with $\pm 0.75 \cdot 2^{-i}$. Since $a_i \approx b_i \approx 2^{-i}$,
the remainder satisfies $|x| \lesssim 2^{-i}$ before step $i$ for any
decision taken on an approximation of $x$ within about $0.2 \cdot
2^{-i}$ (to first order, a remainder below $0.75 \cdot 2^{-i}$ is kept,
and one above it drops by about $2^{-i}$; the tolerance of $0.25 \cdot
2^{-i}$ shrinks by the second-order terms). The decisions therefore run
on a 64-bit model of the remainder, refreshed every 52 indices from a
short head, and never branch on the full-precision operand; about $r/3$
digits are nonzero, against $r/2$ for the unsigned greedy reduction.

The full-precision operands are held in radix $2^{56}$ in the 64-bit
lanes of `vec4n` vectors (see [](machine_vectors)) in carry-save form:
every addition of a table entry and every shift-and-add of the
reconstruction is a lanewise vector operation, and the eight bits of
headroom per lane admit about 60 operations between carry
normalisations. The same scheme serves the logarithm (driving $1 + x$
multiplicatively to 1), the arctangent (signed rotations) and the
half-angle reconstruction. Measured on x86-64 with AVX2, this is 1.05
to 1.5 times faster than the limb version, and the cheaper reduction
moves the optimal $r$ up (to 768 at several hundred limbs).

The reduction parameter $r$ trades the table and reduction work against
the length of the remainder's series. Its default is tuned per function
and size: per size at four to ten limbs, where the candidates are the
table depths $7L - 1$, and from a ladder of values beyond (from 64 at a
dozen limbs to 768).

## Bit-burst evaluation {#bit-burst}

At high precision the Taylor series of the reduced argument has many
terms, and evaluating it with full-precision arithmetic costs
$O(\sqrt{N})$ full products at best. Brent's bit-burst algorithm
[@Bre1976b] instead splits the argument into slices of increasing
length,
$$
t = t_0 + t_1 + \cdots + t_K, \qquad t_k = \text{bits } 2^k\ldots 2^{k+1}
\text{ of } t \text{ (roughly)},
$$
and computes $\exp(t) = \prod_k \exp(t_k)$, where each $\exp(t_k)$ is
the sum of a series with a short rational argument, evaluated exactly
by binary splitting ([](algorithms/floating-point/hypergeometric-series)).
A slice of $2^k$ bits starting at bit $2^k$ needs about $p/2^k$ terms
whose integers grow by $2^k$ bits per term, so each slice costs
$O(M(p)\log p)$ and the whole $O(M(p)\log^2 p)$.

For the exponential the splitting computes
$T/(QB^{QE}) = \sum_k x^k/(k!\,B^{Dk})$ for a slice $x B^{-D}$. For the
sine and cosine, both are obtained from one tree over the series in
$y = x^2$:
$$
\sin x' = x'\,\frac{Q - B'}{Q}, \qquad 1 - \cos x' = \frac{A}{Q}, \qquad
-\frac{B'}{Q} = \sum_{k \ge 1}\frac{(-y)^k}{(2k+1)!}, \quad
-\frac{A}{Q} = \sum_{k \ge 1}\frac{(2k+1)(-y)^k}{(2k+1)!},
$$
using $1/(2k)! = (2k+1)/(2k+1)!$ so that both share the denominator
$Q = \prod (2i)(2i+1)$; the second track costs two extra products per
merge, half the work of a splitting of $\exp(ix)$ over all $2N$ terms.
Past a term threshold the cosine track is dropped and recovered from
the slice's denominator by one square root. The slice factors
$\exp(ix_k)$ are multiplied by `mp_real_mul_complex`.

These splittings have fixed midpoints, so every merge at a given level
multiplies by the same power $x^s$ of the argument; the powers are
tabulated, and their `fft_small` transforms are computed once per tree
for the entries used by at least two merges
([](algorithms/floating-point/hypergeometric-series#powtab)).

`_mp_real_exp_reduced` and `_mp_real_sin_cos_reduced` choose between the
series, one bit-burst step (the leading slice by binary splitting and
the rest by the series at the reduced size) and the full cascade, by
tuned thresholds that depend on the number of threads (a slice
parallelises where the series does not). The slices of a cascade run as
tasks on the available threads.

## Table-free evaluation {#notab}

A single evaluation at huge precision cannot amortise tables. The
table-free kernels halve the argument $h$ times, so that the reduced
argument is below $2^{-r}$ for a tuned depth $r$, evaluate the reduced
function by bit-burst, and undo the halving: the exponential by $h$
squarings, the sine and cosine by $h$ angle doublings on
$g = 1 - \cos$ alone,
$$
g(2\theta) = 2g(2 - g), \qquad \sin x = \sqrt{2g - g^2},
$$
with every intermediate a pure fraction and one square root at the end.

::: proposition {#prop-notab}
Squaring a value with relative error $\varepsilon$ gives relative error
at most $2\varepsilon + \varepsilon^2$; a doubling of $g$ with absolute
error $\delta$, for $g < 1/2$, gives an absolute error at most $4\delta +
2\delta^2$; and $\sqrt{2g - g^2}$ divides an absolute error of $g$ by
$2\sin x \cdot$ (to first order) $1/(1 - g)$.
:::

::: proof
$(y(1 + \varepsilon))^2 = y^2(1 + 2\varepsilon + \varepsilon^2)$. For
the doubling, $\frac{d}{dg} 2g(2 - g) = 4 - 4g \in (2, 4]$ and the
second derivative is $-4$. For the square root,
$\frac{d}{dg}\sqrt{2g - g^2} = (1 - g)/\sqrt{2g - g^2} = \cos x/\sin x$,
and $\cos x \le 1$.
:::

Hence $h$ squarings amplify the relative error by $2^h$, and $h$
doublings followed by the square root amplify the absolute error by at
most $2^{2h}/\sin x$; since any path that doubles at all has
$x \ge 2^{-z-1}$ for an argument with $z$ leading zero bits, $h + \log_2
n + 8$ resp. $2h + z + \log_2 n + 8$ guard bits keep the amplified
error below one output ulp, and the bound written is 128 ulps.

`mp_real_exp_notab_squaring` applies the same idea to balls, with $s$
ball squarings, and `mp_real_exp_agm` improves a starting value computed
this way by one Newton–Taylor step on the [AGM logarithm](#agm).

## Diophantine reduction {#diophantine}

From 600 to 65536 limbs, the exponential and the trigonometric
functions subtract integer combinations of logarithms of small primes
(resp. arguments of Gaussian primes) and multiply by a rational (resp.
Gaussian rational) number, which needs far cheaper tables than the
bitwise reduction at these sizes.

::: summary algorithms/floating-point/diophantine-reduction :::

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

## Sines of rational multiples of π {#chebyshev}

The values $\cos(\pi R/Q)$ are algebraic of degree at most $Q$, and at
high precision it is cheaper to refine one as a root of a polynomial
than to evaluate the cosine.

::: theorem Chebyshev step {#thm-chebyshev}
Let $\varphi^* = \pi R/Q \in (0, \pi/2)$ in lowest terms,
$x^* = \cos\varphi^*$, $\varepsilon = (-1)^R$, and let $x_0 = \cos\varphi_0$
with $\varphi_0 = \varphi^* + \delta$, $|Q\delta| < \pi/2$. Put
$C = T_Q(x_0)$, $A = x_0 T_Q(x_0) - T_{Q+1}(x_0)$, $a = \varepsilon A$,
$y = 1 - C^2$ and $b = 1/Q$. Then
$$
x^* = x_0 + x_0\, G(y) + \frac{a}{Q} F(y),
$$
$$
G(y) = {}_2F_1\bigl(-\tfrac b2, \tfrac b2; \tfrac12; y\bigr) - 1, \qquad
F(y) = {}_2F_1\bigl(\tfrac{1-b}{2}, \tfrac{1+b}{2}; \tfrac32; y\bigr).
$$
:::

::: proof
Since $T_k(\cos\varphi) = \cos(k\varphi)$, $C = \cos(Q\varphi^* +
Q\delta) = \varepsilon\cos(Q\delta)$, and
$x_0 T_Q(x_0) - T_{Q+1}(x_0) = \cos\varphi_0\cos(Q\varphi_0) -
\cos((Q+1)\varphi_0) = \sin\varphi_0 \sin(Q\varphi_0)$, so
$a = s_0\sigma$ with $s_0 = \sin\varphi_0$ and $\sigma = \sin(Q\delta)$,
and $y = \sigma^2$. Then $x^* = \cos(\varphi_0 - \delta) = x_0\cos\delta +
s_0\sin\delta$, and with $\delta = b\operatorname{asin}\sigma$ (valid on
the principal branch, $|Q\delta| < \pi/2$) the classical expansions
$\cos(b\operatorname{asin}\sigma) = {}_2F_1(-\frac b2, \frac b2; \frac12;
\sigma^2)$ and $\sin(b\operatorname{asin}\sigma) = b\sigma\,{}_2F_1(\frac{1-b}2,
\frac{1+b}2; \frac32; \sigma^2)$ give the claim, $s_0\sin\delta$
becoming $s_0\sigma\, b F(y) = (a/Q)F(y)$.
:::

No division or square root enters: $s_0$ appears only through $a$ and
$y$. The pair $(T_k(x_0), T_{k+1}(x_0))$ is computed by a Montgomery
ladder on $V_k = 2T_k$, $V_{2k} = V_k^2 - 2$, $V_{2k+1} = V_kV_{k+1} -
V_1$, one squaring and one product per bit of $Q$. The coefficients of
$G$ and $F$ decrease in magnitude from $1/(2Q^2)$ resp. 1, so for
$|y| \le 1/2$ the tails after $y^J$ are below $|y|^J/Q^2$ resp.
$2|y|^J$; in $y/(4Q^2)$ both have integer coefficients over common
denominators, so the series cost a few products and one division by a
word each. With $x_0$ correct to $e$ bits, $\sigma \approx Q2^{-e}$, and
summing to the target precision makes a step of any order $r$ from $x_0$
at about $(p + (r-1)\log_2 Q)/r$ bits; the order used is 12. The
recursion runs in ball arithmetic and checks the branch condition
$|Q\delta| < \pi/2$ from the rigorous radius of the previous level. The
cost is $O(M(n)\log Q)$, against the $O(M(n)\log^2 n)$ of the general
kernel; the iteration is used for denominators up to a bit length that
grows with the precision (none below 48 limbs, 12 bits at 512 limbs, 48
at 65536 limbs, and all denominators beyond). Of $\cos(\pi t/q)$ and
$\sin(\pi t/q) = \cos(\pi(q - 2t)/(2q))$ the one with the smaller
denominator is iterated and the other is $\sqrt{1 - x^{*2}}$.

## Use in FLINT {#use}

- The ball functions `mp_real_exp_bits`, `mp_real_log_bits`,
  `mp_real_atan_bits`, `mp_real_sin_cos_bits`, `mp_real_tan_bits`, the
  functions of $\pi x$, and the compositions `mp_real_log1p_bits`,
  `mp_real_expm1_bits`, `mp_real_sinh_cosh_bits`, `mp_real_tanh_bits`,
  `mp_real_asin_bits`, `mp_real_acos_bits`, `mp_real_asinh_bits`,
  `mp_real_acosh_bits` and `mp_real_atanh_bits` implement the layers
  above.
- The kernels are `_mp_real_exp_bitwise_rs` and its siblings,
  `_mp_real_exp_reduced`, `_mp_real_sin_cos_reduced`,
  `_mp_real_exp_notab`, `_mp_real_sin_cos_notab`,
  `_mp_real_exp_diophantine`, `_mp_real_sin_cos_diophantine`,
  `_mp_real_neglog_newton`, `_mp_real_atan_newton`,
  `_mp_real_neglog_agm` and `_mp_real_sin_cos_pi_ui_div_ui`; the
  table-free ball exponentials are `mp_real_exp_notab_squaring` and
  `mp_real_exp_agm`.
- `nfloat` calls the kernels and the ball functions for its elementary
  functions ([](nfloat)), and the partition function uses the table-free
  exponentials. The `arb` elementary functions have their own
  implementations of the same algorithms ([](algorithms/special-functions/elementary)).

::: functions mp_real_exp_bits mp_real_sin_cos_bits _mp_real_exp_bitwise_rs _mp_real_exp_reduced _mp_real_neglog_newton _mp_real_neglog_agm _mp_real_sin_cos_pi_ui_div_ui :::
