---
title: Table-based argument reduction
summary: The kernels of the elementary functions from one limb to 600 limbs - table-driven evaluation in registers, dense tables of seven-bit reduction steps with exact one-limb multipliers, and bitwise (BKM-type) reduction with signed digits on vector lanes and the half-angle reconstruction of the sine and cosine.
complexity: $O(n^2)$ word operations per evaluation at $n$ limbs with tables of $O(n)$ entries (bitwise: $O(r)$ table entries of $n$ limbs for a reduction depth $r$), plus a series of $O(n/r)$ terms.
order: 11
---

From one limb to a few hundred limbs, the kernels of the elementary
functions ([](algorithms/special-functions/elementary)) reduce their
argument against precomputed tables until a short series finishes the
evaluation: at one to three limbs from tables of a few hundred entries
in registers, at four to 32 limbs by seven bits per step against dense
tables, and up to 600 limbs by one bit per step against tables of
$\log(1 + 2^{-i})$ and $\operatorname{atan}(2^{-i})$. The remainder is
evaluated by the series of [](algorithms/special-functions/elementary-series).
Throughout, $B = 2^{\bits}$ and $n$ is a number of limbs.

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

## Use in FLINT {#use}

- The table-driven kernels of one to three limbs are inline functions
  (`mp_real/small.h`), also called directly by the register paths of
  `arf`, `nfloat`, `arb` and `acf` at one and two limbs.
- The dense kernels are `_mp_real_exp_dense_rt`,
  `_mp_real_sin_cos_dense_rt` and their siblings.
- The bitwise kernels are `_mp_real_exp_bitwise_rs`,
  `_mp_real_log1p_bitwise_rs`, `_mp_real_sin_cos_bitwise_rs`,
  `_mp_real_atan_bitwise_rs` and `_mp_real_tan_bitwise_rs`.

::: functions _mp_real_exp_dense_rt _mp_real_exp_bitwise_rs _mp_real_log1p_bitwise_rs _mp_real_sin_cos_bitwise_rs _mp_real_atan_bitwise_rs _mp_real_tan_bitwise_rs :::
