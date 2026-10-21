---
title: Elementary functions in double-word to quad-double precision
summary: Table-driven argument reduction, polynomials evaluated in several chains with tapered precision, and static error bounds for exp, log, sin, cos and atan on floating-point expansions, with the propagation of ball radii.
complexity: A constant number of $N$-term expansion operations per function value: one or two table lookups, a polynomial of degree 3 to 15 with tapered precision, and a few full-precision products.
order: 87
---

The `dfloat` module evaluates $\exp$, $\operatorname{expm1}$, $\log$,
$\operatorname{log1p}$, $\sin$, $\cos$ and $\operatorname{atan}$ on
expansions of $N = 1, \ldots, 4$ doubles (see
[](algorithms/floating-point/double-word)) as pure floating-point
kernels with *static* error bounds: each kernel is a fixed sequence of
expansion operations whose total error has been bounded once and for
all, over the whole argument range, by a model of every operation in
`dev/dfloat_exp_bound.py`. The ball versions then add the static bound
to the radius instead of tracking the rounding errors operation by
operation, and propagate the input radius with a bound on the
derivative. The other functions (`tan`, `sinh`, `cosh`, `tanh`, `asin`,
`acos`, `atan2`, `pow`, and the complex functions) are compositions of
these kernels with the arithmetic.

The design follows the classical table-driven methods for double
precision [@Tan1989; @Mul2016], adapted to expansions in three ways:

- **Exact argument reduction.** The reductions are arranged so that the
  cancelling operations on the head are exact (by Sterbenz's lemma and
  ulp counting), and the reduced argument keeps relative accuracy.
- **Tapered precision.** In Horner's rule for a polynomial in a small
  argument $w$, the inner coefficients are multiplied by powers of $w$,
  so their errors are scaled down; the inner steps are therefore done
  with fewer components (for example the schedule $[4, 3, 3, 2, 1, 1]$
  for the exponential at $N = 4$). The schedules are chosen by the bound
  script.
- **Several chains.** A polynomial is split into its even and odd parts,
  and these again in even and odd parts, so that four independent chains
  in $w^2$ (or $w^4$) run in the four lanes of a SIMD register; the final
  combination is a single five-term sum. This shortens the dependent
  chain of expansion operations, which determines the latency of a
  scalar call. The vector functions (`_dN_vec_exp`, ...) instead use the
  lanes for four different arguments.

Since tapered chains make the degree of a polynomial almost irrelevant
to the latency, what matters is the number of full-precision operations,
which a table lookup keeps at one or two. For $N = 1$ the kernels are
ordinary double code (one table level, polynomials by FMA) and about as
fast as the C library.

## Exponential {#exp}

Let $x = x_0 + x_1 + \cdots$ be canonicalised. With
$k = \operatorname{round}(x_0/\log 2)$, $L_0$ the 41-bit truncation of
$\log 2$ and $L_1, L_2, \ldots$ the following terms of its expansion,
the head is reduced to $h = x_0 - kL_0$, then to
$s_0 = h - i/256 - j/65536$ with $|i| \le 89$, $|j| \le 128$ and
$|s_0| \le 2^{-17}$.

::: lemma {#lem-exp-reduce}
For $|x_0| \le 746$, the three subtractions computing $h$, $d = h - i/256$
and $s_0 = d - j/65536$ are exact.
:::

::: proof
$kL_0$ is exact since $L_0$ has 41 bits and $|k| \le 1078 < 2^{11}$.
For $|k| \ge 2$, $x_0$ and $kL_0$ are within a factor 2 of each other
($|x_0 - k\log 2| \le 0.35$ and $|k| \log 2 \ge 1.38$), so the
difference is exact by Sterbenz's lemma; for $|k| = 1$ both are
multiples of $2^{-54}$ and their difference is below $1/2$, hence fits
in 53 bits; for $k = 0$ there is nothing to subtract. Similarly, $d$
and $s_0$ are differences of multiples of $\operatorname{ulp}(h)$ and
$\operatorname{ulp}(d)$ of magnitude at most $1/512$ and $2^{-17}$,
which is at most $2^{53}$ ulps when $|h| \ge 2^{-10}$ and
$|d| \ge 2^{-18}$, while $i = 0$ and $j = 0$ otherwise.
:::

The reduced argument
$s = s_0 + (x_1 + \cdots + x_{N-1}) - k(L_1 + \cdots + L_N)$ is
assembled with the expansion kernels in an order that keeps every term
significant. With $T = e^{i/256} e^{j/65536}$ (the product of two table
entries, 179 and 257 entries of 5-term expansions) and $w = s^2$,
$$e^x = 2^k\bigl(T + (Ts) E(w) + (Tw) O(w)\bigr),
\quad E(w) = \sum_j \frac{w^j}{(2j+1)!},\quad O(w) = \sum_j \frac{w^j}{(2j+2)!},$$
where $E$ and $O$ are the even and odd parts of the Taylor polynomial of
$\operatorname{expm1}$, of degree 3, 6, 8 and 11 for $N = 1, \ldots, 4$.
The scalar kernel splits $E$ and $O$ once more in $v = w^2$, so that the
four chains $E_0, E_1, O_0, O_1$ run in parallel and
$(Ts)E_0 + (Tsw)E_1 + (Tw)O_0 + (Tw^2)O_1$ is one lane-wise product and
one five-term sum. The result has relative error at most
`DFLOAT_EXP_EPS_N` (about $2^{-51}$, $2^{-100.2}$, $2^{-151}$ and
$2^{-201.6}$) for $-746 \le x_0 \le 709.9$.

For `expm1` with $|x| \le 0.34$ (where $k = 0$), the constant term $T$
is replaced by $T - 1 = E_1 + T_1 E_2$ from tables of $e^{i/256} - 1$
and $e^{j/65536} - 1$, so that the result is relative when $T - 1$ is
small or zero ($T - 1 = 0$ and $s = x$ exactly for $|x| < 2^{-17}$);
beyond 0.34 it is $e^x - 1$, which loses nothing since
$|e^x - 1| \ge 0.29$ there.

## Sine and cosine {#sin-cos}

With $k = \operatorname{round}(2x_0/\pi)$ and $P_0$ the 33-bit
truncation of $\pi/2$, the head is reduced to $h = x_0 - kP_0$; this is
exact for $|k| \le 2^{20}$ by the argument of [](#lem-exp-reduce), since
$kP_0$ then has at most 53 bits. The reduced argument
$t = h + (x_1 + \cdots) - k(P_1 + \cdots + P_N)$ is assembled and
canonicalised (so that its tail is bounded by $u|t|$ rather than
$u|x|$, which matters for large $x$) and split as
$t = i/64 + j/2^{13} + s$ with $|s| \le 2^{-14}$. With $A = i/64 + j/2^{13}$,
$S_A = \sin A$ and $C_A = \cos A$ from the products of two table levels
(51 and 65 entries, using the sign symmetry) and $w = s^2$,
$$\sin t = S_A + (S_A w) Q(w) + (C_A s) P(w), \qquad
\cos t = C_A + (C_A w) Q(w) - (S_A s) P(w),$$
where $\sin s = sP(s^2)$ and $\cos s = 1 + s^2 Q(s^2)$ are Taylor
polynomials. The results for $x$ are $\pm\sin t$ and $\pm\cos t$
according to $k \bmod 4$. As for the exponential, the four chains
$P_0, P_1, Q_0, Q_1$ in $w^2$ run in the SIMD lanes, as do the table
products and the multipliers.

Since the argument is reduced with $\pi/2$ to about $53N + 48$ bits and
the reduced argument has the absolute accuracy $u^N|x|$ of the input,
the error bounds `DFLOAT_TRIG_EPS_N` are absolute and scaled by $|x|$:
$\varepsilon(|x| + \min(|x|, 1))$ for the sine and $\varepsilon(|x| + 1)$
for the cosine. Relative accuracy near the zeros at multiples of
$\pi/2$ (other than 0) is not attempted. Arguments with
$|x| > 2^{20}\pi/2 - 1$, where $kP_0$ would no longer be exact, are
evaluated through `arb_t`.

## Logarithm {#log}

Write $x = 2^e m$ with $m \in [M_0, 2M_0)$, $M_0 = 0.70703125$, the
exponent and table index read from the bits of the head. The reduction
uses two tabulated reciprocals: $r \approx 1/m$ from 128 entries indexed
by the leading bits of $m$, giving $|mr - 1| \le 2^{-7}$ (with $r = 1$
on the two intervals adjacent to 1, which keeps $\log x$ relative near
$x = 1$), and $r' \approx 1/(1 + j/2^{13})$ with
$j = \operatorname{round}(2^{13}(mr - 1))$, giving
$z = m(rr') - 1$ with $|z| \le 2^{-14}(1 + 2^{-7})$. The product $rr'$
is exact as a double-word $R_h + R_l$, and
$z = (h_0 - 1) + l_0 + (m_1 + \cdots)R_h + mR_l$ with
$m_0 R_h = h_0 + l_0$ is computed with the head cancelled exactly, so
that its error is $u^N|z|$ rather than $u^N$; no division is needed.
The tables hold $L = -\log r$ for the rounded entries, so that
$$\log x = e\log 2 + L + L' + \operatorname{log1p}(z)$$
exactly. With $w = z^2$,
$\operatorname{log1p}(z) = z + wE(w) + zwO(w)$ with
$E(w) = -1/2 - w/4 - \cdots$ and $O(w) = 1/3 + w/5 + \cdots$ (degree 8,
8, 12 and 15 in $z$), evaluated as four chains in $w^2$. The relative
error bound `DFLOAT_LOG_EPS_N` holds for $2^{-1022} \le x < 2^{1022}$.

For `log1p` with $|x| \le 1/2$, the argument $1 + x$ is carried exactly
into the reduction: the head $\operatorname{RN}(1 + x_0)$ with its
rounding error and the tail of $x$ form an $(N + 1)$-term argument,
whose reduction is computed to $u^N|z|$ as above. For $|x| > 1/2$,
$1 + x$ is rounded to $N$ terms, and $|\log(1 + x)| \ge 0.4$ absorbs the
rounding.

A reduction by division, $z = (m - c)/(m + c)$ with
$\log m = \log c + 2\operatorname{atanh} z$, would halve the number of
coefficients by using only odd powers; but the even and odd parts of the
polynomial already run in parallel lanes, so nothing would be gained for
the price of a division, which costs much more than a reduction level by
multiplication.

## Arctangent {#atan}

For $x \ge 0$ (odd symmetry) and $c = i/64$: if $x \le 1$, with
$i = \operatorname{round}(64x)$,
$\operatorname{atan} x = \operatorname{atan} c + \operatorname{atan} z$,
$z = (x - c)/(1 + xc)$; if $x > 1$, with $i = \operatorname{round}(64/x)$,
$\operatorname{atan} x = (\pi/2 - \operatorname{atan} c) - \operatorname{atan} z$,
$z = (1 - cx)/(x + c)$. In both cases $|z| \le 2^{-7}$, the head of the
numerator is exact (Sterbenz), and the division is the long division
kernel; for $x \le 1/128$ no division is done. With $w = z^2$,
$\operatorname{atan} z = z + zwE(w)$, $E(w) = -1/3 + w/5 - \cdots$ with
3, 7, 11 and 14 coefficients, evaluated as four chains in $w^4$. One
division is unavoidable here; a second reduction level would save only
two tapered steps at the cost of a second division. The relative error bound `DFLOAT_ATAN_EPS_N` holds for all
finite $x$.

## Balls {#balls}

For a ball $[m \pm r]$ and the kernel value $y$ at the midpoint with
static relative bound $\varepsilon$:

- $\exp$: radius $|y|(1 + 2\varepsilon)\bigl((e^r - 1) + 2\varepsilon\bigr)$,
  since $|e^{x'} - e^m| \le e^m(e^r - 1)$ for $|x' - m| \le r$ and
  $e^m \le |y|/(1 - \varepsilon)$; $e^r - 1$ is bounded by $r + r^2$ for
  $r \le 1$ and by a power of two otherwise.
- $\log$: radius $\varepsilon|y| + r/(m - r)$, the derivative $1/t$
  being largest at the left end of the ball (a lower bound for $m - r$
  is computed from the head, the tail and the radius).
- $\operatorname{atan}$: radius $\varepsilon|y| + r$ (the derivative is
  at most 1).
- $\sin$, $\cos$: radius $\varepsilon(|m| + a) + r$ with
  $a = \min(|m|, 1)$ for the sine and $a = 1$ for the cosine; a radius
  of 1 or more is replaced by the enclosure $[0 \pm 1]$.

All radii are rounded up with the slack factor of
[](algorithms/floating-point/double-word#radius) and bumped by
$2^{-1069}$ when tiny. Arguments outside the kernel ranges, infinite
radii and nonfinite midpoints are handled through `arb_t`, so the ball
functions are rigorous for every input.

## Use in FLINT {#use}

The kernels are `dN_exp`, `dN_expm1`, `dN_sin`, `dN_cos`,
`dN_sin_cos`, `dN_log`, `dN_log1p` and `dN_atan` for $N = 1, \ldots, 4$,
their ball versions `dNb_exp`, etc., and the vector versions
`_dN_vec_exp`, `_dN_vec_sin_cos`, `_dN_vec_log`, `_dN_vec_atan`, which
run the same kernels on four arguments at a time (the static bounds
apply to every element). The tables are in `src/dfloat/exp_tables.c`,
generated by `dev/gen_dfloat_exp_tables.py`. The vector functions are
used by `gr_powsum_sieved` for the prime terms of the Riemann–Siegel
main sum ([](dfloat/powsum)).

::: functions d2_exp d2_sin_cos d2_log d2_atan d2b_exp _d2_vec_exp :::
