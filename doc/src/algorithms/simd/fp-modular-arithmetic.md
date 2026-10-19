---
title: Modular arithmetic in floating point
summary: Exact arithmetic modulo $n < 2^{50}$ on integers stored in doubles, using a precomputed $1/n$ and a fused multiply-add to split the product, as in FLINT's small-prime FFT.
complexity: A modular product costs two multiplications, two fused multiply-adds, one rounding and one addition, all vectorisable; a reduction costs one multiplication, one rounding and one fused multiply-add.
order: 10
---

Vector units multiply doubles much faster than they multiply 64-bit
integers: AVX2 has no instruction for a full $64 \times 64$-bit product,
while it performs four double-precision fused multiply-adds per
instruction. FLINT therefore does the modular arithmetic of its
small-prime FFT (`fft_small`) on integers stored exactly in doubles. The
functions are `vec4d_mulmod`, `vec4d_reduce_to_pm1n` and their
relatives in [`machine_vectors.h`](machine_vectors), and the ranges
stated there are proved below. The approach is of the kind studied by
van der Hoeven, Lecerf and Quintin [@HLQ2016]; exact linear algebra libraries
such as FFLAS have long relied on integer computations in floating point
[@DGP2008].

## Setting {#setting}

Doubles have a 53-bit significand. Let $u = 2^{-53}$ be the unit
roundoff and write $\mathrm{fl}(x)$ for $x$ rounded to the nearest
double, with ties to even. We assume round to nearest throughout and
ignore overflow and underflow, which cannot occur with the magnitudes
below. The standard facts we need [@MBD2018] are:

- Every integer $x$ with $|x| \le 2^{53}$ is a double. Consequently, an
  operation whose exact result is such an integer returns it exactly.
- For any real $x$ in range, $\mathrm{fl}(x) = x(1 + \delta)$ with
  $|\delta| \le u$.
- If $h = \mathrm{fl}(ab)$ for doubles $a$ and $b$, then $ab - h$ is a
  double, so a fused multiply-add computes it exactly:
  $\mathrm{fma}(a, b, -h) = ab - h$.

Fix an integer modulus $n \ge 1$ and let
$\mathit{ninv} = \mathrm{fl}(1/n) = (1 + \delta_0)/n$ with
$|\delta_0| \le u$. Residues are represented by integers in doubles, not
necessarily reduced: a typical range is $(-n, n)$, and intermediate sums
may grow to a few multiples of $n$.

## Rounding the quotient {#rounding}

Every function computes an approximate quotient $q \approx x/n$ and
returns $x - qn$. The quotient is an integer obtained by rounding the
product of a double $h$ with $\mathit{ninv}$. Two ways are used. The
`mulmod` and `nmulmod` functions, and all functions on NEON, always use
native rounding; `round_mul`, and hence `mulmod_fast` and the reductions
on AVX2, use the method selected per target by
`FLINT_FFT_SMALL_ROUND_USE_NATIVE`:

- *Native rounding*: $q = \mathrm{round}(\mathrm{fl}(h \cdot \mathit{ninv}))$
  with a rounding instruction (`vroundpd` on x86, `frintn` on ARM).
  This rounds twice.
- *Magic constant*: with $R = 3 \cdot 2^{51}$,
  $q = \mathrm{fl}(\mathrm{fl}(h \cdot \mathit{ninv} + R) - R)$, where
  the inner operation is one fused multiply-add. On Intel cores since
  Skylake the rounding instruction is slow, and this is faster.

::: lemma Magic constant {#lem-magic}
Let $x$ be a real number with $|x| \le 2^{51}$ and $R = 3 \cdot 2^{51}$.
Then $\mathrm{fl}(\mathrm{fl}(x + R) - R)$ is the integer nearest to $x$,
with ties to even.
:::

::: proof
The sum $x + R$ lies in $[2^{52}, 2^{53}]$, where consecutive doubles are
exactly 1 apart, so $s = \mathrm{fl}(x + R)$ is an integer nearest to
$x + R$. Since $R$ is an even integer, $s - R$ is an integer nearest to
$x$, and a tie is broken towards an even $s$, that is, an even $s - R$.
The subtraction is exact because $|s - R| \le 2^{51}$.
:::

For $|x| > 2^{51}$ the lemma fails: the spacing of the doubles near
$x + R$ exceeds 1 when $x > 2^{51}$, and is below 1 when $x < -2^{51}$, in
which case $s - R$ need not even be an integer. Native rounding has no
such limit. The difference between the two methods is in the error of
the quotient:

::: lemma {#lem-quotient}
Let $x$ and $y$ be real numbers with $y = x(1 + \epsilon)$.

1. If $q$ is the integer nearest to $y$, then
   $|x - q| \le 1/2 + |\epsilon x|$.
2. If $q$ is the integer nearest to $\mathrm{fl}(y)$, then
   $|x - q| \le 1/2 + |x| \left((1 + |\epsilon|)(1 + u) - 1\right)$.
:::

::: proof
In the first case $|x - q| \le |x - y| + |y - q| \le |\epsilon x| + 1/2$.
The second case is the first with $\mathrm{fl}(y) = x(1 + \epsilon)(1 + \delta)$,
$|\delta| \le u$.
:::

## Reduction {#reduction}

::: algorithm Reduction {#alg-reduce}
**Input:** an integer $a$ with $|a| < 2^{53}$, an integer $n \ge 4$ and
$\mathit{ninv} = \mathrm{fl}(1/n)$.

**Output:** an integer $r \equiv a \pmod n$ with $|r| < n/2 + 2$.

1. $q \gets$ the rounded value of $a \cdot \mathit{ninv}$, by either
   method.
2. Return $r = \mathrm{fma}(-q, n, a)$.
:::

This is `vec4d_reduce_to_pm1n` (and the identical
`vec4d_reduce_to_pm1no`). Following it with "add $n$ if negative" gives
`vec4d_reduce_to_0n`.

::: theorem {#thm-reduce}
[](#alg-reduce) is correct. More precisely, $|r| \le n/2 + u|a| < n/2 + 1$
with the magic constant and $|r| \le n/2 + (2u + u^2)|a| < n/2 + 2$
with native rounding.
:::

::: proof
With the magic constant, the fused multiply-add rounds the exact value
$a \cdot \mathit{ninv} + R$ once, and
$|a \cdot \mathit{ninv}| \le (2^{53} - 1)(1 + u)/4 < 2^{51}$ because
$n \ge 4$. So $q$ is the integer nearest to $y = a \cdot \mathit{ninv}$
by [](#lem-magic), and [](#lem-quotient) with $x = a/n$ and
$\epsilon = \delta_0$ gives $|a/n - q| \le 1/2 + u|a|/n$. With native
rounding, the product is rounded once more, and the second case of
[](#lem-quotient) gives $|a/n - q| \le 1/2 + (2u + u^2)|a|/n$.
Multiplying by $n$ gives the bounds on $a - qn$, and
$(2u + u^2)|a| \le (2u + u^2)(2^{53} - 1) < 2$. The integer $a - qn$
is thus smaller than $n/2 + 2 < 2^{53}$ in absolute value, so it is a
double and the fused multiply-add returns it exactly.
:::

The condition $n \ge 4$ is needed by the magic constant: for $n = 2$
or $3$ and $|a|$ close to $2^{53}$ the quotient exceeds $2^{51}$ and the
result is wrong. It also gives $|r| < n/2 + 2 \le n$, so the result lies
in $(-n, n)$, inside the range $[-n, n]$ promised by the name `pm1n`.

## Multiplication {#multiplication}

::: algorithm Modular multiplication {#alg-mulmod}
**Input:** integers $a, b$ stored as doubles with $|ab| < 2^{103}$,
an integer $1 \le n < 2^{50}$ and $\mathit{ninv} = \mathrm{fl}(1/n)$.

**Output:** an integer $r \equiv ab \pmod n$.

1. $h \gets \mathrm{fl}(ab)$.
2. $q \gets$ the rounded value of $h \cdot \mathit{ninv}$.
3. $\ell \gets \mathrm{fma}(a, b, -h)$, so that $h + \ell = ab$ exactly.
4. Return $r = \mathrm{fl}(\mathrm{fma}(-q, n, h) + \ell)$.
:::

`vec4d_mulmod` rounds natively in step 2 and `vec4d_mulmod_fast` uses
`round_mul`, which is the magic constant on most Intel targets. The
`nmulmod` variants compute $-r$ the same way: they set
$\ell = \mathrm{fma}(-a, b, h) = h - ab$ and return
$\ell - \mathrm{fma}(-q, n, h)$.

::: theorem {#thm-mulmod}
With native rounding, [](#alg-mulmod) returns an integer
$r \equiv ab \pmod n$ with

$$
|r| \le \frac{n}{2} + \eta |ab|, \qquad \eta = (1 + u)^3 - 1 < 3.0001\,u.
$$

With the magic constant, if moreover $|ab| \le (2^{51} - 1) n$, the same
holds with $\eta = (1 + u)^2 - 1 < 2.0001\,u$.
:::

::: proof
Write $x = ab/n$. With native rounding, $h = ab(1 + \delta_1)$ and
$\mathrm{fl}(h \cdot \mathit{ninv}) = x(1 + \delta_1)(1 + \delta_0)(1 + \delta_2)$
with all $|\delta_i| \le u$, so the second case of [](#lem-quotient)
gives $|x - q| \le 1/2 + \eta |x|$, that is,
$|ab - qn| \le n/2 + \eta|ab|$. With the magic constant, the fused
multiply-add rounds the exact value $y + R$ where
$y = h \cdot \mathit{ninv} = x(1 + \delta_1)(1 + \delta_0)$. By the
hypothesis, $|y| \le (2^{51} - 1)(1 + u)^2 < 2^{51}$, so $q$ is the
integer nearest to $y$ by [](#lem-magic), and the first case of
[](#lem-quotient) gives the bound with the smaller $\eta$.

It remains to see that $r$ is computed exactly. The double $h$ is an
integer, because $ab$ is an integer and every double of absolute value
at least $2^{53}$ is an integer. So $\ell = ab - h$ and $h - qn$ are
integers. From $|ab| < 2^{103}$ we get $|\ell| \le u|ab| < 2^{50}$ and
$|ab - qn| \le 2^{49} + 3.0001 \cdot 2^{50}$, so
$|h - qn| \le |ab - qn| + |\ell| < 5 \cdot 2^{50} < 2^{53}$. Hence the
fused multiply-add in step 4 returns $h - qn$ exactly, and the final
addition returns the integer $(h - qn) + \ell = ab - qn$ exactly.
:::

The theorem holds for any $n < 2^{50}$, prime or not, odd or even. The
ranges used by the transforms follow from it.

::: corollary {#cor-ranges}
Let $n < 2^{49}$. With either rounding method, $|ab| < 2n^2$ implies
$|r| < n$, and with native rounding $|ab| < 4n^2$ implies $|r| < 3n/2$.
:::

::: proof
With $|ab| < 2n^2$ and $n < 2^{49}$,
$\eta|ab| < 3.0001 \cdot 2^{-53} \cdot 2 \cdot 2^{49} n < 0.376\,n$, so
$|r| < 0.876\,n$. With $|ab| < 4n^2$ the same computation gives
$|r| < 1.251\,n$. The magic constant applies to the first case because
$|ab| < 2n^2 < (2^{51} - 1) n$.
:::

For $2^{49} \le n < 2^{50}$ the theorem gives only $|r| < 1.25\,n$ for
$|ab| < 2n^2$, and a finer analysis of the rounding in step 2, in the
comments of `fft_small/mulmod_satisfies_bounds.c`, gives $9n/8$. FLINT
therefore checks each 50-bit prime individually with
`fft_small_mulmod_satisfies_bounds`, which evaluates that finer bound
with the actual error $|n \cdot \mathit{ninv} - 1|$ in place of its
maximum $u$. A prime passes if the bound gives $|r| < n$ for
$|ab| < 2n^2$ and $|r| < 3n/2$ for $|ab| < 4n^2$, and
`sd_fft_ctx_init_prime` throws for a prime that fails. According to the
same comments about one random 50-bit prime in 900 fails, and none of 49
bits or less, in agreement with [](#cor-ranges).

::: warning
`mulmod_fast` with the magic constant gives wrong results when
$|ab| > 2^{51} n$. For example `mulmod_fast(-(3n - 1), n - 1)` with
$n = 1086317488242689$ returns a non-integer. With native rounding, and
always on ARM, `mulmod_fast` is `mulmod` and such inputs go unnoticed.
Defining `FLINT_ROUND_MUL_USE_MAGIC` selects the magic constant on any
x86-64 target, for testing.
:::

## Cheap reductions {#cheap}

Between multiplications, values are brought back into a smaller range by
comparisons and conditional additions, which are exact on integers below
$2^{53}$:

| Function | Input | Output | Method |
|---|---|---|---|
| `reduce_pm1no_to_0n` | $(-n, n)$ | $[0, n)$ | add $n$ if negative |
| `reduce_2n_to_n` | $[0, 2n)$ | $[0, n)$ | subtract $n$ if not negative afterwards |
| `reduce_0n_to_pmhn` | $[0, n]$ | $(-n/2, n/2]$ | subtract $n$ if above $n/2$ |
| `reduce_pm1n_to_pmhn` | $[-3n/2, 3n/2]$ | $[-n/2, n/2]$ | add or subtract $n$ if $\lvert a \rvert > n/2$ |

The vectorised versions use `blendv` on a sign bit or a comparison mask
instead of branches. Testing the sign bit treats $-0.0$ as negative, so
`vec4d_reduce_pm1no_to_0n` on AVX2 maps $-0.0$ to $n$, outside the stated
output range; the result is still congruent to the input.

## Use in FLINT {#use}

- **Small-prime FFT.** `fft_small` computes number-theoretic transforms
  modulo primes $p < 2^{50}$ with $p \equiv 1 \pmod{2^k}$, with all
  arithmetic in `vec8d` and `vec4d`. Twiddle factors are stored in
  $(-p/2, p/2)$, the butterflies use `mulmod_fast` and `reduce_to_pm1n`,
  and the bound analysis in `fft_small/mulmod_satisfies_bounds.c` shows
  that the values stay within a few multiples of $p$, so that all sums are
  exact. That analysis needs the actual ranges of
  [](#thm-reduce) and [](#thm-mulmod) in addition to the documented ones
  for the truncated inverse transform; it is carried out in
  [](algorithms/fft/small-prime-fft#arithmetic).
- **Polynomial and integer multiplication.** The conversions in and out of
  the transform use `mulmod`, `reduce_to_0n` and `reduce_0n_to_pmhn`. Through
  `fft_small`, this arithmetic underlies the multiplication of large
  integers (`flint_mpn_mul`) and of large polynomials in `nmod_poly` and
  `fmpz_poly`.
- **Precomputed transforms.** `mpn_mod` and the `gr` implementations of
  transformed polynomials multiply by fixed scalars with `vec8d_mulmod`.
- **Matrices.** `nmod_mat_mul_u8` reduces its results with
  `vec4d_reduce_to_0n` on x86. The floating-point kernel
  `nmod_mat_mul_fp50` uses its own versions of `round_mul`, `mulmod_fast`
  and `reduce_to_pm1n` with a fixed lane count per instruction set
  (`nmod_mat/mul_fp_vec.h`), with the magic constant except on NEON, and
  refers to the same analysis.
