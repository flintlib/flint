---
title: High, low and middle products of integers
summary: Computing only the low half, an approximation of the high half, or an arbitrary window of the limbs of an integer product, with hard-coded kernels, Mulders' recursive short products and proved error bounds.
complexity: About half the time of the full product for quadratic multiplication; about $0.81 M(n)$ with Karatsuba (Mulders); $M(n)$ for the balanced middle product by transposition.
order: 15
---

Many algorithms need only part of a product of two integers: Hensel
lifting and 2-adic Newton iteration need the low half, floating-point
multiplication and Newton division need the high half, and Newton
iteration for the reciprocal needs a middle part of a $2n \times n$
product whose low part is known to cancel. Computing only the required
part saves up to half of the work. The interfaces are described in
[](mpn_extras/mulhigh). FLINT's high products, their error analysis and
their kernels for small sizes are due to Ahlbäck and Johansson
[@AJ2025].

Throughout, $B = 2^{\bits}$ and $u, v$ are $n$-limb integers with limbs
$u_i, v_j$, so that $uv = \sum_{d=0}^{2n-2} \sum_{i+j=d} u_i v_j B^d$.
The partial products $u_i v_j$ with $i + j = d$ form the *diagonal* $d$.

## Low products {#low}

The low product $uv \bmod B^n$ needs only the diagonals $d < n$, that is,
$n(n+1)/2$ of the $n^2$ limb products. In the schoolbook algorithm one
simply truncates each row, and the high halves of the products on the
diagonal $n - 1$ are not needed either.

::: algorithm Mulders' short product {#alg-mulders-low}
**Input:** $n$-limb $u, v$, a split point $k$ with $n/2 \le k \le n$.

**Output:** $uv \bmod B^n$.

1. Write $u = u_1 B^k + u_0$ and $v = v_1 B^k + v_0$ with $u_0, v_0 < B^k$.
2. $p \gets u_0 v_0$, a full $k \times k$ product.
3. $s \gets (u_1 \cdot v_0 + u_0 \cdot v_1) \bmod B^{n-k}$, computed as two
   recursive short products of length $n - k$ of $u_1$ with the low
   $n - k$ limbs of $v_0$ and of $v_1$ with the low $n - k$ limbs of $u_0$.
4. Return $(p + s B^k) \bmod B^n$.
:::

::: theorem {#thm-mulders-low}
[](#alg-mulders-low) is correct. If a full $n \times n$ product costs
$M(n) = n^\alpha$ and $k = \beta n$, then the cost $S(n)$ of the short
product satisfies $S(n) \approx c\, M(n)$ with
$c = \beta^\alpha / (1 - 2(1 - \beta)^\alpha)$.
:::

::: proof
$uv = u_1 v_1 B^{2k} + (u_1 v_0 + u_0 v_1) B^k + u_0 v_0$, and
$B^{2k} \equiv 0 \pmod{B^n}$ since $2k \ge n$. Modulo $B^n$, the middle
term only depends on $(u_1 v_0 + u_0 v_1) \bmod B^{n-k}$, and $u_1 v_0
\bmod B^{n-k}$ only on $u_1$ and $v_0 \bmod B^{n-k}$ ($u_1$ has $n - k$
limbs). For the cost, $S(n) = M(k) + 2S(n - k)$; substituting
$S(n) = c n^\alpha$ gives $c = \beta^\alpha + 2c(1-\beta)^\alpha$.
:::

For schoolbook multiplication ($\alpha = 2$) the best choice is
$\beta = 1/2$, with $c = 1/2$; for Karatsuba ($\alpha = \log_2 3$) Mulders
found $\beta \approx 0.694$ and $c \approx 0.808$ [@Mul2000; @HZ2011]. In
practice the best split depends on the implementation of the full
product, and FLINT reads $k$ from a table tuned for each $n$ below 2048
(`flint_mpn_mulhigh_k_tab`), and uses $k = 3n/4$ above.

## High products {#high}

The high half $\lfloor uv / B^n \rfloor$ depends on all diagonals through
the carries, so it cannot be computed exactly faster than the full
product in the worst case. A high product computes the diagonals
$d \ge n - 1$ (or $d \ge n - 2$), which determine the top $n$ limbs up to
a small error, and drops the others. Following Ahlbäck and Johansson
[@AJ2025], write

$$
H'_n(u, v) = \sum_{i + j \ge n - 1} u_i v_j B^{i+j}, \qquad
H''_n(u, v) = H'_n(u, v) + B^{n-1} \sum_{i + j = n - 2}
  \left\lfloor \frac{u_i v_j}{B} \right\rfloor .
$$

Both use $n(n+1)/2$ limb products ($H''_n$ also needs the high halves of
$n - 1$ more), both are multiples of $B^{n-1}$ and so have $n + 1$
significant limbs. $H'_n$ is the classical truncated product, as used for
instance by MPFR's `mpfr_mulhigh_n`; $H''_n$ is what FLINT computes.

::: theorem {#thm-hprime}
$H'_1(u, v) = H''_1(u, v) = uv$, and for $n \ge 2$,
$uv - (n - 1) B^n < H'_n(u, v) \le uv$ and
$uv - (2n - 3) B^{n-1} < H''_n(u, v) \le uv$ [@AJ2025, Theorem 1].
:::

::: proof
The upper bounds hold because only nonnegative parts of $uv$ are
summed. For $H'_n$, the dropped part is
$\sum_{d \le n - 2} \sum_{i+j=d} u_i v_j B^d \le
\sum_{d=0}^{n-2} (d+1)(B-1)^2 B^d \le (n-1)(B-1)(B^{n-1} - 1) <
(n-1)B^n$. The bound for $H''_n$ is proved in
[](#thm-high).
:::

So $H'_n$ may be wrong by up to $n - 1$ units in the lowest of the top
$n$ limbs, and a caller that needs $n$ almost exact limbs must work with
$n + 1$ limbs (MPFR pads its operands by one limb), while $H''_n$ errs
by less than $2n - 3$ units of the limb below, about $B/2$ times less,
for the price
of $n - 1$ high halves of limb products.

::: algorithm Precise basecase high product {#alg-high}
**Input:** $n$-limb $u, v$ with $n \ge 2$.

**Output:** an $(n+1)$-limb $h$ approximating $\lfloor uv/B^{n-1} \rfloor$.

1. $t \gets \sum_{i + j = n - 2} \lfloor u_i v_j / B \rfloor$ (the high
   halves of the products on diagonal $n - 2$).
2. $C \gets t B^{n-1} + \sum_{i + j \ge n - 1} u_i v_j B^{i+j}$.
3. Return $h = C / B^{n-1}$.
:::

::: theorem {#thm-high}
The value $h$ of [](#alg-high) satisfies
$\lfloor uv/B^{n-1} \rfloor - (2n - 4) \le h \le \lfloor uv/B^{n-1} \rfloor$.
:::

::: proof
Every term of $C$ is a multiple of $B^{n-1}$ and $C \le uv$, so $h$ is an
integer and $h \le \lfloor uv / B^{n-1} \rfloor$. The missing amount is
$E = uv - C = \sum_{i+j = n-2} (u_i v_j \bmod B) B^{n-2} +
\sum_{d \le n - 3} \sum_{i+j=d} u_i v_j B^d$. The first sum has $n - 1$
terms, each less than $B$, so it is less than $(n-1) B^{n-1}$. Diagonal
$d$ has $d + 1$ terms, each at most $(B-1)^2$, so the second sum is at
most

$$
\sum_{d=0}^{n-3} (d+1)(B-1)^2 B^d \le (n-2)(B-1)^2 \frac{B^{n-2} - 1}{B - 1}
< (n-2) B^{n-1}.
$$

Hence $E < (2n - 3) B^{n-1}$, and since $C$ is a multiple of $B^{n-1}$,
$\lfloor uv/B^{n-1} \rfloor - h = \lfloor E / B^{n-1} \rfloor \le 2n - 4$.
:::

All algorithms used by `flint_mpn_mulhigh_n` and `flint_mpn_sqrhigh`
return an $h$ with the same bounds ([](#thm-recursive-high) below).

::: corollary {#cor-exact}
Write $h = H B + r$ with $0 \le r < B$, so that $H$ is the computed
approximation of $\lfloor uv / B^n \rfloor$. If $r \le B - 2n + 3$, then
$H = \lfloor uv/B^n \rfloor$.
:::

::: proof
$\lfloor uv/B^n \rfloor = \lfloor (h + e)/B \rfloor$ with
$0 \le e \le 2n - 4$, and $r + e < B$ under the assumption.
:::

So a high product almost always gives the exact top limbs, and the rare
exceptions are detected from the returned limb [@AJ2025, Corollary 4]:
for random inputs the test fails with probability about $2n/B$, which is
$3 \cdot 10^{-18}$ for $n = 32$ and $B = 2^{64}$. Interpreting $u$ and $v$ as fixed-point numbers in
$[0, 1)$ with $n$ limbs, the top $n$ limbs of $H_n'$ have an error of
less than $n$ units in the last place, and those of $H_n''$ an error of
less than 2 units (when $2n < B$) [@AJ2025, Corollary 2].

### Normalised high products {#normalised}

For floating-point multiplication, the mantissas are normalised (top bit
set, $u, v \ge B^n/2$), so that $B^{2n}/4 \le uv < B^{2n}$ and the
product has its top bit in one of two positions.
`flint_mpn_mulhigh_normalised` and `flint_mpn_sqrhigh_normalised`
compute the high product and, when the top bit of the result is zero,
shift the $n + 1$ limbs left by one bit; they return the low limb and the
shift (0 or 1), which the caller subtracts from the exponent. The shift
doubles the error measured in units of the last place, so the normalised
$n$-limb result has an error of less than 2 units in the last place when
$4n < B$ [@AJ2025, Corollary 3]. On x86-64 with ADX, these functions
have assembly routines for up to 9 limbs (squares up to 8); otherwise
they call the plain high product and shift.

### Recursive high products {#recursive-high}

Above the hard-coded sizes, $H''_n$ is assembled from smaller products
by exact identities [@AJ2025, Theorem 5]. Write $\mathrm{hi}(x) =
\lfloor x/B \rfloor$ for a two-limb $x$.

::: lemma Peeling and balanced splitting {#lem-split-high}
Let $n \ge 3$. Write $u = u_0 + B u'$ and $v = v' + B^{n-1} v_{n-1}$ with
$u', v'$ of $n - 1$ limbs. Then

$$
H''_n(u, v) = B\, H''_{n-1}(u', v') + B^n u' v_{n-1} + B^{n-1}
\bigl(u_0 v_{n-1} + \mathrm{hi}(u_0 v_{n-2})\bigr).
$$

If $n = p + q$ with $p, q \ge 1$, write $u = u_0 + B^p u_1$ and
$v = v_0 + B^q v_1$ with $u_0$ and $v_1$ of $p$ limbs and $u_1$ and $v_0$
of $q$ limbs. Then

$$
H''_n(u, v) = B^{n} u_1 v_1 + B^p H''_q(u_1, v_0) +
B^q H''_p(u_0, v_1) + B^{n-1} \mathrm{hi}(u_{p-1} v_{q-1}).
$$
:::

::: proof
Both sides are sums of terms $u_i v_j B^{i+j}$ and $B^{n-1}
\mathrm{hi}(u_i v_j)$; one checks that each pair $(i, j)$ with
$i + j \ge n - 1$ appears exactly once with its full product and each
pair with $i + j = n - 2$ exactly once with its high half. In the first
identity, $H''_{n-1}(u', v')$ covers the pairs with $i \ge 1$ and
$j \le n - 2$ (its diagonals $d' \ge n - 2$ and $n - 3$ are the diagonals
$d = d' + 1 \ge n - 1$ and $n - 2$), the second term the pairs with
$i \ge 1$ and $j = n - 1$, and the last the pairs $(0, n-1)$ and
$(0, n-2)$. In the second, $u_1 v_1$ covers $i \ge p$, $j \ge q$ (all
on diagonals $\ge n$), $H''_q(u_1, v_0)$ the pairs with $i \ge p$,
$j < q$ (its diagonals $d' \ge q - 1$ and $q - 2$ are $d = d' + p \ge
n - 1$ and $n - 2$), $H''_p(u_0, v_1)$ symmetrically the pairs with
$i < p$, $j \ge q$, and the pairs with $i < p$, $j < q$ reach diagonal
$n - 2$ only at $(p-1, q-1)$.
:::

In Mulders' recursion for the high product, the *top* $k$ limbs of $u$
and $v$ give a full product and the two cross terms are recursive high
products of length $\ell = n - k$, mirroring [](#alg-mulders-low).

::: algorithm Mulders' precise high product {#alg-mulders-high}
**Input:** $n$-limb $u, v$, a split point $k$ with $(n + 4)/2 \le k < n$;
let $\ell = n - k$.

**Output:** an $(n+1)$-limb $\tilde h$ approximating
$\lfloor uv/B^{n-1} \rfloor$.

1. Write $u = \hat u_0 + B^\ell \hat u_1$, $v = \hat v_0 + B^\ell \hat v_1$
   with $\ell$-limb $\hat u_0, \hat v_0$ and $k$-limb $\hat u_1, \hat v_1$.
2. $F \gets \lfloor B^{2\ell} \hat u_1 \hat v_1 / B^{n-1} \rfloor$ (a full
   $k \times k$ product, of which the limbs below $B^{n-1}$ are dropped).
3. $h_1 \gets$ the recursive high product of the top $\ell$ limbs of $u$
   and the bottom $\ell$ limbs of $v$; $h_2 \gets$ that of the bottom
   $\ell$ limbs of $u$ and the top $\ell$ limbs of $v$.
4. Return $\tilde h = F + h_1 + h_2 + \mathrm{hi}(u_{k-1} v_{\ell-1}) +
   \mathrm{hi}(u_{\ell-1} v_{k-1})$.
:::

The last two terms are the products of diagonal $n - 2$ that lie in
none of the three subproducts; they are what turns Mulders' original
algorithm, which computes $H'_n$ in the same way, into an approximation of
$H''_n$ [@AJ2025, Theorem 6].

::: theorem {#thm-recursive-high}
Let $C = H''_n(u, v)$, and call an $(n+1)$-limb value $\tilde h$
*admissible* if $C \le B^{n-1} \tilde h \le uv$. The exact $H''_n / B^{n-1}$
is admissible, and so is $\lfloor uv/B^{n-1} \rfloor$ (from a full
product). If the recursive high products in [](#alg-mulders-high) or in
the identities of [](#lem-split-high) are admissible, so is the result.
In particular, [](#thm-high) and [](#cor-exact) hold for every
combination of these methods.
:::

::: proof
For the full product, $\lfloor uv/B^{n-1} \rfloor B^{n-1} \ge C$ because
$C$ is a multiple of $B^{n-1}$ and $C \le uv$. For [](#lem-split-high),
replacing the exact $H''$ terms by admissible values keeps both bounds,
since the identities are sums of nonnegative terms. For
[](#alg-mulders-high): the pairs $(i, j)$ with
$i + j \ge n - 2$ split into those with $i, j \ge \ell$ (in
$\hat u_1 \hat v_1$), those with $i \ge k$, $j < \ell$ (in $h_1$), those
with $i < \ell$, $j \ge k$ (in $h_2$), and the two pairs
$(k-1, \ell-1)$ and $(\ell-1, k-1)$, since $k \ge (n+4)/2$ implies that
all other pairs lie on diagonals below $n - 2$. For the first group,
$B^{n-1} F \le B^{2\ell} \hat u_1 \hat v_1$, and $B^{n-1} F$ is at least
the contribution of these pairs to $C$, because the floor of a sum is at
least the sum of the floors. For $h_1$, admissibility for the
subproduct (which is shifted by $B^k$, and $B^k B^{\ell - 1} =
B^{n-1}$) gives a value between the contribution of its pairs to $C$
and their exact product, and likewise for $h_2$; the two corrections are
exactly their contributions to $C$ and at most $u_i v_j B^{n-2}$.
Summing gives $C \le B^{n-1} \tilde h \le uv$.
:::

`flint_mpn_sqrhigh` uses the same recursion for squares, where the two
cross terms coincide: $\tilde h = F + 2 h_1 + 2\,\mathrm{hi}(u_{k-1}
u_{\ell-1})$, with a full square of the top $k$ limbs and a single high
product, and the theorem holds unchanged. FLINT tunes $k$ for each
$n < 2048$ (`flint_mpn_mulhigh_k_tab` for products, a separate table for
squares) and uses $k = 3\lfloor n/4 \rfloor$, respectively
$\lfloor (n + 4)/2 \rfloor$, above; a table entry $k = n$ means a full
product, $k = 0$ the basecase.

## Windowed middle products {#window}

The windowed middle product computes limbs $[\ell, h)$ of a product of
operands of lengths $m$ and $n$, dropping the diagonals below $\ell$.

::: theorem {#thm-window}
Let $W = \sum_{i+j \ge \ell} u_i v_j B^{i+j-\ell}$ be computed from the
diagonals $d \ge \ell$ only. Then
$\lfloor uv/B^\ell \rfloor - \mu B < W \le \lfloor uv/B^\ell \rfloor$
with $\mu = \min(m, n, \ell)$.
:::

::: proof
The dropped part is $E = \sum_{d < \ell} \sum_{i+j=d} u_i v_j B^d$. The
diagonal $d$ has at most $\min(m, n, d+1) \le \mu$ terms, so
$E \le \mu (B - 1)^2 (B^\ell - 1)/(B - 1) < \mu B^{\ell + 1}$, and
$\lfloor uv/B^\ell \rfloor - W = \lfloor E / B^\ell \rfloor < \mu B$.
:::

Reducing $W$ modulo $B^{h - \ell}$ gives the window; an error can only
reach the limbs above the lowest two through a carry. FLINT's
`flint_mpn_mulmid` chooses between the schoolbook window (all diagonals
in $[\ell, h)$, at cost about $(h - \ell)\min(m, n)$ limb products), a low
or high product when the window is at one end of the product, a full
product, and two transform-based methods:

- **The transposition principle.** The balanced middle product of a
  $(2n - 1)$-limb and an $n$-limb operand (the diagonals $n - 1$ to
  $2n - 2$, the window needed by a Newton step for the reciprocal) is the
  transpose of a full $n \times n$ product, and can be computed by any
  bilinear multiplication algorithm in the same number of operations
  [@HQZ2004]. GMP's internal `mpn_mulmid_n` does so with Karatsuba's
  method; FLINT uses it when available.
- **Cyclic convolution.** With the small-prime FFT of `fft_small`, the
  operands are split into coefficients and multiplied modulo $x^L - 1$,
  which wraps the coefficient of index $c \ge L$ onto $c - L$. For a
  product with coefficients of indices $[0, N)$ and a window $[\ell, h)$
  (in coefficients), a cyclic length $L \ge \max(h, N - \ell)$ suffices:
  the wrapped coefficients then land below the window. For the balanced
  middle product, $L$ is about $2n$ instead of $3n$. Only the
  coefficients that reach the window are reconstructed from their
  residues.

### Middle products by transposed small-prime FFTs {#transposed-fft}

`fft_small` applies the transposition principle [@BLS2003] to its own
FFT product
([](algorithms/polynomials/middle-product#thm-transposed-fft)), as
[@HQZ2004] do for the FFT middle product of power series. On the
polynomials of $\beta$-bit chunks of `_mpn_ctx_mpn_mul_range`, let the
window need the chunk coefficients $[c_0, c_0 + m)$ of the product and
let the shorter operand have $r$ chunks. These coefficients are
$\sum_i \rho_i v_{i+j}$, $0 \le j < m$, with $\rho$ the chunks of the
shorter operand reversed and $v$ those of the longer operand from index
$c_0 - (r - 1)$ on: the transpose of a product of lengths $r$ and $m$,
which the transposed transforms compute with transforms of length about
$r + m$ instead of the length of the full product. With $c_0 \beta$ and
$(r - 1)\beta$ rounded to multiples of 64 (and $c_0$ to a multiple of 8
on the two-prime path), $v$ starts on a limb boundary and the
reconstruction sees the same coefficients at the same bit positions as
for the full product, so the result is the same, bit for bit, and
[](#thm-window) applies unchanged. `_mpn_ctx_mpn_mul_range`, and so
`flint_mpn_mulmid` above its `fft_small` threshold, uses this when the
window has fewer chunks than the longer operand and the transform is
strictly shorter (than the full product on the two-prime path up to
40000 limbs, which has no cyclic shortcut, and than the cyclic length
above), except for squarings. For the balanced middle product above, the
transform shrinks from about $3n$ to about $2n$ limbs' worth of chunks
for every $n$, not only when a power of two fits. `radix_mulmid_fft_small`
does the same for products in radix $b^e$
([](radix)), and `flint_mpn_mulhigh_n` and `flint_mpn_mullow_n`, whose
windows are longer than both operands, never use it.

## Known low limbs {#known-low}

In Hensel division and 2-adic square roots, a product $xy$ is needed
whose low $\ell$ limbs are known in advance (they equal a part of the
dividend). Then the window $[\ell, h)$ is computed exactly as follows: a
windowed product starting at $\ell - 3$ has, by [](#thm-window), a deficit
of less than $\mu B < B^2$ units of limb $\ell - 3$; comparing its lowest
three limbs with the known limbs $\ell - 3$, $\ell - 2$, $\ell - 1$ tells
exactly which carry into limb $\ell$ was lost. This is
`_flint_mpn_mulhigh_known_low`, which turns the high half of a Hensel
step into a high product at half the cost of a full product.

## Algorithm selection and performance {#selection}

The small sizes use the techniques of
[](algorithms/integers/basecase-multiplication): tables of straight-line
kernels, one per size, with two carry chains on x86-64. On x86-64 with
ADX, `flint_mpn_mulhigh_n` uses

- for $n \le 9$, assembly kernels computing $H''_n$, generated by
  `dev/gen_mulhigh_basecase.jl`;
- for $10 \le n \le 13$, peeling ([](#lem-split-high)) down to the
  9-limb kernel, with rows of `mpn_addmul_1`;
- for $14 \le n \le 22$, the splitting of [](#lem-split-high), with
  $p = n - 9$, $q = 9$ up to 18 limbs, so that one half is the 9-limb
  kernel, and a balanced split above;
- for $23 \le n \le 40$, a quadratic assembly loop
  (`_flint_mpn_mulhigh_basecase`, derived from GMP's
  `mpn_mullo_basecase` for Broadwell);
- for $40 < n \le 500$ (`FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF`), Mulders'
  recursion [](#alg-mulders-high) with $k$ from the tuned table, whose
  leaves are the methods above;
- above, the window $[n - 2, 2n)$ of a middle product by `fft_small`
  ([](#window)).

`flint_mpn_sqrhigh` has assembly kernels up to 8 limbs, assembly loops
(for even and odd $n$) up to 90 limbs, Mulders' recursion for squares up
to 580 limbs and `fft_small` above. `flint_mpn_mullow_n` has assembly
kernels up to 8 limbs, an assembly loop up to 50 limbs, then Mulders'
recursion [](#alg-mulders-low) (on operands padded to $n + 1$ limbs, so
that limb $n$ of the product is also obtained) up to 500 limbs and
`fft_small` above. On ARM64, high products and high squares have
assembly kernels up to 8 limbs and high products an assembly loop above;
the other cases, and all cases on other machines, use C code (C kernels
for high products up to 16 limbs, rows of `mpn_addmul_1` for the loops).
Without `fft_small`, Mulders' recursion is used up to 2000 limbs and a
full product above, so that the result is then exact.

In the measurements of [@AJ2025] on AMD Zen 3, the best strategy for the
high product was the 9-limb kernel with peeling up to 13 limbs, splitting
up to 40 limbs and Mulders' recursion above. At 9 limbs the high
product was about 1.8 times as fast as the full product. Compared with
MPFR's `mpfr_mulhigh_n`, it was 2 to 3 times as fast at 1, 2 and 9 limbs
and 1.1 times at 64 limbs; floating-point multiplication in `nfloat`,
built on `flint_mpn_mulhigh_normalised`, was 2.7 to 4.4 times as fast as
`mpfr_mul` at 1, 2 and 9 limbs and 1.2 times at 64 limbs.

## Use in FLINT {#use}

- `flint_mpn_mullow_n` and `flint_mpn_mulhigh_n`, `flint_mpn_sqrhigh`
  and their normalised versions, selected as in [](#selection); they are
  used by `nfloat` and `mp_real` arithmetic, by `mpn_mod`, and by the
  preinverted division and modular multiplication of `mpn_extras`.
- `flint_mpn_mulmid` for arbitrary windows, used by the Newton division
  and square root of `mp_real`, by `flint_mpn_sqrtrem` and by
  `flint_mpn_binv`.
- `_flint_mpn_mulhigh_known_low`, used by `flint_mpn_divisible` and the
  2-adic division and square root.
- The complex versions `flint_mpn_mulhigh_n_complex`,
  `flint_mpn_mullow_complex` and `flint_mpn_mulmid_complex`, and the
  radix-$b^e$ analogues `radix_mulmid` and `_radix_mulhigh_known_low`.
- For polynomials, the corresponding operations are the short and middle
  products of [](algorithms/polynomials/middle-product).

::: functions flint_mpn_mulhigh_n flint_mpn_mullow_n flint_mpn_mulmid _flint_mpn_mulhigh_known_low :::
