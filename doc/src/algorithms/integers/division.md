---
title: Integer division and Newton iteration
summary: Euclidean division of multiprecision integers by schoolbook division with a 3-by-2 inverse, Burnikel–Ziegler divide and conquer, and Newton iteration with the Karp–Markstein trick; exact division by 2-adic (Hensel) methods; division by a fixed divisor with a precomputed inverse; and square roots by the same techniques.
complexity: Schoolbook $O(mn)$; divide and conquer $O(M(n) \log n)$; Newton division, exact division and square root $O(M(n))$, for an $n$-limb divisor and quotient.
order: 20
---

Given integers $a \ge 0$ of $m$ limbs and $b > 0$ of $n$ limbs, Euclidean
division computes the quotient $q = \lfloor a/b \rfloor$ of $m - n + 1$
limbs and the remainder $r = a - qb \in [0, b)$. Most of the code is in
`mpn_extras`; it replaces GMP's `mpn_tdiv_qr` throughout FLINT, and the
`fmpz` division functions ([](fmpz/division)) reach it through an
`mpz`-like layer (`flint_mpz_tdiv_qr` and its relatives) that applies the
rounding conventions at the limb level.

Throughout, $B = 2^{\bits}$, $M(n)$ is the cost of an $n \times n$-limb
multiplication ([](algorithms/integers/multiplication)), and a divisor is
*normalised* if its top limb is at least $B/2$, which is achieved by
shifting both operands left by the same number of bits.

## Schoolbook division {#schoolbook}

The schoolbook method (Knuth's Algorithm D [@Knu1997, §4.3.1]) produces
one quotient limb per step: it estimates the next limb from the top limbs
of the current remainder and the divisor, subtracts the product of the
estimate and the divisor, and corrects the estimate if the result is
negative. The estimate is good because the divisor is normalised.

::: lemma Quotient limb estimate {#lem-estimate}
Let $v = \sum_{i<n} v_i B^i$ with $v_{n-1} \ge B/2$, let $0 \le u < vB$
with $u = \sum_{i \le n} u_i B^i$, and let $q = \lfloor u/v \rfloor$
(so $q < B$). Then
$$
\hat q = \min\Bigl(\Bigl\lfloor \frac{u_n B + u_{n-1}}{v_{n-1}} \Bigr\rfloor, B - 1\Bigr)
$$
satisfies $q \le \hat q \le q + 2$.
:::

::: proof
Write $t = u_n B + u_{n-1} = \lfloor u / B^{n-1} \rfloor$. Since
$q v_{n-1} B^{n-1} \le qv \le u < (t + 1) B^{n-1}$, we have $q v_{n-1}
\le t$, so $q \le \lfloor t / v_{n-1} \rfloor$; also $q \le B - 1$, hence
$q \le \hat q$. For $n = 1$ the estimate is exact. Otherwise $v >
B^{n-1}$, and $\hat q \le t / v_{n-1} \le u / (v_{n-1} B^{n-1}) < u / (v -
B^{n-1})$, using $v < (v_{n-1} + 1) B^{n-1}$. Suppose $\hat q \ge q +
3$. Since $q > u/v - 1$,
$$
3 \le \hat q - q < \frac{u}{v - B^{n-1}} - \frac{u}{v} + 1
= \frac{u B^{n-1}}{v (v - B^{n-1})} + 1,
$$
so $u/v > 2 (v - B^{n-1}) / B^{n-1} \ge 2 (v_{n-1} - 1)$. Then $B - 4 \ge
\hat q - 3 \ge q \ge 2 (v_{n-1} - 1)$, so $v_{n-1} < B/2$, a
contradiction.
:::

Using the top two limbs of the divisor and three of the remainder (a
"3-by-2" division) sharpens the estimate to $q \le \hat q \le q + 1$, and
with a precomputed reciprocal of the top two divisor limbs each estimate
costs a few multiplications instead of a hardware division, as in
Möller and Granlund [@GraMol2010] ([](algorithms/arithmetic/moller-granlund)).
The method costs $O((m - n + 1) n)$ limb operations.

In FLINT, divisors of up to `FLINT_MPN_DIV_SMALL_BN` $= 7$ limbs with
short quotients are divided by dedicated functions, one per divisor
length, which keep the normalised divisor and the partial remainder in
registers and normalise the dividend on the fly; one-limb divisors use
GMP's assembly `mpn_divrem_1`. Larger divisors use a port of GMP's
`mpn_sbpi1_div_qr` (schoolbook with a 3-by-2 inverse) with FLINT's
multiplication.

## Divide and conquer {#divconquer}

Burnikel and Ziegler [@BurZie1998] observed that the schoolbook method
works with blocks of limbs in place of limbs: dividing a $2k$-block
number by a $k$-block number takes two steps, each dividing $3k/2$ by
$k$ blocks, and each of those is a recursive division of the top half
followed by one multiplication of half-size blocks.

::: proposition {#prop-dc-cost}
Divide-and-conquer division of a $2n$-limb by an $n$-limb integer costs
$D(n) = 2 D(n/2) + 2 M(n/2) + O(n)$, so $D(n) = O(M(n) \log n)$, and
$D(n) = O(M(n))$ if $M(n) = n^{\alpha}$ with $\alpha > 1$.
:::

::: proof
Each half step does one recursive division of size $n/2$ and one
$n/2 \times n/2$ product (plus an $O(1)$ number of corrections, each
$O(n)$), which gives the recurrence. Unrolling, $D(n) = \sum_{i}
2^{i+1} M(n/2^{i+1}) + O(n \log n)$. Since $M(n)/n$ is nondecreasing,
$2^{i} M(n/2^{i}) \le M(n)$, which gives $O(M(n) \log n)$; if $M(n) =
n^\alpha$, the terms decrease geometrically with ratio $2^{1-\alpha}$ and
the sum is $O(M(n))$.
:::

FLINT uses a port of GMP's `mpn_dcpi1_div_qr` with its own
multiplication from `FLINT_MPN_DIV_DC_CUTOFF` $= 18$ limbs (on recent
x86-64), together with variants that compute only the quotient, or an
approximate quotient (at most one too large), which save the
corrections of the remainder. When only a quotient about as long as the
divisor is needed, a short division in the style of Mulders, which forms
only the high halves of the products, is faster still [@HZ2011]
([](algorithms/integers/short-products)).

## Newton iteration {#newton}

For large operands, division reduces to multiplication by computing an
approximate reciprocal by Newton's method. Viewing $b$ as a fixed-point
number in $[1/2, 1)$, the iteration
$$
x_{k+1} = x_k + x_k (1 - b x_k)
$$
converges quadratically to $1/b$.

::: lemma {#lem-newton}
If $x = (1 - \epsilon)/b$, then $x' = x + x(1 - bx)$ satisfies $x' = (1 -
\epsilon^2)/b$.
:::

::: proof
$1 - bx = \epsilon$, so $b x' = (1 - \epsilon)(1 + \epsilon) = 1 -
\epsilon^2$.
:::

Each step doubles the number of correct limbs, and only the last step
needs full precision. Since $1 - bx$ is small, its top half is known to
be zero (or all ones), and the product $b x$ can be replaced by a middle
product of half the cost. Summing over the precisions $n, n/2, n/4,
\ldots$, the reciprocal costs a small constant times $M(n)$.

Karp and Markstein [@KM1997] noticed that the last step can produce the
quotient directly, so that the reciprocal is only needed to half
precision.

::: theorem Karp–Markstein division {#thm-karp-markstein}
Let $x = (1 - \epsilon)/b$, $y = ax$ and $q = y + x(a - by)$. Then $q =
(a/b)(1 - \epsilon^2)$.
:::

::: proof
$a - by = a - abx = a(1 - bx) = a \epsilon$, so $q = ax + x a \epsilon =
a x (1 + \epsilon) = (a/b)(1 - \epsilon)(1 + \epsilon)$.
:::

With $x$ accurate to $n/2$ limbs, $y$ truncated to $n/2$ limbs, and the
correction $x(a - by)$ needing only $n/2$ limbs, the quotient to $n$
limbs costs a half-precision reciprocal and a few half-size products. In
practice the computed quotient carries one or two guard limbs, and its
error, from the truncations, is bounded by a few units of the last guard
limb.

FLINT implements this in fixed-point arithmetic (`_mp_real_inv_newton`
and `_mp_real_div_newton` of `mp_real`, with some steps of third order,
which triple the precision instead of doubling it). `_flint_mpn_tdiv_qr`
calls it when the divisor and the quotient both have at least
`FLINT_MPN_TDIV_QR_NEWTON_CUTOFF` $= 950$ limbs, or the divisor has at
least 608 limbs and the quotient at least twice as many. It computes the
quotient with two guard limbs; the error is at most 4 units of the first
guard limb, so if that limb lies in $[2, B - 2]$ the integer part is
certain and the remainder follows from a low product. Otherwise the
quotient is corrected by $O(1)$ steps using a full product.

## Precomputed inverses and long quotients {#preinv}

When the quotient is much longer than the divisor, or many dividends are
divided by the same divisor, the reciprocal is computed once and used for
every block of quotient limbs, as in Barrett's method
([](algorithms/arithmetic/barrett)):

- `flint_mpn_preinvn` computes an $n$-limb approximation of the
  reciprocal of a normalised $n$-limb divisor, and
  `flint_mpn_divrem_preinvn` divides by it $n$ limbs at a time: each
  block of the quotient is the high part of the product of the inverse
  with the top limbs of the remainder, followed by a low product and at
  most a few corrections. `fmpz_fdiv_qr_preinvn` exposes it, and
  `fmpz_mod_mul` uses it for moduli of more than 20 limbs.
- For a dividend of at least three times the length of a divisor of
  hundreds of limbs, `_flint_mpn_tdiv_qr` divides from the top down in
  blocks of $2n$ by $n$ limbs with one approximate inverse of $n + 3$
  fraction limbs shared by all blocks
  (`_flint_mpn_tdiv_qr_unbalanced`). For divisors of a few dozen limbs
  and long dividends it uses the `preinvn` division instead, which is up
  to 40% faster than GMP there.

## Exact division {#hensel}

If $b$ is known to divide $a$, the quotient can be computed from the low
limbs instead of the high ones, by 2-adic arithmetic. This is Hensel's
idea applied to integers, introduced for division by Jebelean
[@Jeb1993].

::: proposition {#prop-hensel}
Let $b$ be odd and $0 \le q < B^k$ with $a = qb$. Then $q \equiv a
b^{-1} \pmod{B^k}$, so $q$ is the residue in $[0, B^k)$ of $a b^{-1}$.
:::

::: proof
$b$ is invertible modulo $B^k$ because it is odd, and $qb \equiv a$
determines $q$ modulo $B^k$; since $0 \le q < B^k$, it is that residue.
:::

A power of two and whole zero limbs in $b$ are removed first (they must
divide $a$ equally). The quotient of an $m$-limb by an $n$-limb integer
has $k = m - n + 1$ limbs (one fewer if the top limb of $a$ is below that
of $b$), and by the proposition it depends only on the low $k$ limbs of
$a$ and $b$, so long operands are truncated before dividing. The
quotient is developed from the least significant limb upwards: $q_0 = a_0
\cdot (b_0^{-1} \bmod B)$, then $b q_0$ is subtracted and the process
continues. Compared with Euclidean division, no quotient limb is ever
estimated or corrected, and the unused high limbs are never touched.

The 2-adic inverse needed for large operands is computed by Newton
iteration, which converges as in [](#lem-newton):

::: lemma {#lem-binv}
If $by \equiv 1 + B^k h \pmod{B^{2k}}$, then $y' = y(2 - by)$ satisfies
$b y' \equiv 1 \pmod{B^{2k}}$.
:::

::: proof
$b y' = by (2 - by) = (1 + B^k h)(1 - B^k h) = 1 - B^{2k} h^2$.
:::

The analogue of [](#thm-karp-markstein) also holds: with $y = b^{-1}
\bmod B^{\lceil k/2 \rceil}$, the low half of the quotient is $q_0 = a y
\bmod B^{\lceil k/2 \rceil}$, and the high half is $y (a - b q_0) /
B^{\lceil k/2 \rceil}$ modulo the remaining precision.

`_flint_mpn_divexact` (behind `flint_mpn_divexact` and `fmpz_divexact`)
uses, for a truncated divisor of $d = \min(n, k)$ limbs: register-based
code for $d \le 7$; ports of GMP's schoolbook and divide-and-conquer
Hensel division (`mpn_sbpi1_bdiv_q`, `mpn_dcpi1_bdiv_q`, from 40 and 80
limbs) with FLINT's multiplication; and for $d \ge$
`FLINT_MPN_DIVEXACT_NEWTON_CUTOFF` $= 1187$ limbs, or a quotient of at
least four times the length of a divisor of 531 limbs or more,
`flint_mpn_bdiv_q`, which computes the inverse with `flint_mpn_binv`
(Newton iteration with middle products) and uses Karp–Markstein Hensel
division, or block division for short divisors.

The same tools give:

- `flint_mpn_divisible` (behind `fmpz_divisible`): the 2-adic part, then
  a cheap trial-division screen modulo small primes for large operands,
  then the Hensel remainder, which is zero exactly when $b$ divides $a$;
- `flint_mpn_divexact_preinv`, which stores the 2-adic inverse of a fixed
  divisor, used by `_fmpz_vec_scalar_divexact_fmpz` to divide a whole
  vector by one large integer;
- `fmpz_invmod_2exp`, `fmpz_divmod_2exp`, `fmpz_sqrtmod_2exp` and
  `fmpz_rsqrtmod_2exp`, which expose inversion, division and (inverse)
  square roots modulo $2^N$.

## Square roots {#sqrt}

The integer square root $s = \lfloor \sqrt a \rfloor$ with remainder $r =
a - s^2 \in [0, 2s]$ is computed by the same techniques.

**Divide and conquer.** Zimmermann's "Karatsuba square root" [@Zim1999]
splits a normalised $a = a_3 \beta^3 + a_2 \beta^2 + a_1 \beta + a_0$
into quarters and computes
$$
(s', r') = \operatorname{sqrtrem}(a_3 \beta + a_2), \quad
(q, u) = \operatorname{divrem}(r' \beta + a_1, 2 s'),
$$
$$
s = s' \beta + q, \quad r = u \beta + a_0 - q^2,
$$
followed by one correction ($r \gets r + 2s - 1$, $s \gets s - 1$) if
$r < 0$. Zimmermann proves that the result is correct; the cost is one
recursive square root, one division and one squaring of half size, so
the method costs $O(M(n))$ when $M(n) = n^\alpha$ with $\alpha > 1$ and
$O(M(n) \log n)$ in general.

**Newton iteration.** The iteration $y_{k+1} = y_k + y_k (1 - a y_k^2)/2$
converges quadratically to $1/\sqrt a$: if $y = (1 - \epsilon)/\sqrt a$,
then $\sqrt a\, y_{k+1} = (1 - \epsilon)(1 + \epsilon - \epsilon^2/2) = 1
- \tfrac32 \epsilon^2 + \tfrac12 \epsilon^3$. As for division, the last
step computes the root directly (Karp–Markstein): with $s = ay$, the
value $s + y (a - s^2)/2$ is accurate to twice the precision of $y$.

`_flint_mpn_sqrtrem` uses dedicated code for one to four limbs (a
floating-point square root, a Newton correction, and one divide-and-
conquer step for three and four limbs), FLINT's divide-and-conquer code
(as GMP's `mpn_sqrtrem`, but with FLINT's division and multiplication) up
to `FLINT_MPN_SQRTREM_NEWTON_CUTOFF` $= 7423$ limbs, and the
Newton–Karp–Markstein method (`_mp_real_sqrt_newton`, certified as for
division by inspecting a guard limb) above. `fmpz_sqrt`, `fmpz_sqrtrem`,
`fmpz_is_square` and `fmpz_perfect_sqrt` use it ([](fmpz/roots)).

## Use in FLINT {#use}

`_flint_mpn_tdiv_qr(q, r, a, m, b, n)` (behind `flint_mpn_tdiv_qr`,
[](mpn_extras/tdiv_qr)) chooses, for $k = m - n + 1$ quotient limbs,
with the cutoffs of recent x86-64 processors:

1. $n \le 3$, or $n \le 7$ with $m \le 3n$: register-based schoolbook
   division (`tdiv_qr_small`).
2. $n \ge 950$ and $k \ge 950$, or $n \ge 608$ and $k \ge 2n$: Newton
   (Karp–Markstein) division, or for $m \ge 3n$ block division with a
   shared inverse.
3. $n \ge 64$ and $m \ge 4n$, or $n \ge 512$ and $m \ge 3n$: block
   division with a shared inverse.
4. $n \ge 32$ and $m \ge 4n$, or $n \ge 4$ and $m \ge 32n$: division with
   `flint_mpn_divrem_preinvn`.
5. Otherwise: the ports of GMP's schoolbook and divide-and-conquer
   division, with FLINT's multiplication.

Functions that need only the quotient (`flint_mpn_tdiv_q`), only an
approximate quotient (`flint_mpn_divapprox`, used by `mp_real`
division), or a remainder by a precomputed inverse use variants of the
same algorithms that skip the work they do not need.

::: functions fmpz_tdiv_qr fmpz_divexact fmpz_divisible fmpz_fdiv_qr_preinvn fmpz_sqrtrem :::
