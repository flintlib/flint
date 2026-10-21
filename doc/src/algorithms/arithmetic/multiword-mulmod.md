---
title: Modular multiplication with a multiword modulus
summary: Reduces a product of two residues modulo an $\ell$-limb $n$ with a precomputed $\ell$-limb inverse, using one high and one low product and at most a few subtractions; also the two-limb Barrett reduction of `fmpz_mod`, Shoup's method and the matrix method for a fixed multiplier, and the fused $ab + cd$.
complexity: One $\ell \times \ell$-limb product, one high product and one low product of $\ell$ limbs, and $O(1)$ subtractions per modular product; a fixed multiplier saves the high product of the reduction (Shoup) or replaces both products by an $\ell \times \ell$ matrix-vector product and a two-limb reduction (matrix method).
order: 51
---

FLINT computes $ab \bmod n$ for residues $0 \le a, b < n$ when $n$ has
$\ell \ge 2$ limbs in three places: in `mpn_mod`
(`mpn_mod_mul`), in `mpn_extras` (`flint_mpn_mulmod_preinvn` and its
relatives, which `mpn_mod`, `fmpz` and the primality tests use) and in
`fmpz_mod` for two-limb moduli. All methods follow Barrett's idea
[@Bar1987] of replacing the division by $n$ with multiplications by a
precomputed approximation of $1/n$; for a single limb, the analogous
methods are described in [](algorithms/arithmetic/moller-granlund),
[](algorithms/arithmetic/barrett) and [](algorithms/arithmetic/shoup).
A textbook treatment of multiword Barrett and Montgomery reduction is
[@BZ2011, §2.4].

Throughout, $B = 2^{\bits}$, $n$ has exactly $\ell$ limbs
($B^{\ell - 1} \le n < B^{\ell}$), and $\nu$ is the number of leading zero
bits of its top limb.

## Normalisation {#normalisation}

The reduction works with the *normalised* modulus $D = n 2^{\nu}$, which
satisfies $B^{\ell}/2 \le D < B^{\ell}$. Reducing modulo $D$ is enough,
since for every integer $T$

$$
(2^{\nu} T) \bmod D = 2^{\nu} (T \bmod n).
$$

To reduce a product $ab$, one reduces $T = (2^{\nu} a) b$ modulo $D$ and
shifts the remainder right by $\nu$ bits. `mpn_mod` stores $D$ in its
context (`MPN_MOD_CTX_MODULUS_NORMED`) and shifts $a$ before the
multiplication; `flint_mpn_mulmod_preinvn` instead expects operands that
are already shifted.

## The precomputed inverse {#inverse}

`flint_mpn_preinvn` computes the $\ell$-limb number

$$
v = \left\lfloor \frac{B^{2\ell}}{D + 1} \right\rfloor - B^{\ell}.
$$

Since $B^{\ell}/2 < D + 1 \le B^{\ell}$, the quotient
$\mu = \lfloor B^{2\ell}/(D+1) \rfloor$ lies in $[B^{\ell}, 2B^{\ell})$,
so $0 \le v < B^{\ell}$ and $\mu = B^{\ell} + v$ is an approximation of
$B^{2\ell}/D$ from below whose leading 1 need not be stored. (For
$D = B^{\ell} - 1$, $v = 0$.) This is the $\ell$-limb analogue of the
one-limb inverse of Möller and Granlund [@GraMol2010]. It costs one
division of a $2\ell$-limb number, which `mpn_mod` performs once per
context.

## Reduction of a product {#preinvn}

::: algorithm Multiword reduction with a precomputed inverse {#alg-preinvn}
Input: $D$ normalised, $v$ as above, and $0 \le T < D B^{\ell}$.
Output: $T \bmod D$.

1. Write $T = T_1 B^{\ell} + T_0$ with $0 \le T_0 < B^{\ell}$.
2. Compute $q = T_1 + \lfloor T_1 v / B^{\ell} \rfloor$, the high half of
   $T_1 v$ plus $T_1$.
3. Compute $R = T - qD$ from the low $\ell + 1$ limbs of $T$ and of $qD$
   (a low product).
4. While $R \ge B^{\ell}$, set $R \leftarrow R - D$. Then, if $R \ge D$,
   set $R \leftarrow R - D$. Return $R$.
:::

The input bound holds for $T = (2^{\nu} a) b$ with $a, b < n$, since then
$T < D n < D B^{\ell}$. It implies $T_1 < D$.

::: theorem {#thm-preinvn}
Let $Q = \lfloor T / D \rfloor$. Then $Q - 5 \le q \le Q$. Consequently
$0 \le R = T - qD < 6D$, step 3 computes $R$ exactly, and step 4 performs
at most five subtractions and returns $T \bmod D$.
:::

::: proof
Since $\mu = \lfloor B^{2\ell}/(D+1) \rfloor$ and
$q = \lfloor T_1 \mu / B^{\ell} \rfloor$,

$$
q \le \frac{T_1 \mu}{B^{\ell}} \le \frac{T_1 B^{\ell}}{D + 1} < \frac{T_1 B^{\ell}}{D} \le \frac{T}{D},
$$

so $q \le Q$. In the other direction, $\mu > B^{2\ell}/(D+1) - 1$ and
$T_1 < B^{\ell}$ give

$$
q > \frac{T_1 \mu}{B^{\ell}} - 1 > \frac{T_1 B^{\ell}}{D + 1} - 2,
$$

while $Q \le T/D < (T_1 + 1) B^{\ell} / D$. Hence

$$
Q - q < \frac{T_1 B^{\ell}}{D(D+1)} + \frac{B^{\ell}}{D} + 2 < \frac{B^{\ell}}{D + 1} + \frac{B^{\ell}}{D} + 2 \le 6,
$$

using $T_1 < D$ and $D \ge B^{\ell}/2$. So $0 \le R < (Q - q + 1) D \le 6D$.
In particular $0 \le R < 6 B^{\ell} < B^{\ell + 1}$, so $R$ is determined
by $T - qD \bmod B^{\ell+1}$, which only involves the low $\ell + 1$
limbs of $T$ and of $qD$. Also $q \le Q < B^{\ell}$, so $q$ has $\ell$
limbs and the carry of the addition in step 2 is zero. Each subtraction
in step 4 keeps $R \equiv T \pmod D$ and $R \ge 0$; after the loop
$R < B^{\ell} \le 2D$, and the final subtraction gives $R < D$.
:::

The bound 5 is crude; it only matters that the number of corrections is
bounded by a small constant. The cost is one $\ell \times \ell$ product
(to form $T$), one high product $T_1 v$ and one low product $q D$ of
$\ell + 1$ limbs, which is about three products when all three are
computed in full. `flint_mpn_mul_or_mulhigh_n` and
`flint_mpn_mul_or_mullow_n` compute the high and low products with the
short product algorithms of `mpn_extras` when these are faster
([](algorithms/integers/short-products)). The approximate high product
is a lower bound that is short by at most one unit
(`flint_mpn_mulhigh_n`), so $q$ is still at most $Q$ and the loop in
step 4 runs at most once more.

For $\ell = 2$, `flint_mpn_mulmod_preinvn_2` (used by `mpn_mod_mul` and
the two-limb vector loops) is the same algorithm written with inline
two-limb arithmetic; the low product $qD$ is computed as three limbs.

### Variants {#variants}

- `flint_mpn_fmmamod_preinvn` computes $a_1 b_1 + a_2 b_2 \bmod n$. The
  shifted sum $T = 2^{\nu}(a_1 b_1 + a_2 b_2)$ is less than
  $2 D n < 2 D B^{\ell}$; one conditional subtraction of $D B^{\ell}$
  (from the high half) brings it below $D B^{\ell}$, and the algorithm
  above finishes. This costs two products and one reduction, against two
  products and two reductions for separate modular products.
- `flint_mpn_mod_preinvn` and `flint_mpn_divrem_preinvn` reduce numbers
  of any length $m \ge \ell$ by applying the algorithm to $\ell$ limbs of
  the dividend at a time, from the top.
- `mpn_mod_set_mpn`, which reduces unreduced sums in `mpn_mod`, uses
  `flint_mpn_mod_preinvn` for $\ell \ge 3$ and GMP's `mpn_divrem_2` for
  $\ell = 2$.

## Two-limb Barrett reduction in fmpz_mod {#two-limb}

For $B < n < B^2$, `fmpz_mod_ctx_init` precomputes the three-limb number
$m = \lfloor B^4 / n \rfloor$, and `fmpz_mod_mul` (through `_fmpz_mod_mul2`)
reduces $x = ab < n^2$ by classical Barrett reduction:

$$
q = \left\lfloor \frac{\lfloor x / B \rfloor \, m}{B^3} \right\rfloor, \qquad r = x - qn.
$$

::: proposition {#prop-two-limb}
$0 \le r < 3n$, so at most two subtractions of $n$ give $x \bmod n$.
:::

::: proof
Let $x_1 = \lfloor x / B \rfloor$ and $Q = \lfloor x/n \rfloor$. Then
$q \le x_1 m / B^3 \le x_1 B / n \le x / n$, so $q \le Q$. Also
$m > B^4/n - 1$, $x_1 B > x - B$ and $x_1 < B^3$ (since $x < B^4$), so

$$
\frac{x_1 m}{B^3} > \frac{x_1 B}{n} - \frac{x_1}{B^3} > \frac{x - B}{n} - 1 > \frac{x}{n} - 2,
$$

using $n > B$. Hence $q > x/n - 3$, $q \ge Q - 2$ and $r < 3n$.
:::

The code computes the five-limb product $x_1 m$ with ten word products,
skipping only the low half of the lowest one (which cannot carry into
the limbs that are used), and the low three limbs of $qn$. The remainder
fits in three limbs, so the corrections are three-limb comparisons. The
modulus $n = B$ has its own code (`_fmpz_mod_mul2s`, a wrapping
product), and one-limb moduli use `nmod_mul`.

## A fixed multiplier {#fixed}

When the same $a$ multiplies many $b$, as in scaling a vector or in the
division steps of Euclid's algorithm, precomputation for $a$ pays off.
`flint_mpn_mulmod_want_precond` chooses the method from $\ell$ and the
number of products: none for fewer than 4 products, Shoup's method for
$\ell \le 10$ (and for $\ell \le 12$ with at most 12 products), the
matrix method for larger $\ell$ up to 64 limbs (and beyond with enough
products), so that `mpn_mod` uses Shoup's method up to 10 limbs and the
matrix method for 11 to 16 limbs. For $\ell = 2$ with the top bit of $n$
set, no precomputation is used through this function, but the vector
functions of `mpn_mod` have their own inlined two-limb Shoup code.

### Shoup's method {#shoup}

Shoup's method [@Sho2021; @Har2014] precomputes the $\ell$-limb quotient

$$
a' = \lfloor a B^{\ell} / n \rfloor
$$

(`flint_mpn_mulmod_precond_shoup_precompute`, one division). A product
is then computed as

$$
q = \lfloor a' b / B^{\ell} \rfloor, \qquad r = ab - qn,
$$

where $r$ is obtained from the low halves of $ab$ and $qn$ alone.

::: proposition {#prop-shoup}
$0 \le r < n (1 + b / B^{\ell}) < 2n$.
:::

::: proof
Write $a B^{\ell} = a' n + e$ with $0 \le e < n$ and $a' b = q B^{\ell} + f$
with $0 \le f < B^{\ell}$. Then

$$
ab B^{\ell} = a' b n + e b = (q B^{\ell} + f) n + e b,
$$

so $r = ab - qn = (fn + eb)/B^{\ell}$, which is nonnegative and less than
$(B^{\ell} n + nb)/B^{\ell}$.
:::

The cost per product is one high product and two low products of
$\ell$ limbs, against one full, one high and one low product without
precomputation. If $n < B^{\ell}/2$ ($\nu > 0$), then $r < 2n < B^{\ell}$
fits in $\ell$ limbs and one conditional subtraction finishes.

The implementation computes $q$ with the approximate high product
`flint_mpn_mulhigh_n`, a lower bound which may be $q - 1$. This can only
happen when $f < c B^{\ell - 1}$ for a small constant $c$ (the error of
the approximation, a few units of the limb below the high half). Then the
computed remainder is
$r + n < n (1 + c/B + b/B^{\ell}) < n (3/2 + c/B) < 2n$ when
$b < n < B^{\ell}/2$, so for $\nu > 0$ one subtraction still suffices.
If the top bit of $n$ is set, $r$ may exceed
$B^{\ell}$; the code then computes modulo $B^{\ell + 1}$ and subtracts $n$
up to twice (for $\ell = 2$ it uses an exact high product and one
subtraction).

### The matrix method {#matrix}

The matrix method precomputes the $\ell$ residues

$$
a_i = (2^{\nu} a B^{i}) \bmod D, \qquad 0 \le i < \ell
$$

(`flint_mpn_mulmod_precond_matrix_precompute`, $\ell - 1$ one-limb
shifts and reductions). For $b = \sum_i b_i B^i$,

$$
2^{\nu} a b \equiv \sum_{i=0}^{\ell - 1} b_i a_i \pmod D,
$$

and the right-hand side is a matrix-vector product of $\ell^2$ word
products, less than $\ell B D < B^{\ell + 2}$. It is reduced modulo $D$
with two steps of one-limb quotient estimation (`flint_mpn_mod_preinv1`)
and shifted right by $\nu$ bits. Compared with the full product followed
by the reduction of [](#alg-preinvn), the high and low products of the
reduction are replaced by $O(\ell)$ work, at the price of $\ell^2$ limbs
of precomputed data. The fused version `flint_mpn_fmmamod_precond_matrix`
adds two such matrix-vector products before reducing.

## Inversion and powering {#inv-pow}

Inverses modulo a multiword $n$ are computed with GMP's extended
Euclidean algorithm (`mpn_gcdext` in `mpn_mod_inv`, `fmpz_gcdinv` in
`fmpz_mod_inv`), which is subquadratic for large operands. Powers with
large exponents use GMP's `mpz_powm` (in `mpn_mod_pow_fmpz`, and through
`fmpz_powm` in `fmpz_mod_pow_fmpz` for moduli below 70000 bits). It uses
Montgomery's representation [@Mon1985], in which no quotient estimate is
needed, with sliding windows; `mpn_mod` uses it for exponents with more
than 6 bits. For moduli of 70000 bits and more (and when
`fft_small` is available), `fmpz_mod_pow_fmpz` switches to generic
sliding-window exponentiation with `fmpz_mod_mul`, whose products use
the FFT.

## Use in FLINT {#use}

- `mpn_mod_mul`, `mpn_mod_sqr`, `_mpn_mod_vec_mul`: [](#alg-preinvn),
  for $\ell = 2$ through `flint_mpn_mulmod_preinvn_2`.
- `mpn_mod_fmma` and the length-2 dot products:
  `flint_mpn_fmmamod_preinvn`.
- `_mpn_mod_vec_mul_scalar`, `_mpn_mod_vec_addmul_scalar`,
  `_mpn_mod_vec_submul_scalar`: [Shoup's method](#shoup) or the
  [matrix method](#matrix) for vectors of at least 4 (or 6) elements.
- `_mpn_mod_poly_divrem_q1_preinv1_fmma_precond` and
  `_mpn_mod_poly_divrem_q1_preinv1_karatsuba_precond`: the matrix method
  for the two quotient coefficients of a Euclidean step.
- `fmpz_mod_mul` for $B < n < B^2$: [](#prop-two-limb).
- `flint_mpn_mulmod_preinvn`, `flint_mpn_powmod_preinvn` and
  `flint_mpn_sqrtmod_preinv` in `mpn_extras`, used by `mpn_mod_sqrt` and
  the integer square root modulo a prime.
- Sums of many products are not reduced one product at a time; see
  [](algorithms/arithmetic/multiword-delayed-reduction).

::: functions mpn_mod_mul flint_mpn_mulmod_preinvn flint_mpn_preinvn flint_mpn_mulmod_precond_shoup flint_mpn_mulmod_precond_matrix :::
