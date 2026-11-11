---
title: Barrett and Lemire reduction
summary: Reduce a single word modulo $n$ with two products and at most one correction, using $\lfloor 2^{\bits}/n \rfloor$ (Barrett) or $\lceil 2^{\bits}/n \rceil$ (Lemire, Kaser and Kurz).
complexity: Barrett: one high and one low product and one conditional subtraction. Lemire: one low and one high product, no correction.
order: 40
---

The methods of Barrett and of Lemire, Kaser and Kurz reduce a *single*
word $x$ modulo $n$, rather than a two-word product, using a precomputed
approximation of $B/n$. They are cheaper than a two-word
[Möller–Granlund](algorithms/arithmetic/moller-granlund) reduction when
the value to be reduced is known to fit in one word, as in hashing,
sieving and table lookups. Barrett's method also reduces two and three
words for $n < B/2$ with a two-word constant ([](#two-words)), and two-word
values below $2^{\nu + \bits - 2}$, where $\nu$ is the bit length of $n$,
with the one-word constant of `NMOD_RED` ([](#bounded)). Throughout,
$B = 2^{\bits}$.

## Barrett reduction {#barrett}

Barrett's method [@Bar1987] estimates the quotient as follows.

::: proposition {#prop-barrett}
Let $1 \le n < B$, let $m$ be an integer with $B/n - 1 \le m \le B/n$
and $0 \le x < B$, and let $\tilde q = \lfloor x m / B \rfloor$. Then
$\lfloor x/n \rfloor - 1 \le \tilde q \le \lfloor x/n \rfloor$, so
$0 \le x - \tilde q n < 2n$.
:::

::: proof
From $B/n - 1 \le m \le B/n$ we get
$$
\frac{x}{n} - 1 < \frac{x}{n} - \frac{x}{B} \le \frac{x m}{B} \le \frac{x}{n},
$$
using $x < B$. Taking floors gives the claim.
:::

The value $x - \tilde q n$ is at most $x$, so it is computed exactly with
word arithmetic even when $2n > B$. `n_mod_barrett_lazy` returns it, and
`n_mod_barrett` subtracts $n$ once if needed. `n_barrett_precomp` computes
$m = \lfloor B/n \rfloor$ with one division, for $n \ge 2$;
`n_barrett_precomp_preinv` computes $m = \lfloor (B - 1)/n \rfloor$
without a division ([below](#nmod)), for $n \ge 1$. Both satisfy the
hypothesis: $\lfloor (B - 1)/n \rfloor \ge ((B - 1) - (n - 1))/n = B/n - 1$,
since the remainder of $B - 1$ by $n$ is at most $n - 1$. They differ only
when $n$ is a power of two, where the second one is $B/n - 1$.

## Barrett reduction with an `nmod_t` {#nmod}

`NMOD_RED` reduces a word with the same kind of quotient estimate, but its
constant comes from the `nmod_t`, so it needs no division and no extra
field. Let $\nu$ be the number of bits of $n$, that is
$2^{\nu - 1} \le n < 2^{\nu}$ and $\mathit{norm} = \bits - \nu$. The
`ninv` field is the [Möller–Granlund reciprocal](algorithms/arithmetic/moller-granlund#reciprocal)
of $n 2^{\mathit{norm}}$.

::: lemma {#lem-nmod-m}
Let $m = 2^{\bits - 1} + \lfloor \mathit{ninv}/2 \rfloor$. Then
$m = \lfloor (2^{\nu + \bits - 1} - 1)/n \rfloor$, and
$2^{\bits - 1} \le m < B$.
:::

::: proof
By definition $B + \mathit{ninv} = \lfloor (B^2 - 1)/(n 2^{\mathit{norm}}) \rfloor
= \lfloor (2^{\nu + \bits} - 1)/n \rfloor$, using
$\lfloor \lfloor x \rfloor / k \rfloor = \lfloor x/k \rfloor$ for integers $k \ge 1$ and
$\lfloor (B^2 - 1)/2^{\mathit{norm}} \rfloor = 2^{\nu + \bits} - 1$. Hence
$m = \lfloor (B + \mathit{ninv})/2 \rfloor = \lfloor (2^{\nu + \bits} - 1)/(2n) \rfloor$.
With $y = 2^{\nu + \bits - 1}$, an integer $k$ satisfies $2nk \le 2y - 1$ if
and only if $2nk \le 2y - 2$, since $2nk$ is even; so this is
$\lfloor (y - 1)/n \rfloor$. The bounds follow from $\mathit{ninv} < B$.
:::

::: theorem {#thm-nmod-red}
Let $1 \le n < B$, $0 \le a < B$, $m$ as in [](#lem-nmod-m), and
$\tilde q = \lfloor a m / 2^{\nu + \bits - 1} \rfloor$. Then
$\lfloor a/n \rfloor - 1 \le \tilde q \le \lfloor a/n \rfloor$, so
$0 \le a - \tilde q n < 2n$.
:::

::: proof
Write $mn = 2^{\nu + \bits - 1} - 1 - \rho$ with $0 \le \rho < n$
([](#lem-nmod-m)). Then
$$
\frac{a m}{2^{\nu + \bits - 1}} = \frac{a}{n} - \frac{a (1 + \rho)}{n\, 2^{\nu + \bits - 1}},
$$
and the subtracted term is in $[0, 1)$ because $a < B$, $1 + \rho \le n$
and $\nu \ge 1$. Taking floors gives the claim.
:::

So `NMOD_RED` computes $\tilde q$ as the high word of $a m$ shifted right
by $\bits - 1 - \mathit{norm}$ bits, subtracts $\tilde q n$ (at most $a$, so
there is no overflow) and corrects once. This holds for every modulus
$1 \le n < B$, including powers of two and $n \ge 2^{\bits - 1}$, and every
word $a$. The cost is one high product, one low product and one comparison,
against a full product, a low product and two corrections for a
[Möller–Granlund](algorithms/arithmetic/moller-granlund) step on the
two-word number $(0, a)$.

When $n$ has `FLINT_BITS` bits (`mod.norm` $= 0$), $a < 2n$ for every word $a$,
so `NMOD_RED` skips the product and only corrects once. The test on `mod.norm`
is moved out of loops by the compiler, and the remaining loop of comparisons
can then be trivially vectorized.

`n_mod2_preinv`, `n_divrem2_preinv`, `n_div2_preinv` and `n_divrem_preinv`
use the same step with $\mathit{norm}$ computed from $n$, and the last
three also return $\tilde q$, corrected. Before the correction, the proof
of [](#thm-nmod-red) gives
$a - \tilde q n < n + a (1 + \rho)/2^{\nu + \bits - 1} < n + 2$, since
$a < B$ and $1 + \rho \le n < 2^{\nu}$; `n_mod2_preinv_lazy` returns this
value.

Shifting $m$ gives the constants of the same family for other input sizes:
$\lfloor m / 2^k \rfloor = \lfloor (2^{\nu + \bits - 1 - k} - 1)/n \rfloor$
for $0 \le k < \nu + \bits - 1$. Indeed, with $e = \nu + \bits - 1$,
$\lfloor m/2^k \rfloor = \lfloor (2^e - 1)/(2^k n) \rfloor
= \lfloor \lfloor (2^e - 1)/2^k \rfloor / n \rfloor$ by the nested floor
identity used in [](#lem-nmod-m), and $\lfloor (2^e - 1)/2^k \rfloor = 2^{e - k} - 1$.
For instance $k = \nu - 1$ gives $\lfloor (B - 1)/n \rfloor$, which
`n_barrett_precomp_preinv` computes this way, without a division, and on
64-bit machines $k = 32$ gives the 32-bit constant
$\lfloor m/2^{32} \rfloor = 2^{31} + \lfloor \mathit{ninv}/2^{33} \rfloor$ of the
SIMD kernels of `_nmod_vec_reduce` ([](nmod_vec/simd#reduce)).

## Two and three words {#two-words}

For $n < B/2$, the same estimate reduces a two-word $x < B^2$ with the
two-word constant $Q = \lfloor B^2/n \rfloor$: [](#prop-barrett) with
$B^2$ in place of $B$ shows that $\tilde q = \lfloor x Q / B^2 \rfloor$
satisfies $0 \le x - \tilde q n < 2n < B$. So the remainder is computed
from the low words of $x$ and $\tilde q n$ alone, and one subtraction of
$n$ finishes. `n_ll_mod_barrett` takes the low word of $\tilde q$, the
third word of $Q x$, from the low three words of this product
(`FLINT_MPN_MUL_3P2X2`, three word products).

A three-word $y = y_2 B^2 + y_1 B + y_0$ is first folded into two words,
$x = y_2 b_2 + y_1 b_1 + y_0 \equiv y \pmod n$ with $b_1 = B \bmod n$
and $b_2 = B^2 \bmod n$. Since $b_1, b_2 \le n - 1$,
$x \le (B - 1)(2n - 1) < B^2$ for $n < B/2$, whatever the words $y_i$
(`n_lll_mod_barrett`).

The three constants are the quotient and the remainders of the division of
$2^{\mathit{norm}} B^2$ by the normalized $d = 2^{\mathit{norm}} n$, in two
steps of [](algorithms/arithmetic/moller-granlund#alg-div) with the inverse
`ninv`: the first gives the high word $\lfloor B/n \rfloor$ of $Q$ and the
remainder $2^{\mathit{norm}} b_1$, the second the low word of $Q$ and
$2^{\mathit{norm}} b_2$ (`n_ll_barrett_precomp_preinv`). This also holds
for $n \ge B/2$, for which only $b_1$ and $b_2$ are used (next section).

## Three words, full-word moduli {#fold}

For $n \ge B/2$ the estimate with $Q$ no longer leaves a remainder that
fits a word, and FLINT keeps the fold above, followed by one
[Möller–Granlund](algorithms/arithmetic/moller-granlund) step on the two
words of $x$ (`NMOD_RED3_FOLD_FULLWORD`). Here $b_1 \le B - n$ and
$b_2 \le n - 1$, so $b_1 + b_2 + 1 \le B$ and
$x \le (B - 1)(b_1 + b_2 + 1) \le (B - 1) B$: $x$ fits two words, and its
high word, below $B < 2n$, is brought below $n$ by one conditional
subtraction, as the Möller–Granlund step requires. This is cheaper than
reducing $y_2$ and then taking two Möller–Granlund steps (`NMOD_RED3`).

The subtraction can be skipped when the top words $y_2$ are known to be at
most some $Y$, typically for sums of products of residues. The high word
of $x$ is below $n$ for every such input exactly when
$Y b_2 + (B - 1)(b_1 + 1) < nB$, since the left-hand side is the largest
fold, reached at $y_2 = Y$ and $y_1 = y_0 = B - 1$
(`_nmod_red3_fold_can_skip_sub`, `NMOD_RED3_FOLD_FULLWORD_NOSUB`). It
holds for every $n \ge B/2$ as soon as $2Y(Y + 1) \le B$, which covers the
sums of up to about $2^{\bits/2}$ products. Write $n = B/2 + \delta$; the
condition reads $Y b_2 < 2\delta B - n + 1$. If $2\delta \ge Y + 1$, it
follows from $b_2 \le n - 1$. Otherwise $\delta \le Y/2 < 2^{\bits/2 - 1}$,
so $b_1 = B/2 - \delta$, $b_2 = 4\delta^2$ (or $b_1 = b_2 = 0$ when
$\delta = 0$), and $Y \le B/(2(Y + 1)) < B/(4\delta)$ gives
$Y b_2 \le \delta B < 2\delta B - n + 1$. For larger $Y$ it can fail, for
$n$ just above $B/2$.
`_nmod_red3_precomp_init` evaluates the test once per modulus and bound,
and `_nmod_red3_precomp` then chooses between `n_lll_mod_barrett` and the
two folds.

## Two words below a bound {#bounded}

The constant $m$ of [](#nmod) reduces two-word integers below
$2^{\nu + \bits - 2}$, which covers the products of two residues when
$\nu \le \bits - 2$. The estimate uses the top $\bits$ bits of the input.

::: theorem {#thm-bounded}
Let $3 \le \nu \le \bits - 1$, $m$ as in [](#lem-nmod-m),
$0 \le a < 2^{\nu + \bits - 2}$, $x = \lfloor a / 2^{\nu - 2} \rfloor$ and
$\tilde q = \lfloor x m / 2^{\bits + 1} \rfloor$. Then $x < B$ and
$\lfloor a/n \rfloor - 1 \le \tilde q \le \lfloor a/n \rfloor$, so
$0 \le a - \tilde q n < 2n$.
:::

::: proof
The bound on $a$ gives $x < 2^{\bits}$. From $m \le 2^{\nu + \bits - 1}/n$
and $x \le a / 2^{\nu - 2}$ we get $x m / 2^{\bits + 1} \le a/n$. For the
other side, write $mn = 2^{\nu + \bits - 1} - 1 - \rho$ with $0 \le \rho < n$
([](#lem-nmod-m)) and use $x > a/2^{\nu - 2} - 1$:
$$
\frac{x m}{2^{\bits + 1}} > \frac{a m}{2^{\nu + \bits - 1}} - \frac{m}{2^{\bits + 1}}
= \frac{a}{n} - \frac{a (1 + \rho)}{n\, 2^{\nu + \bits - 1}} - \frac{m}{2^{\bits + 1}}
> \frac{a}{n} - \frac{1}{2} - \frac{1}{2},
$$
since $a < 2^{\nu + \bits - 2}$, $1 + \rho \le n$ and $m < B$. So $\tilde q$,
an integer greater than $a/n - 1$, is at least $\lfloor a/n \rfloor - 1$.
:::

Since $2n \le 2^{\bits}$, the value $a - \tilde q n$ is the low word of
$a_{lo} - \tilde q n$ with $a = a_{hi} B + a_{lo}$, and one correction
gives $a \bmod n$; also $\tilde q < 2^{\bits - 1}$. `NMOD_RED2_BOUNDED`
computes $x$ with a double-word shift, $\tilde q$ as the high word of $x m$
shifted right by one bit, then $a - \tilde q n$ and the correction: one high
product, one low product and one comparison, as in `NMOD_RED`, against a
full product, a low product and two corrections for a
[Möller–Granlund](algorithms/arithmetic/moller-granlund) step (`NMOD_RED2`,
which instead takes any $a < nB$). `NMOD_DIVREM2_BOUNDED` also gives the
quotient $\lfloor a/n \rfloor$, which is $\tilde q$ or $\tilde q + 1$ and
less than $2^{\bits - 1}$, with a branch-free correction; the
[radix](radix#multiplication) multiplications use it for their carries.

The bound holds for a sum of $k$ products of residues when
$\nu + \lceil \log_2 k \rceil \le \bits - 2$, since such a sum is below
$k\, 2^{2\nu}$. For $\nu = \bits$ the method fails, as $a - \tilde q n$
may not fit in a word; for $\nu \le 2$ the shift by $\nu - 2$ is not
defined.

A three-word integer $a < 2^{\nu + 2\bits - 2}$, for instance a sum of
fewer than $2^{2\bits - 2 - \nu}$ products of residues, is reduced in two
steps (`NMOD_RED3_BOUNDED`). Its top two words $\lfloor a/B \rfloor$ are
below $2^{\nu + \bits - 2}$, so the step above reduces them to some
$r < n$, with no reduction of the top word first; then $rB + a_0$, which
is congruent to $a$, is reduced by one
[Möller–Granlund](algorithms/arithmetic/moller-granlund) step
(`NMOD_RED2`). `NMOD_RED3` instead takes two Möller–Granlund steps.

These reductions apply where the length, the field size or the modulus
bounds the value: the final reductions of two- and three-limb
[dot products](algorithms/linear-algebra/dot-products#scalar), hence the
classical polynomial and matrix products; the unpacking of
[Kronecker substitution](algorithms/polynomials/kronecker#ks); in
`nmod_mat`, the row reductions, the triangular solves with delayed
reductions and the SIMD products by few columns; the basecase polynomial
division; and the loops that call `_nmod_mul_bounded`,
`_nmod_addmul_bounded` and `_nmod_fmma_bounded` ([](nmod#arithmetic)),
such as Horner evaluation, Taylor shifts and the inversion of a vector.

## Lemire–Kaser–Kurz reduction {#lemire}

The method of Lemire, Kaser and Kurz [@LKK2019] computes the remainder
directly from the fractional part of $x/n$. With $c = \lceil B/n \rceil$,
the low word $xc \bmod B$ approximates $B \{x/n\}$, where $\{\cdot\}$
denotes the fractional part; multiplying it by $n$ and keeping the high
word gives $x \bmod n$. The approximation is exact when $n$ and $x$ are
below $2^{\bits/2}$, which is the condition required by `n_mod_lemire`.
The method needs no correction step, but its range restriction makes it
useful only for half-word values. `n_lemire_precomp` computes
$c = \lfloor (B - 1)/n \rfloor + 1$ with one division;
`n_lemire_precomp_preinv` adds one to the constant of
`n_barrett_precomp_preinv` instead, which gives the same word for every
$n \ge 1$. For $n \ge 2$ this is $c$, since
$\lfloor (B - 1)/n \rfloor + 1 = \lceil B/n \rceil$ for all integers
$B, n \ge 1$.

::: functions n_mod_barrett n_mod_barrett_lazy n_barrett_precomp n_barrett_precomp_preinv n_mod2_preinv n_mod2_preinv_lazy n_ll_barrett_precomp_preinv n_ll_mod_barrett n_lll_mod_barrett n_mod_lemire n_lemire_precomp n_lemire_precomp_preinv NMOD_RED NMOD_RED2_BOUNDED NMOD_DIVREM2_BOUNDED NMOD_RED3_BOUNDED NMOD_RED3_FOLD_FULLWORD _nmod_red3_precomp :::
