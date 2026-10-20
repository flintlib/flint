---
title: Newton iteration for power series and division
summary: Inverse, quotient, square root, reciprocal square root, exponential and logarithm of power series by Newton iteration with middle products and the Karp–Markstein trick, and polynomial division with remainder via reversed polynomials, as implemented for `nmod_poly`.
complexity: With products of cost $M(n)$ and middle products: inverse in about $2M(n)$, power series quotient in about $\tfrac{5}{2}M(n)$, division with remainder of length $2n$ by length $n$ in about $\tfrac{7}{2}M(n)$; square roots and exponentials in a small constant times $M(n)$.
order: 20
---

Newton's method doubles the number of correct terms of a power series
approximation with each step. Since the last step dominates, computing
an inverse, a quotient, a square root or an exponential to $n$ terms
costs a small constant times one multiplication of length $n$. Power
series inversion by Newton iteration is due to Sieveking [@Sie1972] and
Kung [@Kun1974]; Brent and Kung [@BrentKung1978] treated general
equations and elementary functions. Textbook treatments are
[@GG2003, Ch. 9] and [@BZ2011, §4.2]; Hanrot and Zimmermann
[@HZ2004] and Bernstein [@Ber2008] discuss the constants.

Throughout, the coefficients lie in $R = \Z/p\Z$ for a word-size
modulus $p$ (we write $p$ rather than $n$ to keep $n$ for the
precision), $\beta$ is the bit length of $p$, "to precision $n$" means
modulo $x^n$, and $M(n)$ is the cost of multiplying two
polynomials of length $n$. We assume, as is true for all methods in
[](algorithms/polynomials/multiplication), that $M(n)/n$ is
nondecreasing, so $\sum_{i \ge 1} M(\lceil n/2^i \rceil) \le M(n) + O(\log n)$.
With the FFT, a low product of length $n$ and a middle product with
output length $n$ cost about $M(n)$ each
([](algorithms/polynomials/middle-product)).

## Generic version {#generic}

The iterations described here are implemented once in `gr_poly` and are
valid over any commutative ring, not only $\Z/p\Z$: each step uses ring
operations and one inversion of the constant coefficient (the leading
coefficient of the divisor for division with remainder), and the
algorithms return `GR_DOMAIN` if that coefficient is not invertible.
Over real and complex balls they give enclosures of the exact result,
since every operation does. The generic dispatchers use the basecase
for short inputs or when the ring has only classical multiplication
(inverse and quotient of power series below length 9, division when the
quotient or the divisor has length at most 10, exponential below length
10), Newton iteration otherwise, and middle products when the ring
provides them ([](algorithms/polynomials/overview#division)). When a
leading coefficient is not invertible, as for division by $2x + 1$ over
$\Z$, division with remainder and power series division fall back to
algorithms that divide exactly at each step and succeed whenever the
quotient lies in the ring. The rings `nmod`, `fmpz_mod`, `mpn_mod` and
the `fq` family override the dispatchers to use their own basecase
kernels and cutoffs, as `nmod` does with the tables quoted below; the
ball rings divide through power series with error bounds.

## The general step {#principle}

::: lemma Quadratic convergence {#lem-newton}
Let $\Phi(x, y) \in R[[x]][y]$, let $y^* \in R[[x]]$ satisfy
$\Phi(x, y^*) = 0$ with $\Phi_y(x, y^*)$ invertible, and let $m \ge 1$.
If $y \equiv y^* \pmod{x^m}$ and $u \equiv \Phi_y(x, y)^{-1} \pmod{x^m}$,
then $y - u \, \Phi(x, y) \equiv y^* \pmod{x^{2m}}$.
:::

::: proof
Let $\delta = y - y^* = O(x^m)$. By Taylor expansion,
$\Phi(x, y) = \Phi_y(x, y^*) \delta + O(\delta^2)$ and
$\Phi_y(x, y) = \Phi_y(x, y^*) + O(\delta)$, so
$\Phi(x, y) = \Phi_y(x, y) \delta + O(x^{2m})$. Since
$\Phi(x, y) = O(x^m)$, replacing $\Phi_y(x, y)^{-1}$ by $u$ changes
$u \Phi(x, y)$ by $O(x^{2m})$. Hence $y - u \Phi(x, y) = y - \delta +
O(x^{2m}) = y^* + O(x^{2m})$.
:::

FLINT runs Newton iterations from the target precision downwards: it
records $a_0 = n$, $a_{i+1} = \lceil a_i / 2 \rceil$ until $a_k$ is below
a cutoff, computes the solution to precision $a_k$ by a basecase
(quadratic) algorithm, and then lifts from $a_{i+1}$ to $a_i$. This
avoids overshooting: no step computes more terms than are needed. The
same scheme lifts $p$-adic solutions, with $p$ in place of $x$
([](algorithms/p-adic/hensel-lifting)).

## Inverse {#inverse}

Let $Q(0)$ be a unit and $g \equiv 1/Q \pmod{x^m}$. Then
$e = 1 - Qg \equiv 0 \pmod{x^m}$, and since $Qg = 1 - e$,

$$
Q \cdot g (1 + e) = 1 - e^2, \qquad Q \cdot g (1 + e + e^2) = 1 - e^3,
$$

so $g(1 + e)$ is correct to $2m$ terms and $g(1 + e + e^2)$ to $3m$.
The first is the Newton step of [](#lem-newton) for $\Phi = 1/y - Q$.
The low $m$ coefficients of $Qg$ are $1, 0, \ldots, 0$, so only the
window $w = (Qg)_{[m, n)}$ has to be computed, which is a middle
product, and $e = -x^m w$.

::: algorithm Power series inverse {#alg-inv}
**Input:** $Q$ with $Q(0)$ a unit, a precision $n$.

**Output:** $g = Q^{-1} \bmod x^n$.

1. If $n$ is below the cutoff, compute $g$ by the recurrence
   $g_0 = Q_0^{-1}$, $g_i = -g_0 \sum_{j=1}^{i} Q_j g_{i-j}$ and return.
2. Compute $g \gets Q^{-1} \bmod x^m$ recursively, with $m = \lceil n/2 \rceil$
   (or $m = \lceil n/3 \rceil$ in a third-order step).
3. $w \gets$ coefficients $m, \ldots, n-1$ of $(Q \bmod x^n) \cdot g$
   (middle product).
4. Second order ($n \le 2m$): $g \gets g - x^m (g w \bmod x^{n-m})$.

   Third order ($n \le 3m$): $g \gets g - x^m \bigl(g (w - x^m w^2) \bmod x^{n-m}\bigr)$.
:::

::: proposition {#prop-inv-cost}
With middle products, [](#alg-inv) with second-order steps costs about
$2M(n)$; computing the full product $Qg$ instead of the middle product
gives the classical bound of about $3M(n)$ [@GG2003, Th. 9.4].
:::

::: proof
The step from $m$ to $n = 2m$ costs one middle product with output length
$m$ and one low product of length $m$, about $2M(n/2)$, and the sum over
all steps is at most about $2M(n)$. Forming the low $2m$ coefficients of
the product of $Q \bmod x^{2m}$ and $g$ instead costs about
$M(n)$ (a product with output length $2m$), giving $M(n) + M(n/2)$ per
step and about $3M(n)$ in total.
:::

Exploiting the structure of the FFT further (reusing transforms between
steps and between the two products) lowers the constant; the best known
value is $\tfrac{13}{9}M(n)$, due to Harvey [@Har2011], who also gives
$\tfrac{4}{3}M(n)$ for the square root. FLINT does not use these
FFT-specific methods.

The third-order step costs one more short squaring but triples the
precision, so there are fewer steps. FLINT uses it for $n \ge 8192$ in
rings of finite characteristic, where it is 5–10% faster with the FFT. (Over $\Z$ and $\mathbb{Q}$ the coefficients of
$e^2$ can be much larger than those of the result.)

**In FLINT.** `nmod_poly_inv_series` uses the recurrence (with
`_nmod_vec_dot_rev` and delayed reduction) for $\operatorname{len}(Q)
\le 10$, and otherwise calls the `gr_poly` code with the `nmod` context,
which uses the recurrence if $n$ is below a cutoff depending on $\beta$
(`NMOD_BITS(mod)`), or if $Q$ is short, and `_gr_poly_inv_series_newton`
otherwise. With `fft_small` the cutoff is about 40 for $\beta \le 7$,
80–300 for $8 \le \beta \le 21$, about 400 for $22 \le \beta \le 28$,
and 230–450 above. The recurrence is used for any $n$ when
$\operatorname{len}(Q)$ is below a second table (between 39 and 141 for
$\beta \ge 23$; smaller moduli do not use this rule), since it costs
only $O(n \operatorname{len}(Q))$. If $Q$ is a binomial $Q_0 + c x^k$, the inverse is a geometric series and is
written down directly. `nmod_poly_inv_series_newton` is the same function.

## Quotient of power series {#div-series}

To compute $A/B \bmod x^n$, one could invert $B$ to $n$ terms and
multiply. The Karp–Markstein trick [@KM1997] instead merges the
multiplication by $A$ into the last Newton step.

::: lemma Karp–Markstein {#lem-karp-markstein}
Let $g \equiv B^{-1}$ and $T \equiv A g \pmod{x^m}$. Then
$T + g (A - BT) \equiv A/B \pmod{x^{2m}}$.
:::

::: proof
$A - B(T + g(A - BT)) = (A - BT)(1 - Bg)$, and both factors are
$O(x^m)$ because $T \equiv A/B$ and $g \equiv 1/B \pmod{x^m}$. Dividing
by $B$ gives the claim.
:::

::: algorithm Power series quotient {#alg-div-series}
**Input:** $A$, $B$ with $B(0)$ a unit, a precision $n$; $m = \lceil n/2 \rceil$.

**Output:** $A/B \bmod x^n$.

1. $g \gets B^{-1} \bmod x^m$ by [](#alg-inv).
2. $T \gets A g \bmod x^m$ (low product).
3. $u \gets A_{[m, n)} - (BT)_{[m, n)}$ (middle product; the low $m$
   coefficients of $A - BT$ vanish).
4. Return $T + x^m (g u \bmod x^{n-m})$.
:::

This costs an inverse to $n/2$ terms plus three products of length
about $n/2$, about $M(n) + \tfrac{3}{2}M(n) = \tfrac{5}{2}M(n)$,
instead of $2M(n) + M(n) = 3M(n)$.

**In FLINT.** `nmod_poly_div_series` uses the recurrence
$q_i = B_0^{-1}(A_i - \sum_{j \ge 1} B_j q_{i-j})$ for
$\operatorname{len}(B) \le 5$, and otherwise `_gr_poly_div_series_newton`
([](#alg-div-series)) above a cutoff that, with `fft_small`, is about
80–340 for $\beta \le 20$, about 600 for $21 \le \beta \le 28$ and
320–490 above, with a separate table for short $B$.

## Division with remainder {#divrem}

Let $A = QB + R$ with $\deg A = a$, $\deg B = b$, $\deg R < b$ and
$\deg Q = a - b$. Write $\operatorname{rev}_k(f) = x^k f(1/x)$.

::: lemma {#lem-reverse}
$\operatorname{rev}_a(A) = \operatorname{rev}_{a-b}(Q) \operatorname{rev}_b(B)
+ x^{a-b+1} \operatorname{rev}_{b-1}(R)$. Consequently, if the leading
coefficient of $B$ is a unit,
$\operatorname{rev}_{a-b}(Q) \equiv \operatorname{rev}_a(A) /
\operatorname{rev}_b(B) \pmod{x^{a-b+1}}$.
:::

::: proof
Substitute $1/x$ in $A = QB + R$ and multiply by $x^a = x^{a-b} x^b$;
the term $R(1/x)$ is multiplied by $x^a = x^{a - b + 1} x^{b-1}$. The
constant term of $\operatorname{rev}_b(B)$ is the leading coefficient
of $B$.
:::

So the quotient is a power series quotient of length
$\ell_Q = a - b + 1$, which needs only the top $\ell_Q$ coefficients of
$A$ and of $B$, and the remainder is $R = (A - QB) \bmod x^b$, a low
product.

::: algorithm Division with remainder {#alg-divrem}
**Input:** $A$ of length $\ell_A$, $B$ of length $\ell_B \le \ell_A$ with
invertible leading coefficient; $\ell_Q = \ell_A - \ell_B + 1$.

**Output:** $Q$, $R$ with $A = QB + R$, $\operatorname{len}(R) < \ell_B$.

1. Let $A'$ be the reversal of the top $\ell_Q$ coefficients of $A$ and
   $B'$ the reversal of the top $\min(\ell_B, \ell_Q)$ coefficients of $B$.
2. $Q \gets$ reversal of $A'/B' \bmod x^{\ell_Q}$ ([](#alg-div-series)).
3. $R \gets A - (QB \bmod x^{\ell_B - 1})$ (low product, then a
   subtraction of $\ell_B - 1$ coefficients).
:::

By [](#lem-reverse) the algorithm is correct. For $\ell_A = 2n$ and
$\ell_B = n$ it costs about $\tfrac{5}{2}M(n) + M(n) = \tfrac{7}{2}M(n)$,
compared with $n^2$ operations for schoolbook division.

If many dividends are reduced modulo the same $B$, the power series
inverse of $\operatorname{rev}(B)$ can be computed once (it is usually
computed to $\ell_B$ terms; $\ell_Q$ terms suffice). Then each
quotient with $\ell_Q \le \ell_B - 1$ is one low product, and the
division costs two products: this is how `_nmod_poly_div_newton_n_preinv`,
`_nmod_poly_divrem_newton_n_preinv` and everything built on them
(`nmod_poly_mulmod_preinv`, `nmod_poly_powmod_ui_binexp_preinv`,
`nmod_poly_compose_mod_brent_kung_preinv`) work. For a sparse divisor (at
most five nonzero coefficients below the leading one, and
$\ell_B > 20$), `_nmod_poly_divrem_newton_n_preinv` instead eliminates
directly, at a cost of a few operations per coefficient of the quotient.

**In FLINT.** `nmod_poly_divrem` and `nmod_poly_rem` use schoolbook
division (`_nmod_poly_divrem_basecase_preinv1`, with delayed reduction)
if $\ell_A - \ell_B < 80$ or $\ell_B$ is less than 400 (for
$\beta \le 27$ on x86-64 with AVX2; 200 without AVX2) or 80 (for
$\beta > 27$), and [](#alg-divrem) otherwise. Without `fft_small`,
divide-and-conquer division is used instead for $\beta \ge 16$,
$\ell_B \ge 1024$ and $\ell_A \le 16384$. `nmod_poly_div` uses schoolbook
division if $\ell_B \le 15$ or $\ell_A - \ell_B \le 15$, and steps 1 and 2
of [](#alg-divrem) otherwise.

## Square roots {#sqrt}

Let $2$ be invertible and let $h(0)$ have a square root in $R$.
For the reciprocal square root, let $g \equiv h^{-1/2} \pmod{x^m}$ and
$e = 1 - hg^2 = O(x^m)$. Since $h = g^{-2}(1 - e)$,

$$
h^{-1/2} = g (1 - e)^{-1/2} = g \bigl(1 + \tfrac{1}{2} e + \tfrac{3}{8} e^2\bigr) + O(x^{3m}).
$$

FLINT uses this third-order iteration throughout
(`_gr_poly_rsqrt_series_newton`): each step computes the square $g^2$
(low product), the window of $h g^2$ (middle product), the short square
$e^2$ and one low product with $g$, and triples the precision. When
$n \le 2m$ the term in $e^2$ is dropped.

The square root itself is obtained from the reciprocal square root to
half precision by a Karp–Markstein step:

::: lemma {#lem-sqrt-km}
Let $g \equiv h^{-1/2}$ and $s = gh \bmod x^m$, so that
$s \equiv h^{1/2} \pmod{x^m}$. Then
$s + \tfrac{1}{2} g (h - s^2) \equiv h^{1/2} \pmod{x^{2m}}$.
:::

::: proof
$h^{1/2} - s = (h - s^2)/(h^{1/2} + s)$. Here $h - s^2 = O(x^m)$ and
$h^{1/2} + s \equiv 2h^{1/2}$, so $1/(h^{1/2} + s) \equiv g/2 \pmod{x^m}$.
:::

For short $h$, FLINT uses instead J. C. P. Miller's recurrence for
powers [@Knu1997, §4.7]: $y = h^\alpha$ satisfies $h y' = \alpha h' y$,
which gives each coefficient of $y$ from the previous ones in
$O(\operatorname{len}(h))$ operations. It divides by $1, \ldots, n - 1$,
so it is used only when $n$ is less than the modulus, and only for
$\operatorname{len}(h)$ below a tuned bound (with `fft_small`, about 40
to 190 depending on $\beta$, and only for $\beta \ge 19$) and $8 \operatorname{len}(h) < n$
(square root) or $2 \operatorname{len}(h) < n$ (reciprocal square root).

`nmod_poly_sqrt_series` and `nmod_poly_invsqrt_series` call these
algorithms. `nmod_poly_sqrt`, the square root of a polynomial that is a
perfect square, removes an even power of $x$, computes the series square
root to half the length plus one, and checks the result by comparing
the high half of its square (`_nmod_poly_mulhigh`) with the input.

## Exponential and logarithm {#exp}

For $h(0) = 0$, $\log$ and $\exp$ of power series are defined if
$1, 2, \ldots, n - 1$ are invertible, which for prime modulus $p$ means
$n \le p$. The logarithm is

$$
\log f = \log f(0) + \int \frac{f'}{f},
$$

one power series division and an integration (`_gr_poly_log_series`,
called by `nmod_poly_log_series`; if $f$ is a binomial
$f_0 + c x^k$ the coefficients are written down directly).

For the exponential, FLINT uses a Newton iteration that computes
$f = \exp(h)$ and $g = \exp(-h) = 1/f$ together [@HZ2004; @BZ2011, §4.2.5].

::: lemma {#lem-exp}
Let $f \equiv \exp(h)$ and $g \equiv 1/f \pmod{x^m}$, with $f$ of length
$m$. Then $h' f - f' \equiv 0 \pmod{x^{m-1}}$, and with

$$
\varepsilon = \int g \, (h' f - f'),
$$

we have $\varepsilon = O(x^m)$ and $f (1 + \varepsilon) \equiv \exp(h) \pmod{x^{2m}}$.
:::

::: proof
Since $f \equiv \exp(h)$, $f' \equiv h' f \pmod{x^{m-1}}$. Note
$f(0) = 1$. Let $\eta = h - \log f = \int (h' f - f')/f$, which is
$O(x^m)$. Then $\eta' - \varepsilon' = (h' f - f')(1/f - g)$ is the
product of an $O(x^{m-1})$ and an $O(x^m)$ series, so
$\varepsilon = \eta + O(x^{2m})$. Finally $f = \exp(h - \eta)$, so
$f (1 + \varepsilon) = \exp(h) e^{-\eta} (1 + \eta) + O(x^{2m}) =
\exp(h) + O(x^{2m})$.
:::

::: algorithm Exponential {#alg-exp}
**Input:** $h$ with $h(0) = 0$, a precision $n$.

**Output:** $f = \exp(h) \bmod x^n$ (and optionally $g = \exp(-h)$).

1. If $n$ is small, compute $f$ by the basecase and $g = 1/f$ by
   [](#alg-inv), and return.
2. Recursively compute $f$, $g$ to precision $m = \lceil n/2 \rceil$.
3. $t \gets$ coefficients $m-1, \ldots, n-2$ of $h' f$ (middle product;
   these equal those of $h' f - f'$, as $f'$ has degree $m - 2$).
4. $u \gets g t \bmod x^{n-m}$, integrated as the coefficients of
   $x^m, \ldots, x^{n-1}$: $u_k \gets u_k / (m + k)$.
5. $f \gets f + x^m (f u \bmod x^{n-m})$.
6. If $g$ is needed at this precision: $w \gets$ coefficients $m, \ldots,
   n - 1$ of $f g$ (middle product) and $g \gets g - x^m (g w \bmod x^{n-m})$.
:::

Each step costs two middle products and three low products of length
about $n/2$; the last step skips the update of $g$ unless $\exp(-h)$ is
requested. Faster but more complicated methods are known: $\tfrac{7}{3}M(n)$
by van der Hoeven [@vdH2010] and $\tfrac{13}{6}M(n)$ by Harvey
[@Har2009b], both relying on FFT-specific techniques.

The basecase uses the recurrence $k f_k = \sum_{j=1}^{k} j h_j f_{k-j}$
from $f' = h' f$, costing $O(n \operatorname{len}(h))$. An intermediate
method (`_gr_poly_exp_series_basecase_mul`) splits $h = a + x^m b$ with
$m = \lceil n/3 \rceil$ and uses
$\exp(h) \equiv \exp(a)(1 + x^m b + x^{2m} b^2/2) \pmod{x^n}$.
`nmod_poly_exp_series` uses the basecase if $\operatorname{len}(h)$ is
below a tuned bound (230–490 with `fft_small`), the split method below a
second bound (470–1491), and [](#alg-exp) otherwise, with the second bound
as its basecase cutoff. `nmod_poly_exp_expinv_series` computes the
exponential and then inverts it; `nmod_poly_sinh_series` and
`nmod_poly_cosh_series` use it. The other trigonometric and inverse
trigonometric series functions are computed by their `gr_poly`
counterparts.

## Use in FLINT {#use}

All the Newton iterations live in `gr_poly` and are called by
`nmod_poly` through the `gr` context for $\Z/n\Z$, which supplies the
`nmod_poly` multiplication and middle product and the tuned cutoffs
quoted above (tables indexed by $\beta$).

| Operation | `nmod_poly` function | Iteration |
|---|---|---|
| inverse | `nmod_poly_inv_series` | `_gr_poly_inv_series_newton` |
| quotient | `nmod_poly_div_series` | `_gr_poly_div_series_newton` |
| division with remainder | `nmod_poly_divrem`, `nmod_poly_div`, `nmod_poly_rem` | `_gr_poly_divrem_newton`, `_gr_poly_div_newton` |
| with precomputed inverse | `nmod_poly_divrem_newton_n_preinv`, `nmod_poly_div_newton_n_preinv` | one low product per quotient |
| reciprocal square root | `nmod_poly_invsqrt_series` | `_gr_poly_rsqrt_series_newton` |
| square root | `nmod_poly_sqrt_series`, `nmod_poly_sqrt` | `_gr_poly_sqrt_series_newton` |
| exponential | `nmod_poly_exp_series` | `_gr_poly_exp_series_newton` |
| logarithm | `nmod_poly_log_series` | via `_gr_poly_div_series` |

The `nmod_poly` cutoffs between the basecase algorithms and Newton
iteration are listed in [](nmod_poly/inv_series#algorithms) for series
inversion and division and in [](nmod_poly/divrem#algorithms) for
division with remainder.

Power series composition and reversion, which also use Newton iteration
for long series, are treated in [](algorithms/polynomials/composition).

::: functions nmod_poly_inv_series nmod_poly_div_series nmod_poly_divrem nmod_poly_divrem_newton_n_preinv nmod_poly_sqrt_series nmod_poly_exp_series nmod_poly_log_series :::
