---
title: Arithmetic in unramified extensions
summary: How qadic represents the unramified extension Q_q = Q_p[X]/(f) and computes products with reduction modulo a sparse defining polynomial, inverses, powers, traces by Newton's identities, norms by resultants or as exp(Tr(log x)), and square roots including the Artin–Schreier step for p = 2.
complexity: Multiplication: $O(M(d \cdot N\log p))$ by Kronecker substitution plus $O(dt)$ for the reduction ($t$ terms in $f$); trace: $O(dt)$ operations modulo $p^N$; norm: one $(2d-1) \times (2d-1)$ determinant modulo $p^N$, or a logarithm, a trace and an exponential.
order: 40
---

An unramified extension of $\mathbf{Q}_p$ of degree $d$ is
$\mathbf{Q}_q = \mathbf{Q}_p[X]/(f)$, $q = p^d$, for a monic
$f \in \mathbf{Z}_p[X]$ whose reduction $\bar f$ modulo $p$ is
irreducible over $\mathbf{F}_p$. It depends only on $d$, its ring of
integers is $\mathbf{Z}_q = \mathbf{Z}_p[X]/(f)$, and its residue field
is $\mathbf{Z}_q / p\mathbf{Z}_q \cong \mathbf{F}_p[X]/(\bar f) =
\mathbf{F}_q$. In particular $p$ is still a uniformiser: every nonzero
element is $p^v u$ with $u \in \mathbf{Z}_q^\times$, and
$u \in \mathbf{Z}_q$ is a unit if and only if its reduction in
$\mathbf{F}_q$ is nonzero [@Rob2000, Ch. 3]. The Frobenius automorphism
and Teichmüller lifts are treated in [](algorithms/p-adic/frobenius).

## Representation {#representation}

A `qadic_ctx_t` stores $f$ with integer coefficients in $[0, p)$ as a
sparse list of terms $f = \sum_k a_k X^{j_k}$, $j_0 < \cdots <
j_{t-1} = d$, $a_{t-1} = 1$. By default (`qadic_ctx_init`) $f$ is the
Conway polynomial for $(p, d)$ if the database has it, otherwise the
polynomial chosen by `nmod_poly_minimal_irreducible` for word-size $p$,
or a random sparse irreducible polynomial for larger $p$; any monic
polynomial irreducible modulo $p$ can be supplied instead
([](algorithms/finite-fields/defining-polynomials)).

An element is a `padic_poly_t` $p^v g(X)$ with $\deg g < d$, where the
coefficients of $g$ are reduced modulo $p^{N - v}$ and not all
divisible by $p$; this is the fixed-precision model of
[](algorithms/p-adic/precision#fixed). Since $g$ has a coefficient that
is a unit of $\mathbf{Z}_p$, its image in $\mathbf{F}_q$ is nonzero, so
$v$ is the valuation of the element and $g$ represents a unit of
$\mathbf{Z}_q$.

## Multiplication and reduction {#mul}

The product of $p^{v} g$ and $p^{w} h$ is $p^{v + w} gh$; as
$\bar g \bar h \ne 0$ in the field $\mathbf{F}_q$, the product is again
canonical and no factor of $p$ has to be removed. `qadic_mul` multiplies
$g$ and $h$ in $\Z[X]$ (`_fmpz_poly_mul`, which uses Kronecker
substitution or other methods depending on the sizes,
[](algorithms/polynomials/multiplication-integers)) and reduces the
product modulo $f$ and $p^{N - v - w}$.

::: algorithm Reduction modulo a sparse monic polynomial {#alg-reduce}
**Input:** $R = \sum_{i < \ell} r_i X^i \in \Z[X]$ and
$f = X^d + \sum_{k < t - 1} a_k X^{j_k}$.

**Output:** $R \bmod f$, of length at most $d$.

1. For $i = \ell - 1, \ell - 2, \ldots, d$: for each $k < t - 1$,
   subtract $r_i a_k$ from $r_{i - d + j_k}$; then set $r_i = 0$.
:::

Each step removes the term $r_i X^i$ by subtracting
$r_i X^{i - d} f$, so the result is congruent to $R$ modulo $f$. The
reduction costs $(\ell - d)(t - 1)$ multiply-subtracts, which is small
compared with the product when $f$ is sparse, as Conway polynomials and
the minimal-weight polynomials used by default usually are. The
coefficients are reduced modulo $p^{N - v}$ once, after the reduction
modulo $f$ (`_fmpz_mod_poly_reduce`). The same functions serve
`qadic_reduce` and the powering, inversion, Frobenius and series
evaluation routines of `qadic`. The finite field modules use the same
sparse elimination ([](algorithms/finite-fields/reduction)).

`qadic_pow` uses left-to-right binary exponentiation with this
reduction.

## Inversion {#inverse}

The inverse of $p^v g$ is $p^{-v} g^{-1}$, and $g^{-1}$ is computed by
the Newton iteration $z \gets 2z - gz^2$ in $\mathbf{Z}_q$
([](algorithms/p-adic/hensel-lifting#inverse)), starting from the
inverse of $\bar g$ in $\mathbf{F}_q$, which `_qadic_inv` computes by
the extended Euclidean algorithm in $\mathbf{F}_p[X]$ modulo $\bar f$.
The total cost is a constant number of multiplications in $\mathbf{Z}_q$
at precision $N$.

## Trace {#trace}

The trace $\operatorname{Tr}(x)$ of $x \in \mathbf{Q}_q$ is the trace of
multiplication by $x$ as a $\mathbf{Q}_p$-linear map, equivalently the
sum of the images of $x$ under the $d$ embeddings into an algebraic
closure. For $x = \sum_{i < d} c_i X^i$ it is
$\sum_i c_i S_i$, where $S_i = \sum_{f(\theta) = 0} \theta^i$ are the
power sums of the roots of $f$.

::: proposition Newton's identities {#prop-newton}
For monic $f = X^d + \sum_{m = 1}^{d} b_m X^{d - m}$ the power sums
satisfy $S_0 = d$ and, for $1 \le i < d$,

$$
S_i = -\Bigl(i\, b_i + \sum_{m = 1}^{i - 1} b_m S_{i - m}\Bigr).
$$
:::

`_qadic_trace` evaluates this recurrence on the sparse representation,
in which only the $m = d - j_k$ with $a_k \ne 0$ contribute, so that
$S_0, \ldots, S_{\ell-1}$ for an element of length $\ell$ cost
$O(\ell t)$ operations modulo $p^N$, and then forms
$\sum c_i S_i$. As $\mathbf{Q}_q/\mathbf{Q}_p$ is unramified,
$\operatorname{Tr}(\mathbf{Z}_q) = \mathbf{Z}_p$; the trace of
$p^v g$ is $p^v \operatorname{Tr}(g)$, which may have valuation larger
than $v$, and `qadic_trace` canonicalises it.

## Norm {#norm}

The norm $N(x)$ is the determinant of multiplication by $x$, or the
product of the conjugates of $x$; $N(p^v g) = p^{dv} N(g)$, and
$N(g)$ is a unit for a unit $g$. FLINT has two methods.

The *resultant method* uses $N(g) = \operatorname{Res}(f, g)$ for
monic $f$, which holds because $\operatorname{Res}(f, g) =
\prod_{f(\theta) = 0} g(\theta)$. `_qadic_norm_resultant` builds the
Sylvester matrix of $f$ and $g$, of dimension at most $2d - 1$, and
computes its determinant modulo $p^N$ with `fmpz_mod_mat_det`, which
uses Gaussian elimination when the pivots it meets are invertible and
falls back to the division-free Berkowitz algorithm otherwise.

The *analytic method* applies when $g \equiv 1 \pmod p$ (modulo 4 if
$p = 2$).

::: proposition {#prop-norm-analytic}
If $x \in 1 + p\mathbf{Z}_q$ ($1 + 4\mathbf{Z}_q$ if $p = 2$), then
$N(x) = \exp(\operatorname{Tr}(\log x))$.
:::

::: proof
The automorphisms $\sigma$ of $\mathbf{Q}_q/\mathbf{Q}_p$ are
continuous, so they commute with the power series:
$\log(\sigma x) = \sigma(\log x)$. Hence $\log N(x) = \log \prod_\sigma
\sigma(x) = \sum_\sigma \sigma(\log x) = \operatorname{Tr}(\log x)$.
Both $\log x$ and its trace lie in $p\mathbf{Z}_p$ (in $4\mathbf{Z}_2$
for $p = 2$), where $\exp$ converges and inverts $\log$, and
$N(x) \in 1 + p\mathbf{Z}_p$ (respectively $1 + 4\mathbf{Z}_2$), so
$\exp(\log N(x)) = N(x)$.
:::

`_qadic_norm_analytic` computes the logarithm in $\mathbf{Z}_q$, its
trace and the exponential in $\mathbf{Z}_p$
([](algorithms/p-adic/exp-log)); its cost is quasi-linear in $d$ and
$N$. `_qadic_norm` chooses the analytic method when it applies and
$4 \lfloor \log_2 N \rfloor^2 \lfloor \log_2 d \rfloor < d^3$, and the
resultant method otherwise; for elements of $\mathbf{Q}_p$ (length 1) it
returns $g^d$.

## Square roots {#sqrt}

A nonzero $p^v u$ is a square in $\mathbf{Q}_q$ if and only if $v$ is
even and $u$ is a square in $\mathbf{Z}_q$. For odd $p$, Hensel's lemma
reduces the latter to $\bar u$ being a square in $\mathbf{F}_q$.
`qadic_sqrt` computes a square root in $\mathbf{F}_q$ by
$x^{(q+1)/4}$ if $q \equiv 3 \pmod 4$, by Atkin's formula if
$q \equiv 5 \pmod 8$ (with $y = x^{(q-1)/4} = \pm 1$, the root is
$x^{(q+3)/8}$ or $2^{(q-1)/4} x^{(q+3)/8}$), and otherwise by the
Tonelli–Shanks algorithm with a quadratic non-residue found by
enumerating the elements of $\mathbf{F}_q$ in a fixed order
([](algorithms/finite-fields/frobenius-norm-trace)). It then lifts the
inverse square root as in [](algorithms/p-adic/hensel-lifting#sqrt).

For $p = 2$ every element of $\mathbf{F}_q$ is a square, $\bar x^{q/2}$
being the square root of $\bar x$, but not every unit of $\mathbf{Z}_q$
is a square. Let $r = u^{-1}$ and let $s$ lift the square root of
$\bar r$. A lift $s + 2X$ is an inverse square root of $u$ modulo 8 if
and only if $r - s^2 \equiv 4(X^2 + sX) \pmod 8$. This requires
$r \equiv s^2 \pmod 4$, a condition independent of the choice of $s$,
and then, with $c = (r - s^2)/4$ and $X = sY$, the Artin–Schreier
equation

$$
Y^2 + Y = c\,\bar u \quad \text{in } \mathbf{F}_q .
$$

The map $Y \mapsto Y^2 + Y$ is $\mathbf{F}_2$-linear with kernel
$\mathbf{F}_2$, so the equation is solved by linear algebra over
$\mathbf{F}_2$: `qadic_sqrt` builds the $d \times d$ matrix of the map
and its LUP decomposition, and `_qadic_char2_sqrt_precomp_init` stores
that decomposition so that `_qadic_char2_sqrt_with_precomp` can reuse it
for many square roots in the same field. If the equation has no
solution, $u$ is not a square. Otherwise the lift is refined by the
2-adic inverse square root iteration of
[](algorithms/p-adic/hensel-lifting#lem-rsqrt).

## Use in FLINT {#use}

- Reduction: `_fmpz_poly_reduce`, `_fmpz_mod_poly_reduce`,
  `qadic_reduce`; multiplication and powering: `qadic_mul`,
  `_qadic_pow`, `qadic_pow`.
- Inversion: `_qadic_inv`, `qadic_inv`.
- Trace and norm: `_qadic_trace`, `qadic_trace`, `_qadic_norm`,
  `qadic_norm`, `_qadic_norm_resultant`, `qadic_norm_resultant`,
  `_qadic_norm_analytic`, `qadic_norm_analytic`.
- Square roots: `qadic_sqrt`, `_qadic_char2_sqrt_precomp_init`,
  `_qadic_char2_sqrt_with_precomp`, `_qadic_char2_sqrt_precomp_clear`.
- Exponential and logarithm: [](algorithms/p-adic/exp-log).
