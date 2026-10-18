---
title: Roots of polynomials over finite fields
summary: Roots in $\mathbb{F}_q$ by squarefree factorisation, the split $\gcd(f, x^{(q-1)/2} \pm 1)$ and random shifts, with powers of $x + c$ computed through Taylor shifts; roots modulo prime powers by Hensel lifting of all roots, and modulo composite $n$ by the Chinese remainder theorem.
complexity: Over $\mathbb{F}_q$, degree $n$: $O(M(n) \log q \log n)$ expected. Over $\Z/p^k\Z$: one root finding modulo $p$ plus $O(\log k)$ lifting steps per root, each one evaluation of $f$ and $f'$. Over $\Z/n\Z$: the product of the numbers of roots modulo the prime powers, which can be exponential in the number of prime factors.
order: 80
---

A root $a$ of $f$ in a field $K$ is the same as a linear factor $x - a$,
so finding all roots of $f \in \mathbb{F}_q[x]$ is the part of the
factorisation that concerns factors of degree 1
([](algorithms/polynomials/factoring-finite-fields)). It is much cheaper
than a complete factorisation, because only one distinct-degree step is
needed and the splitting does not require the trace map. Over the ring
$\Z/n\Z$ with $n$ composite the situation is different: a polynomial can
have more roots than its degree, and the roots are found prime by prime.

## Roots in a finite field {#finite-fields}

By [](algorithms/polynomials/factoring-finite-fields#thm-xqd) with $d =
1$,
$$
x^q - x = \prod_{a \in \mathbb{F}_q} (x - a),
$$
so for squarefree $g$ the product of the linear factors of $g$ is
$\gcd(g, x^q - x)$. For odd $q$ the nonzero elements split further:

::: lemma {#lem-squares}
For odd $q$,
$$
x^{(q-1)/2} - 1 = \prod_{a \in (\mathbb{F}_q^\times)^2} (x - a), \qquad
x^{(q-1)/2} + 1 = \prod_{a \in \mathbb{F}_q^\times \setminus
(\mathbb{F}_q^\times)^2} (x - a).
$$
For $q = 2^k$ and $T(x) = x + x^2 + \cdots + x^{2^{k-1}}$, $x^q - x = T(x)
(T(x) + 1)$, and $T(x)$, $T(x) + 1$ are the products of $x - a$ over the
elements of absolute trace 0 and 1 respectively.
:::

::: proof
$\mathbb{F}_q^\times$ is cyclic of even order $q - 1$, so exactly half of
its elements are squares, and these are the roots of $x^{(q-1)/2} - 1$
(Euler's criterion). The other half are roots of the second factor of
$x^{q-1} - 1 = (x^{(q-1)/2} - 1)(x^{(q-1)/2} + 1)$. For even $q$, $T(x)^2
= T(x) + x^q - x$ since squaring is additive, so $T (T + 1) = x^q - x$;
$T(a) \in \mathbb{F}_2$ is the trace of $a$.
:::

The same splitting applies to the shifted roots $a + c$, which gives a
randomised equal-degree factorisation for $d = 1$, due to Berlekamp
[@Ber1970] and Rabin [@Rab1980b]:

::: algorithm Roots in $\mathbb{F}_q$ {#alg-roots}
**Input:** $f \in \mathbb{F}_q[x]$ nonzero.

**Output:** the roots of $f$ with multiplicities.

1. Compute the squarefree factorisation $f = c \prod_i g_i^{\,i}$
   ([](algorithms/polynomials/squarefree)). For each $g = g_i$:
2. If $x \mid g$, output the root 0 and remove the factor $x$.
3. For odd $q$: $t \gets x^{(q-1)/2} \bmod g$; $u_1 \gets \gcd(g, t - 1)$
   and $u_2 \gets \gcd(g, t + 1)$. For $q = 2^k$: $t \gets T(x) \bmod g$;
   $u_1 \gets \gcd(g, t)$ and $u_2 \gets \gcd(g, t + 1)$.
4. Split $u_1$ and $u_2$ into linear factors: while a piece $u$ has
   degree at least 2, choose a random $c \in \mathbb{F}_q$ and replace $u$
   by $\gcd(u, (x + c)^{(q-1)/2} - 1)$ and its cofactor if this is a
   proper split (for even $q$, use $T(\alpha x + c)$ with random $\alpha
   \ne 0$).
5. Output the roots of the linear pieces with multiplicity $i$.
:::

Step 3 discards all factors of $g$ of degree at least 2 and already
splits the roots into two classes; this initial split is Rabin's. In step
4 a root $\rho$ of $u$ goes into the gcd if and only if $\rho + c$ is a
nonzero square, which for a uniformly random $c$ happens with probability
about $1/2$, independently enough for distinct roots that a piece of
degree $m$ is split with probability at least $4/9$
([](algorithms/polynomials/factoring-finite-fields#lem-split)) and
the expected depth of the recursion is $O(\log n)$. Each level costs
exponentiations of total degree $O(n)$, so the expected total cost is
$O(M(n) \log q \log n)$, dominated by the exponentiations.

**Powers of $x + c$.** An exponentiation $x^e \bmod u$ is cheaper than
$a^e \bmod u$ for a general $a$, since the multiplications by $x$ in the
square-and-multiply loop are shifts. FLINT uses this for the random
shifts when $q$ has more than 16 bits:

::: lemma {#lem-taylor-shift}
Let $u(x)$ have degree $m$ and let $R = x^e \bmod u(x - c)$. Then
$(x + c)^e \bmod u(x) = R(x + c)$.
:::

::: proof
Write $x^e = Q(x) u(x - c) + R(x)$ with $\deg R < m$ and substitute $x
\mapsto x + c$: $(x + c)^e = Q(x + c) u(x) + R(x + c)$, and $\deg R(x + c)
< m$.
:::

So $(x + c)^{(q-1)/2} \bmod u$ costs two Taylor shifts of degree $m$
([](algorithms/polynomials/composition#taylor-shift)) and one
exponentiation of $x$.

**Distinct nonzero roots.** When $f$ is known to split into distinct
linear factors $x - a$ with $a \ne 0$, as in the decoding of
Reed–Solomon codes and in sparse interpolation, the squarefree
factorisation is unnecessary, and step 3 decides at once whether $f$ has
this form: the degrees of $u_1$ and $u_2$ add up to $\deg f$ if and only
if it does. `nmod_poly_find_distinct_nonzero_roots` and
`fmpz_mod_poly_find_distinct_nonzero_roots` implement this test and the
splitting (see [](nmod_poly/roots)).

## Roots modulo prime powers {#prime-powers}

Let $f \in \Z[x]$ and $k \ge 2$. Every root of $f$ modulo $p^k$ reduces to
a root modulo $p^e$ for every $e \le k$, so the roots modulo $p^k$ can be
found by lifting all roots modulo $p$, doubling the precision in each
step. Unlike the classical form of Hensel's lemma [@Hen1908] (see also
[](algorithms/polynomials/hensel-lifting) for the lifting of
factorisations), this lifts roots that are multiple modulo $p$ as well.

::: lemma {#lem-hensel}
Let $f(a) \equiv 0 \pmod{p^e}$ and $e' = \min(k, 2e)$. For $t \in \Z$,
$$
f(a + p^e t) \equiv 0 \pmod{p^{e'}} \iff f'(a)\, t \equiv
-\frac{f(a)}{p^e} \pmod{p^{e' - e}}.
$$
:::

::: proof
By Taylor's formula $f(a + h) = f(a) + h f'(a) + h^2 r(a, h)$ with $r \in
\Z[x, h]$ (the coefficients $f^{(j)}/j!$ of the expansion are integral).
With $h = p^e t$ the last term is divisible by $p^{2e}$, hence by
$p^{e'}$, so $f(a + p^e t) \equiv f(a) + p^e t f'(a) \pmod{p^{e'}}$.
Dividing by $p^e$ gives the claim.
:::

The linear congruence $\alpha t \equiv \beta \pmod{p^s}$ has no solution
if $p^j = \gcd(\alpha, p^s)$ does not divide $\beta$, and otherwise
exactly $p^j$ solutions modulo $p^s$, an arithmetic progression with step
$p^{s-j}$ (computed by `fmpz_divides_mod_list`). Hence a simple root
modulo $p$ ($f'(a) \not\equiv 0$) lifts to exactly one root modulo $p^k$,
while a multiple root can lift to none or to many. The lifts $a + p^e t$
for $t$ modulo $p^{e'-e}$ are distinct modulo $p^{e'}$, and every root
modulo $p^{e'}$ arises from its reduction modulo $p^e$, so lifting every
root at every step enumerates all roots modulo $p^k$. If $f \equiv 0
\pmod p$, all $p$ residues are roots modulo $p$ and are lifted in the
same way.

The number of roots can be large: the roots of $x^2$ modulo $p^k$ are
the multiples of $p^{\lceil k/2 \rceil}$, of which there are
$p^{\lfloor k/2 \rfloor}$.

## Composite moduli {#composite}

For $n = \prod_j p_j^{k_j}$, the Chinese remainder theorem gives
$\Z/n\Z \cong \prod_j \Z/p_j^{k_j}\Z$, and the roots of $f$ modulo $n$
are the CRT combinations of all tuples of roots modulo the $p_j^{k_j}$.
Their number is the product of the numbers of roots modulo the prime
powers, which is zero as soon as one factor is zero, and can be far
larger than $\deg f$.

In $(\Z/n\Z)[x]$ the multiplicity of a root $a$ is defined as the largest
$e$ such that $(x - a)^e$ divides $f$; since $x - a$ is monic,
divisibility holds modulo $n$ if and only if it holds modulo every
$p_j^{k_j}$, so the multiplicity is the minimum of the multiplicities
modulo the prime powers, and each of these is found by repeated exact
division by $x - a$.

## Use in FLINT {#use}

`gr_poly_roots_finite_field` implements [](#alg-roots) for every finite
field context; step 4 is `_gr_poly_factor_equal_deg_with_frob` with $d =
1$, which also contains the Taylor shift trick. `nmod_poly_roots` and
`fmpz_mod_poly_roots` wrap it and return the roots as linear factors in
an `nmod_poly_factor_t` or `fmpz_mod_poly_factor_t`; the
`fq_*_poly_roots` functions of [](fq_poly_factor) and its variants are
analogous, with the field arithmetic of
[](algorithms/finite-fields/overview).

`nmod_poly_roots_factored` and `fmpz_mod_poly_roots_factored` take the
factorisation of the modulus and implement the prime-power lifting and the
CRT combination. They give up and return 0 when a list of roots would
exceed about $2^{32}$ entries ($2^{25}$ on 32-bit systems);
`fmpz_mod_poly_roots_factored_with_length_limit` takes the limit for the
CRT combination as a parameter. The generic method `gr_poly_roots` over
`nmod` and `fmpz_mod` contexts uses the finite field algorithm when the
modulus is prime and otherwise factors the modulus and calls
`fmpz_mod_poly_roots_factored_with_length_limit` with a limit of $10^6$
roots, returning `GR_UNABLE` if the modulus has more than 20 prime
factors (with multiplicity) or there are too many roots.

::: functions nmod_poly_roots nmod_poly_roots_factored fmpz_mod_poly_roots fmpz_mod_poly_roots_factored fmpz_mod_poly_roots_factored_with_length_limit nmod_poly_find_distinct_nonzero_roots fmpz_mod_poly_find_distinct_nonzero_roots :::

Rational roots of integer polynomials (`fmpz_poly_roots_fmpq`) are found
from the roots modulo a suitable prime by $p$-adic lifting, which is the
same lifting as above for simple roots.
