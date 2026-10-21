---
title: Factoring polynomials over finite fields
summary: Squarefree, distinct-degree and equal-degree factorisation over $\mathbb{F}_q$; the Cantor–Zassenhaus, Berlekamp and Kaltofen–Shoup algorithms, the baby-step giant-step distinct-degree factorisation with modular composition, the trace and minimal-polynomial splittings, and how the generic implementation in `gr_poly` chooses between them.
complexity: For $f$ of degree $n$ over $\mathbb{F}_q$: Cantor–Zassenhaus $O(n M(n) \log q)$ expected; Kaltofen–Shoup with Brent–Kung composition $O(n M(n) + n^{1/2} C(n) + M(n) \log q \log n)$ expected, $C(n) = O(n^{1/2} M(n) + n^{(\omega+1)/2})$; Berlekamp $O(n^\omega + n M(n) \log n + M(n) \log q \log n)$ expected.
order: 75
---

Let $\mathbb{F}_q$ be a finite field with $q = p^k$ elements and let $f
\in \mathbb{F}_q[x]$ have degree $n$. Every factorisation algorithm used
by FLINT reduces the problem to three steps, which go back to Berlekamp
[@Ber1967] [@Ber1970] and to Cantor and Zassenhaus [@CanZas1981]:

1. *Squarefree factorisation* writes $f = c \prod_i g_i^{e_i}$ with
   squarefree, pairwise coprime $g_i$ and distinct $e_i$
   ([](algorithms/polynomials/squarefree)).
2. *Distinct-degree factorisation* (DDF) splits a squarefree $g$ into
   $g = \prod_d g_{[d]}$, where $g_{[d]}$ is the product of the
   irreducible factors of $g$ of degree $d$.
3. *Equal-degree factorisation* (EDF) splits each $g_{[d]}$ into its
   irreducible factors, by a randomised algorithm.

The algorithms differ in how they perform steps 2 and 3. The
*Cantor–Zassenhaus algorithm* computes the DDF with one exponentiation
by $q$ per degree; the *Kaltofen–Shoup algorithm* [@KalSho1998]
[@Sho1995] uses a baby-step giant-step method in which the
exponentiations are replaced by modular compositions; *Berlekamp's
algorithm* skips the DDF and splits $g$ with the kernel of a linear map.
Textbook accounts are [@GG2003, Ch. 14] and [@Sho2009, Ch. 20].

All of them are implemented once, generically, in `gr_poly`
(`gr_poly_factor_finite_field` and the functions it calls), for any
`gr` context that is a finite field: `nmod`, `fmpz_mod`, `fq`,
`fq_nmod` and `fq_zech`. The modules `nmod_poly_factor`,
`fmpz_mod_poly_factor` and the `fq_*_poly_factor` modules are wrappers.
Throughout, $M(n)$ denotes the cost of multiplying polynomials of degree
$n$ over $\mathbb{F}_q$, counted in field operations, and $C(n)$ the cost
of one modular composition $a(b) \bmod f$ for $a, b$ of degree less than
$n$ ([](algorithms/polynomials/composition#modular)).

## Finite fields and the Frobenius map {#frobenius}

::: theorem {#thm-xqd}
For every $d \ge 1$,
$$
x^{q^d} - x = \prod_{\substack{g \text{ monic irreducible} \\ \deg g \mid d}} g
$$
in $\mathbb{F}_q[x]$.
:::

::: proof
The polynomial $x^{q^d} - x$ is squarefree, since its derivative is $-1$.
Its roots in an algebraic closure are the $q^d$ elements of
$\mathbb{F}_{q^d}$. A monic irreducible $g$ of degree $e$ divides it if
and only if a root $\alpha$ of $g$ lies in $\mathbb{F}_{q^d}$, that is, if
and only if $\mathbb{F}_q(\alpha) \cong \mathbb{F}_{q^e}$ is a subfield of
$\mathbb{F}_{q^d}$, which happens exactly when $e \mid d$ [@LN1997, Ch. 2].
:::

::: corollary {#cor-gcd}
If $g$ is squarefree, $\gcd(g, x^{q^d} - x)$ is the product of the
irreducible factors of $g$ whose degree divides $d$. In particular, $g$
of degree $n$ is irreducible if and only if $\gcd(g, x^{q^d} - x) = 1$
for $1 \le d \le n/2$.
:::

The gcd is computed with $x^{q^d} \bmod g$, never with the huge polynomial
$x^{q^d} - x$ itself. The powers $x^{q^d} \bmod g$ are linked by the
Frobenius endomorphism $a \mapsto a^q$ of $\mathbb{F}_q[x]$:

::: lemma {#lem-frob}
For $a \in \mathbb{F}_q[x]$, $a(x)^q = a(x^q)$. Consequently, if $h_i =
x^{q^i} \bmod g$, then $h_{i+j} = h_i(h_j) \bmod g$ and, more generally,
$a^{q^i} \bmod g = a(h_i) \bmod g$.
:::

::: proof
The map $a \mapsto a^q$ is a ring endomorphism in characteristic $p$
($(a + b)^p = a^p + b^p$) and fixes the coefficients, since $c^q = c$ for
$c \in \mathbb{F}_q$. Hence $(\sum c_j x^j)^q = \sum c_j x^{jq}$. Applying
this $i$ times gives $a^{q^i} = a(x^{q^i})$, and reducing modulo $g$
gives $a^{q^i} \equiv a(h_i)$. With $a = x^{q^j}$ this is $h_{i+j} \equiv
h_j(h_i) \equiv h_i(h_j)$ (the two orders are equal by symmetry).
:::

So there are two ways to step from $h_i$ to $h_{i+1}$: an exponentiation
$h_{i+1} = h_i^q \bmod g$, which costs about $\tfrac{3}{2} \log_2 q$
multiplications modulo $g$, and a composition $h_{i+1} = h_i(h_1) \bmod
g$, which costs $C(n)$ independently of $q$. With the Brent–Kung
algorithm and the matrix of powers of $h_1$ precomputed, a composition
costs about $\sqrt{n}$ multiplications modulo $g$ plus a product of an
$\sqrt{n} \times \sqrt{n}$ and an $\sqrt{n} \times n$ matrix, so
composition wins when $\log q$ is large compared to $\sqrt{n}$. (The
Frobenius map on elements of $\mathbb{F}_q$ itself is discussed in
[](algorithms/finite-fields/frobenius-norm-trace).) The
generic code uses composition when $q$ has more than $\tfrac{3}{4}(\lfloor
\sqrt{n} \rfloor + 1)$ bits.

## Squarefree factorisation and deflation {#squarefree}

Squarefree factorisation needs only gcds, divisions, derivatives and, in
characteristic $p$, $p$-th roots of coefficients; it costs $O(M(n) \log
n)$ operations per gcd and is described in
[](algorithms/polynomials/squarefree). It reduces the problem to monic
squarefree polynomials, and all steps below assume this.

Before the squarefree factorisation, `_gr_poly_factor_finite_field`
*deflates* the input: if $F(x) = G(x^m)$ for some $m \ge 2$ (all exponents
of nonzero terms are multiples of $m$), it factors $G$, and then factors
$g_i(x^m)$ for each factor $g_i$ of $G$ (without deflating again). Each
factor of $F$ divides exactly one $g_i(x^m)$, so this splits the problem
into independent problems of smaller size; polynomials like $x^m - a$ or
$x^{2m} + a x^m + b$ benefit most.

## Distinct-degree factorisation {#ddf}

### The classical algorithm {#ddf-classical}

::: algorithm Classical distinct-degree factorisation {#alg-ddf-classical}
**Input:** $g$ monic squarefree of degree $n$.

**Output:** the pairs $(g_{[d]}, d)$ with $g_{[d]} \ne 1$.

1. $v \gets g$, $h \gets x$, $i \gets 0$.
2. While $\deg v \ge 2(i + 1)$:
    1. $i \gets i + 1$; $h \gets h^q \bmod v$.
    2. $u \gets \gcd(v, h - x)$; if $u \ne 1$, output $(u, i)$, set $v
       \gets v / u$ and $h \gets h \bmod v$.
3. If $v \ne 1$, output $(v, \deg v)$.
:::

By induction on $i$, when $u$ is computed $v$ has no factors of degree
less than $i$, so by [](#cor-gcd) $u$ is the product of the factors of
degree exactly $i$. When $\deg v < 2(i + 1)$, the remaining $v$ has no
factor of degree at most $i$ and therefore cannot be a product of two or
more factors; it is irreducible. The loop runs at most $n/2$ times, and
each iteration costs an exponentiation, $O(M(n) \log q)$, and a gcd, for
a total of $O(n M(n) \log q)$.

### Baby steps and giant steps {#ddf-bsgs}

Shoup [@Sho1995] and Kaltofen and Shoup [@KalSho1998] remove the factor
$\log q$ by grouping the degrees into intervals. Choose $\ell \approx
\sqrt{n/2}$ and $m = \lceil n / (2\ell) \rceil$, so that the $m$ intervals
$(\ell j, \ell (j+1)]$, $0 \le j < m$, cover all degrees up to $n/2$.

::: lemma {#lem-interval}
Let $h_i = x^{q^i} \bmod g$ and $H_j = x^{q^{\ell(j+1)}} \bmod g$. An
irreducible factor $u$ of $g$ of degree $d$ divides $H_j - h_i$ if and
only if $d \mid \ell(j + 1) - i$.
:::

::: proof
Modulo $u$, $H_j - h_i \equiv x^{q^{\ell(j+1)}} - x^{q^i} = (x^{q^{\ell(j+1)
- i}} - x)^{q^i}$, since the $q^i$-th power map is additive. Since $u$ is irreducible, it divides a $q^i$-th power if and only if
it divides the base, and by [](#thm-xqd) that happens if and only if $d
\mid \ell(j+1) - i$.
:::

::: algorithm Baby-step giant-step distinct-degree factorisation {#alg-ddf-bsgs}
**Input:** $g$ monic squarefree of degree $n \ge 2$, and $h_1 = x^q \bmod
g$.

**Output:** the pairs $(g_{[d]}, d)$ with $g_{[d]} \ne 1$.

1. *Baby steps:* $h_i \gets x^{q^i} \bmod g$ for $0 \le i \le \ell$.
2. *Coarse DDF:* $v \gets g$. For $j = 0, 1, \ldots, m - 1$:
    1. $H_j \gets h_\ell$ if $j = 0$, else $H_j \gets H_{j-1}(h_\ell) \bmod
       v$ (a modular composition with the precomputed matrix of powers
       of $h_\ell$).
    2. $I_j \gets \prod_{i=0}^{\ell - 1} (H_j - h_i) \bmod v$, omitting
       the terms for degrees $\ell(j+1) - i > \tfrac12 \deg v$.
    3. $F_j \gets \gcd(v, I_j)$ (or $F_j \gets v$ if $I_j = 0$); $v \gets
       v / F_j$.
    4. Stop if $\deg v < 2(\ell(j + 1) + 1)$.
3. If $v \ne 1$, output $(v, \deg v)$.
4. *Fine DDF:* for each $j$ with $F_j \ne 1$: if $j = 0$ or $\deg F_j >
   \ell(j+1)$, then for $i = \ell - 1, \ldots, 0$: $u \gets \gcd(F_j,
   H_j - h_i)$; if $u \ne 1$, output $(u, \ell(j+1) - i)$ and set $F_j
   \gets F_j / u$. Otherwise output $(F_j, \deg F_j)$.
:::

::: theorem {#thm-ddf-bsgs}
[](#alg-ddf-bsgs) is correct. With Brent–Kung modular composition it
uses $\ell + m \approx \sqrt{2n}$ modular compositions (or
exponentiations, see below), at most $\ell m \approx n/2$ multiplications
modulo $g$, and $O(\sqrt{n})$ gcds, in total
$$
O\bigl(n M(n) + \sqrt{n}\, C(n)\bigr)
$$
operations, plus $O(M(n) \log q)$ for $h_1$.
:::

::: proof
By [](#lem-interval), the factors of $v$ that divide $I_j$ are those
whose degree divides some integer in $(\ell j, \ell(j+1)]$. By induction
on $j$, when $F_j$ is computed every factor of $v$ has degree greater
than $\ell j$ (smaller degrees $d$ divide $d$ itself, which lies in an
earlier interval), so $F_j$ is the product of the factors with degree in
$(\ell j, \ell(j+1)]$. Omitting the terms for degrees above $\tfrac12 \deg
v$ loses nothing, because if $v$ had a factor of such a degree $d$ and
were reducible, it would also have a factor of degree $\le \deg v - d <
d$, already removed. The stopping rule and step 3 are justified as in the
classical algorithm. Reducing $H_{j-1}$ and the matrix of powers of
$h_\ell$ modulo the new $v$ after a split is valid since $v$ divides the
old modulus.

In the fine phase, for $j \ge 1$ the factors of $F_j$ have degrees in
$(\ell j, \ell(j+1)]$, an interval in which no degree divides another
(since $2(\ell j + 1) > \ell(j+1)$); so $\gcd(F_j, H_j - h_i)$ is exactly
the product of the factors of degree $\ell(j+1) - i$. If $\deg F_j \le
\ell(j+1)$, $F_j$ has a single factor, since two factors would have
total degree greater than $2\ell j \ge \ell(j+1)$. For $j = 0$ the
degrees $1, 2, \ldots, \ell$ are processed in increasing order and the
factors found are removed, so the gcd for degree $d$ only contains
factors of degree exactly $d$, as in the classical algorithm.

The cost bound counts the compositions of steps 1 and 2.1, the $\ell$
multiplications per interval of step 2.2, and one gcd per interval and
per degree in the fine phase, $O(\sqrt{n} M(n) \log n)$ in total, which
is dominated by the compositions.
:::

The baby steps are computed in `_gr_poly_iterated_frobenius_preinv`,
by repeated exponentiation by $q$ if $q$ is small (see
[](#frobenius)) and otherwise by modular composition using a doubling
scheme: given $h_1, \ldots, h_s$, the vector composition
$h_{s+i} = h_i(h_s)$, $1 \le i \le s$, is a single matrix product with the
powers of $h_s$ (`gr_poly_preinv_compose_mod_brent_kung_vec`). The giant
steps always use composition with the precomputed matrix of powers of
$h_\ell$.

With $C(n) = O(\sqrt{n} M(n) + n^{(\omega+1)/2})$ the bound of
[](#thm-ddf-bsgs) is $O(n M(n) + n^{(\omega+2)/2})$ and no longer
depends on $q$ except through $h_1$; this is what makes the algorithm
fast for large $q$. Kaltofen and Shoup also show how to make the
interval products subquadratic, which together with fast matrix
multiplication gives a factorisation algorithm with $O(n^{1.815} \log q)$
operations [@KalSho1998]; the asymptotically fastest known method, based
on the modular composition algorithm of Kedlaya and Umans, uses
$O(n^{1.5 + o(1)} + n^{1 + o(1)} \log q)$ bit operations per field
operation [@KU2011]. FLINT uses neither of these refinements.

## Equal-degree factorisation {#edf}

Let now $g$ be monic, squarefree, of degree $n = rd$, the product of $r
\ge 2$ irreducible factors $u_1, \ldots, u_r$ of degree $d$. By the
Chinese remainder theorem
$$
R = \mathbb{F}_q[x]/(g) \cong \prod_{i=1}^r \mathbb{F}_q[x]/(u_i) \cong
\mathbb{F}_{q^d}^r, \qquad a \mapsto (a \bmod u_1, \ldots, a \bmod u_r).
$$
The *Berlekamp subalgebra* $B = \{ b \in R : b^q = b \}$ consists of the
elements whose components all lie in $\mathbb{F}_q$, so $B \cong
\mathbb{F}_q^r$. Every $b \in B$ whose components are not all equal gives
a proper factor $\gcd(g, b - c)$ for some $c \in \mathbb{F}_q$; the
algorithms below produce such $b$ at random and then find a $c$, or a
polynomial in $b$ playing the same role.

::: lemma Trace map {#lem-trace}
The map $\operatorname{Tr}(a) = a + a^q + \cdots + a^{q^{d-1}}$ is an
$\mathbb{F}_q$-linear surjection $R \to B$. If $a \in R$ is uniformly
random, then $\operatorname{Tr}(a)$ is uniformly random in $B$.
:::

::: proof
Componentwise, $\operatorname{Tr}$ is the trace of $\mathbb{F}_{q^d}$ over
$\mathbb{F}_q$, which is $\mathbb{F}_q$-linear and surjective [@LN1997,
§2.3]. A surjective linear map sends the uniform distribution to the
uniform distribution, since all fibres are cosets of the kernel.
:::

::: lemma Splitting probability {#lem-split}
Let $q$ be odd and let $b \in B$ be uniformly random. Then $\gcd(g,
b^{(q-1)/2} - 1)$ is a proper factor of $g$ with probability at least
$4/9$. For $q = 2^k$ and the absolute trace $T(b) = b + b^2 + \cdots +
b^{2^{k-1}}$, $\gcd(g, T(b))$ is a proper factor with probability $1 -
2^{1-r} \ge 1/2$.
:::

::: proof
For odd $q$, a component $c \in \mathbb{F}_q$ of $b$ satisfies
$c^{(q-1)/2} = 1$ if and only if $c$ is a nonzero square, which happens
with probability $s = (q - 1)/(2q)$, independently for the $r$
components. The gcd is the product of the $u_i$ for which this holds, so
it is trivial with probability $s^r + (1 - s)^r \le s^2 + (1 - s)^2 = 1 -
2s(1 - s)$, and $2s(1-s) = (q^2 - 1)/(2q^2) \ge 4/9$ for $q \ge 3$. For
even $q$, $T$ is the trace of $\mathbb{F}_q$ over $\mathbb{F}_2$ on each
component, so the components of $T(b)$ are independent and uniform in
$\{0, 1\}$, and the gcd is trivial only if they are all 0 or all 1.
:::

::: algorithm Equal-degree splitting {#alg-edf}
**Input:** $g$ as above, $h_1 = x^q \bmod g$.

**Output:** the irreducible factors of $g$.

1. If $\deg g = d$, return $g$.
2. Choose $a \in R$ at random and compute $b \gets \operatorname{Tr}(a)$;
   repeat if $b$ is constant.
3. $t \gets b^{(q-1)/2} - 1 \bmod g$ for odd $q$, $t \gets T(b) \bmod g$
   for even $q$; $u \gets \gcd(g, t)$.
4. If $u$ is a proper factor, recurse on $u$ and $g / u$ (with $h_1$
   reduced modulo each); otherwise go to step 2.
:::

By [](#lem-split), each attempt succeeds with probability at least $4/9$,
and because each factor goes to either side independently, the expected
depth of the recursion is $O(\log r)$. The cost of an attempt is that of
the trace, an exponentiation, $O(M(n) \log q)$, and a gcd. FLINT computes
the trace in one of three ways (`_gr_poly_trace`):

- for small $q$ (the same test as in [](#frobenius)), as $a + a^q + \cdots$
  with $d - 1$ exponentiations by $q$;
- for $d < 16$, with $d - 1$ compositions $a^{q^{i+1}} = a^{q^i}(h_1)$
  using the precomputed matrix of powers of $h_1$;
- otherwise by the doubling method of von zur Gathen and Shoup
  [@GatSho1992]: with $T_s(a) = a + a^q + \cdots + a^{q^{s-1}}$, the
  identities $T_{2s}(a) = T_s(a) + T_s(a)(x^{q^s})$, $T_{s+1}(a) = a +
  T_s(a)(h_1)$ and $x^{q^{2s}} = x^{q^s}(x^{q^s})$ compute $T_d(a)$ along
  the binary expansion of $d$ with $O(\log d)$ compositions.

**Roots ($d = 1$).** Every element of $R$ is then in $B$, so no trace is
needed, and $b = x + c$ with random $c$ suffices (for even $q$, $b =
\alpha x + c$ with random $\alpha \ne 0$): the components of $b$ are the
values $\rho + c$ at the roots $\rho$ of $g$, so $\gcd(g, (x +
c)^{(q-1)/2} - 1)$ collects the roots for which $\rho + c$ is a nonzero
square. For $q$ of more than 16 bits, the
power is computed as $(x + c)^e \bmod g = R(x + c)$ where $R = x^e \bmod
g(x - c)$, because powers of $x$ modulo a polynomial are cheaper (the
multiplications by $x$ in the exponentiation are shifts); this costs two
Taylor shifts. See [](algorithms/polynomials/root-finding).

**Splitting by the minimal polynomial.** When $q$ is large and the
number $r$ of factors is moderate, it is cheaper to split $g$ into all
its factors at once [@Mar2014]. If $b \in B$ has components $c_1, \ldots,
c_r$, its minimal polynomial over $\mathbb{F}_q$ is $\mu(y) = \prod_{c
\in C} (y - c)$, where $C$ is the set of distinct values $c_i$, and
$$
g = \prod_{c \in C} \gcd(g, b - c).
$$
FLINT computes the sequence $s_i = w(b^i \bmod g)$, $0 \le i < 2r$, for a
random linear form $w$, obtains its minimal polynomial with the
Berlekamp–Massey algorithm [@Mas1969] (which equals $\mu$ unless $w$ is
degenerate, in which case the attempt is repeated), finds the roots $C$
of $\mu$ by the $d = 1$ case above, and splits $g$ by divide and conquer:
with $\mu_1 = \prod_{c \in C_1} (y - c)$ for the first half $C_1$ of the
roots, $\gcd(g, \mu_1(b) \bmod g)$ is the product of the $u_i$ with $c_i
\in C_1$, and $\mu_1(b) \bmod g$ is a modular composition. This uses $2r$
multiplications modulo $g$ and no exponentiation by $(q-1)/2$ for the
splitting element. A random $b$ has $r$ distinct components with
probability at least $1 - r^2/(2q)$, in which case the factorisation is
complete after one round. The generic code uses this method when $d >
1$, $q \ge 2^7$ and $r \le 8 \lambda(d + 1) + \lambda(q)$, where
$\lambda$ denotes the bit length.

## Berlekamp's algorithm {#berlekamp}

Berlekamp's algorithm [@Ber1967] [@Ber1970] does not need the DDF. Let
$g$ be monic squarefree of degree $n$ with $r$ irreducible factors (of
arbitrary degrees). The map $\phi(a) = a^q - a$ on $R = \mathbb{F}_q[x]/(g)$
is $\mathbb{F}_q$-linear, and its matrix $Q - I$ in the basis $1, x,
\ldots, x^{n-1}$ has as column $i$ the coefficients of $x^{qi} \bmod g$
minus the unit vector $e_i$.

::: theorem {#thm-berlekamp}
$\ker \phi = B$ has dimension $r$. In particular $g$ is irreducible if
and only if the nullity of $Q - I$ is 1.
:::

::: proof
By the Chinese remainder theorem $R \cong \prod_i \mathbb{F}_{q^{d_i}}$,
and $a^q = a$ holds in $\mathbb{F}_{q^{d_i}}$ exactly for $a \in
\mathbb{F}_q$. So $\ker \phi \cong \mathbb{F}_q^r$.
:::

`gr_poly_factor_berlekamp` computes $x^q \bmod g$, the $n$ columns
$x^{qi} \bmod g$ by repeated multiplication, and a basis of the kernel
with `gr_mat_nullspace`. If the nullity is 1, $g$ is irreducible.
Otherwise it takes random linear combinations $b$ of the basis vectors
(uniform in $B$), tries $\gcd(g, b)$ and then $\gcd(g, b^{(q-1)/2} - 1)$
(or $\gcd(g, T(b))$ for even $q$), and recurses on the two parts of the
first proper split, recomputing the matrix for each part. The kernel
costs $O(n^\omega)$, the matrix $O(n M(n))$, and each split attempt
$O(M(n) \log q)$; since random splits are balanced on average, the
expected total is $O(n^\omega + n M(n) \log r + M(n) \log q \log r)$.
The cubic linear algebra makes the algorithm uncompetitive for large $n$
in FLINT, where it serves mainly as an independent check.

## The complete algorithms and FLINT's choice {#choice}

All three algorithms begin with `gr_poly_factor_squarefree` and treat
each squarefree factor separately.

- `gr_poly_factor_cantor_zassenhaus` runs [](#alg-ddf-classical) on each
  squarefree factor, computing $h \gets h^q \bmod v$ by binary
  exponentiation, and calls the equal-degree factorisation
  (`_gr_poly_factor_equal_deg_with_frob`) on each $g_{[d]}$ with $d <
  \deg g_{[d]}$, passing $x^q \bmod g_{[d]}$ obtained from the first
  exponentiation.
- `gr_poly_factor_kaltofen_shoup` computes $x^q \bmod g$ once, runs
  [](#alg-ddf-bsgs) (`_gr_poly_factor_distinct_deg_with_frob`), and
  calls the equal-degree factorisation on each piece with $x^q$ reduced
  modulo the piece.
- `gr_poly_factor_berlekamp` is described above.

`gr_poly_factor_finite_field` (after deflation) uses Cantor–Zassenhaus
for degree
$$
n < 125 / (\log_2 q)^{1.3}
$$
and Kaltofen–Shoup otherwise. The rule fits measured crossovers for
dense random polynomials over `nmod`: about $n = 125$ for $q = 2$, 70 for
$q = 3$, 32 for $q = 7$, 16 for $q = 31$, 8 for $q \approx 2^8$ and 4 for
$q \approx 2^{16}$. For small $q$ the exponentiation by $q$ is a few
squarings and the classical DDF has the smaller constant; for large $q$
the $\log q$ factor dominates and the modular compositions win.

**Parallelism.** When threads are available (`flint_set_num_threads`)
and the context is thread-safe, polynomials of degree above 200
(`GR_POLY_FACTOR_THREADED_CUTOFF`) use
`_gr_poly_factor_distinct_deg_with_frob_threaded`: the giant steps are
computed in blocks of $B$ (the number of threads) independent
compositions $H_{j+B} = H_j(H_{B-1})$ with a precomputed matrix for
$H_{B-1}$, the interval polynomials of a block in parallel, and the fine
DDF in parallel over the intervals; $v$ is only reduced after each block,
which costs a little extra work when there are factors of small degree.
The equal-degree factorisation processes the pieces on the top of its
work stack in parallel, each with its own random state. Matrix products
use the threads of `gr_mat_mul`.

**Randomness.** Each top-level function initialises its own
`flint_rand_t` with the default seed, so results, including the order of
the factors, are reproducible from run to run (but not across FLINT
versions). The factors are returned in no particular order.

## Use in FLINT {#use}

The generic implementations are `gr_poly_factor_finite_field` (with
`_gr_poly_factor_finite_field` for a fixed algorithm),
`gr_poly_factor_cantor_zassenhaus`, `gr_poly_factor_kaltofen_shoup`,
`gr_poly_factor_berlekamp`, `gr_poly_factor_squarefree`,
`gr_poly_factor_distinct_deg`, `gr_poly_factor_equal_deg` and
`gr_poly_factor_equal_deg_prob`, with the helpers
`gr_poly_iterated_frobenius_preinv`,
`_gr_poly_factor_distinct_deg_with_frob` and
`_gr_poly_factor_equal_deg_with_frob`. The `gr` method
`GR_METHOD_POLY_FACTOR` of the finite field contexts points to
`gr_poly_factor_finite_field`. The type-specific modules wrap these
functions without copying the polynomials:

::: functions nmod_poly_factor nmod_poly_factor_with_cantor_zassenhaus nmod_poly_factor_with_berlekamp nmod_poly_factor_with_kaltofen_shoup nmod_poly_factor_squarefree nmod_poly_factor_distinct_deg nmod_poly_factor_equal_deg fmpz_mod_poly_factor fmpz_mod_poly_factor_cantor_zassenhaus fmpz_mod_poly_factor_berlekamp fmpz_mod_poly_factor_kaltofen_shoup fmpz_mod_poly_factor_distinct_deg fmpz_mod_poly_factor_equal_deg :::

The [`fq_poly_factor`](fq_poly_factor), `fq_nmod_poly_factor` and
`fq_zech_poly_factor` modules are generated from the same template and
wrap the same functions; the arithmetic of the finite fields themselves
is described in [](algorithms/finite-fields/overview). The modular factorisations used by the Zassenhaus algorithm
over $\Z$ ([](algorithms/polynomials/factoring-integers)) and by the
multivariate factorisation code
([](algorithms/multivariate/modp-factoring)) call `nmod_poly_factor` and
`fmpz_mod_poly_factor`. Irreducibility tests are described in
[](algorithms/polynomials/irreducibility), and root finding, which is
the case $d = 1$ with a cheaper first split, on
[](algorithms/polynomials/root-finding).
