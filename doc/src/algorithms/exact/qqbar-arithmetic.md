---
title: Arithmetic with algebraic numbers
summary: How qqbar computes x + y, x - y, xy and x/y exactly: composed polynomials from power sums, factorisation and numerical selection of the right factor, closed forms and certified norms for quadratic operands, half-degree polynomials for conjugates, linear algebra in a common number field or cyclotomic field found by LLL, and guessed minimal polynomials with a proof of correctness.
complexity: The general method computes a composed polynomial of degree $N = \deg x \deg y$ with $O(M(N))$ operations on rational power series and factors it over $\Z$; the structured methods replace the factorisation by LLL searches of dimension about $\max(\deg x, \deg y)$ and linear algebra in a field of that degree.
order: 20
---

For algebraic numbers $x$ and $y$ given by minimal polynomials $P$ and
$Q$ and isolating balls (see [](algorithms/exact/qqbar-representation)),
the number $w = x \circ y$ for $\circ \in \{+, -, \times, /\}$ is again
algebraic, and [](qqbar) computes its minimal polynomial and an
isolating ball. Operations with a rational operand are affine
substitutions $P((t - b)/a)$ and cost little (`qqbar_scalar_op`). What
follows is about two irrational operands, which are handled by
`qqbar_binary_op`. Its default strategy, classical since Loos
[@Loo1983], computes a polynomial of degree $\deg x \deg y$ that
annihilates $w$, factors it and selects the right factor; the other
strategies avoid this factorisation, which is the bottleneck when the
composed polynomial has many factors.

## Composed polynomials {#composed}

Let $P$ have roots $a_1, \ldots, a_m$ and $Q$ roots $b_1, \ldots, b_n$.
The *composed sum* and *composed product*

$$
P \oplus Q = \prod_{i,j} (t - a_i - b_j), \qquad
P \otimes Q = \prod_{i,j} (t - a_i b_j)
$$

have degree $mn$ and vanish at $x + y$ and $xy$. Differences and
quotients reduce to these by replacing $Q(t)$ by $Q(-t)$ (roots $-b_j$)
or by its reversal $t^n Q(1/t)$ (roots $1/b_j$). Up to the factor $L =
\operatorname{lc}(P)^n \operatorname{lc}(Q)^m$ they are resultants,
$L \cdot (P \oplus Q)(t) = \operatorname{Res}_s(P(s), Q(t - s))$, and
therefore have integer coefficients after multiplication by $L$.
Instead of a resultant, FLINT computes them with power sums, following
Bostan, Flajolet, Salvy and Schost [@BFSS2006].

::: lemma Newton sums {#lem-newton}
Let $F \in \mathbb{Q}[t]$ have degree $d$, roots $c_1, \ldots, c_d$
and power sums $p_k = \sum_i c_i^k$. With the reversed polynomial
$\tilde F(t) = t^d F(1/t)$,

$$
\sum_{k \ge 0} p_k t^k = \frac{t^{d-1} F'(1/t)}{\tilde F(t)}, \qquad
\tilde F(t) = \tilde F(0) \exp\Bigl(-\sum_{k \ge 1} \frac{p_k}{k} t^k\Bigr).
$$
:::

::: proof
Up to the constant $\operatorname{lc}(F)$, $\tilde F(t) = \prod_i (1 -
c_i t)$. Its logarithmic derivative is $\tilde F'/\tilde F = -\sum_i c_i
/ (1 - c_i t) = -\sum_{k \ge 1} p_k t^{k-1}$; integrating and
exponentiating gives the second formula. The first follows from
$t^{d-1} F'(1/t) / \tilde F(t) = \sum_i 1/(1 - c_i t)$, which is the
reversal of $F'/F = \sum_i 1/(t - c_i)$.
:::

The power sums of the composed polynomials are obtained from those of
$P$ and $Q$:

$$
\sum_{i,j} (a_i b_j)^k = p_k(P) \, p_k(Q), \qquad
\sum_{i,j} (a_i + b_j)^k = \sum_{l=0}^{k} \binom{k}{l} p_l(P) \, p_{k-l}(Q).
$$

The first is a termwise (Hadamard) product of the power sum series;
the second says that the exponential generating functions $\sum_k p_k
t^k / k!$ multiply. Hence `qqbar_fmpz_poly_composed_op` computes the
first $mn + 1$ power sums of $P$ and $Q$ by a power series division,
combines them by a termwise product or by a product of exponential
generating functions (multiplying by factorials before and after, a
Borel transform), and recovers the composed polynomial by a power
series exponential, all with `fmpq_poly` arithmetic. With fast
multiplication this costs $O(M(mn))$ arithmetic operations on rational
numbers, although the coefficients grow to about $mn$ times the size of
those of $P$ and $Q$.

When $x$ and $y$ are distinct roots of the same irreducible $P$ of degree
$d$, for example complex conjugates, $P \oplus P$ has the factor
$\prod_i (t - 2a_i)$ and every root $a_i + a_j$ with $i \ne j$ twice, so a
polynomial of half the degree suffices. `qqbar_fmpz_poly_symmetric_composed_op`
computes the polynomials of degree $d(d-1)/2$ with roots $a_i + a_j$,
$a_i a_j$ or $(a_i - a_j)^2$ for $i < j$ from the power sums

$$
\sum_{i<j} (a_i a_j)^k = \frac{p_k^2 - p_{2k}}{2}, \qquad
\sum_{i<j} (a_i + a_j)^k = \frac{1}{2} \Bigl(\sum_{l} \binom{k}{l} p_l p_{k-l} - 2^k p_k\Bigr),
$$

and the even part of the exponential generating function $E(t) E(-t)$
for $(a_i - a_j)^2$, where $E(t) = \sum p_k t^k/k!$. The difference $x -
y$ is then a root of $V(t^2)$, where $V$ is the minimal polynomial of
$(x - y)^2$.

## Selecting the factor {#select}

The minimal polynomial of $w$ is an irreducible factor of the composed
polynomial $H$, which is factored over $\Z$ (see
[](algorithms/polynomials/factoring-integers)). `_qqbar_binary_op_select_factor`
then computes a ball $W$ containing $w$ from the enclosures of $x$ and
$y$ and evaluates every irreducible factor $f_k$ on $W$. Distinct
irreducible polynomials have no common root, so at most one $f_k$
vanishes at $w$, and as the enclosures of $x$ and $y$ are refined (at
doubling precision) the values of all other factors on $W$ eventually
exclude zero. When exactly one factor has a value containing zero and
[](algorithms/exact/qqbar-representation#thm-newton) certifies that $W$
isolates a root of it, that factor is the minimal polynomial and the
certified ball is the enclosure of $w$.

## Quadratic operands {#quadratic}

If both operands are quadratic, $x = a_1 + b_1 \sqrt{D_1}$ and
$y = a_2 + b_2 \sqrt{D_2}$ with rational $a_k, b_k$ and principal square roots,
the result lies in $\mathbb{Q}(\sqrt{D_1}, \sqrt{D_2})$, which is a
quadratic field if $D_1 D_2$ is a square and a biquadratic field
otherwise. Its minimal polynomial is written down in closed form, of
degree at most 2 or 4; no factorisation is needed.

If only $y = a + b\sqrt{D}$ is quadratic and $\deg x = d \ge 3$, then $x
= L(w)$ for a polynomial $L$ of degree 1 with coefficients in $K =
\mathbb{Q}(\sqrt{D})$ (for example $L(z) = z - y$ for $w = x + y$).
Write $G(z) = P(L(z)) = U(z) + \sqrt{D} V(z)$ with $U, V \in
\mathbb{Q}[z]$. The norm $H = U^2 - D V^2 = G \overline{G}$ has degree
$2d$ and rational coefficients and vanishes at $w$.

::: lemma Norms of irreducible polynomials {#lem-norm}
If $\sqrt{D} \notin \mathbb{Q}(x)$, then $H$ is a power of an
irreducible polynomial in $\mathbb{Q}[z]$; in particular, $H$ is
irreducible if it is squarefree.
:::

::: proof
Since $[K : \mathbb{Q}] = 2$ and $K \not\subseteq \mathbb{Q}(x)$, the
fields $K$ and $\mathbb{Q}(x)$ are linearly disjoint, so $[K(x) : K] =
d$ and $P$ remains irreducible over $K$. A linear substitution preserves
irreducibility, so $G$ is irreducible over $K$, with root $w$. The
minimal polynomial $f$ of $w$ over $\mathbb{Q}$ is divisible by $G$ in
$K[z]$, hence also by $\overline{G}$, and every root of $H = G
\overline{G}$ is a conjugate of $w$ over $\mathbb{Q}$. Thus $H$ is a
product of copies of $f$ (up to a constant).
:::

The hypothesis is certified modularly. Choose a prime $p$, not dividing
$2D \operatorname{lc}(P) \operatorname{lc}(H)$, such that $P \bmod p$ is
squarefree and has an irreducible factor of odd degree $f$, $D$ is a
quadratic nonresidue modulo $p$, and $H \bmod p$ is squarefree of full
degree. Since $P \bmod p$ is squarefree, $p$ divides neither the
discriminant of $P$ nor the index of $\Z[x]$ in the ring of integers of
$\mathbb{Q}(x)$, so by Dedekind's theorem [@Neu1999, I.8.3] the factor
of degree $f$ corresponds to a prime $\mathfrak{p}$ of $\mathbb{Q}(x)$
with residue field $\mathbb{F}_{p^f}$. If $\sqrt{D}$ were in
$\mathbb{Q}(x)$, its reduction modulo $\mathfrak{p}$ would be a square
root of $D$ in $\mathbb{F}_{p^f}$; but a nonresidue modulo $p$ remains a
nonsquare in an extension of odd degree. Hence $\sqrt{D} \notin
\mathbb{Q}(x)$, and $H$, being squarefree modulo $p$, is squarefree and
by [](#lem-norm) irreducible: it is the minimal polynomial of $w$ (made
primitive). If no such prime is found among a few candidates, $H$ is
factored. Either way, the factorisation of a polynomial of degree $d^2$
or $2d^2$ is replaced by a polynomial of degree $2d$.

## Operands in a common field {#structured}

The composed polynomial of degree $mn$ splits into many factors over
$\Z$, and modulo every prime, when the Galois group of the compositum is
small, for example when $y \in \mathbb{Q}(x)$ or when $x$ and $y$ are
elements of a cyclotomic field. Then the factorisation, with its
combinatorial recombination of modular factors, is very expensive,
while the result can be computed by linear algebra in a field of degree
about $\max(m, n)$. `_qqbar_binary_op_structured` detects and exploits
three situations, once $mn$ (or $d(d-1)/2$ for conjugates) is at least
36.

*One operand in the field of the other.* Suppose $n \mid m$ and $y \in
\mathbb{Q}(x)$, say $y = R(x)$ with $R \in \mathbb{Q}[t]$, $\deg R < m$.
Then $w = h(x)$ with $h = t \circ R(t) \bmod P$ (using an extended GCD
to invert $R$ modulo $P$ for division), and the minimal polynomial of
$h(x)$ is the minimal polynomial of the matrix of multiplication by
$h(x)$ on $\mathbb{Q}(x)$ with the basis $1, x, \ldots, x^{m-1}$, since
a polynomial annihilates the matrix if and only if it annihilates
$h(x)$. This is what `qqbar_evaluate_fmpq_poly` computes, followed by a
numerical selection of the root. The polynomial $R$ is found with LLL
[@LLL1982] as an integer relation between $1, x, \ldots, x^{m-1}, y$ at
increasing precision, and verified exactly with
`qqbar_equal_fmpq_poly_val`, which checks $Q(R(x)) \equiv 0 \bmod P$.
Since a failed LLL search is costly, a necessary condition is tested
first: for a prime $p$ not dividing the leading coefficients and
discriminants, the degrees of the irreducible factors of $P \bmod p$
are the residue degrees of the primes of $\mathbb{Q}(x)$ above $p$
(Dedekind's theorem), and if $y \in \mathbb{Q}(x)$ then $Q$ must have a
root in each of the corresponding residue fields, two distinct roots if
$P = Q$. This is checked for a few dozen primes.

*A common cyclotomic field.* If $x, y \in \mathbb{Q}(\zeta_\ell)$, every
prime $p \equiv 1 \pmod \ell$ splits completely in
$\mathbb{Q}(\zeta_\ell)$ [@Was1997],
so $P$ and $Q$ split into distinct linear factors modulo $p$. Candidate
conductors $\ell \not\equiv 2 \pmod 4$ with $m \mid \varphi(\ell)$ and
$n \mid \varphi(\ell)$ are enumerated in order of $\varphi(\ell)$ (up to
$\varphi(\ell) \le 128$), and the first that passes this splitting test
for primes $p \equiv 1 \pmod \ell$ (four of them) is taken. Then $x$ and
$y$ are expressed with LLL as polynomials in $\zeta_\ell$, or for real
operands in the basis $1, 2\cos(2\pi k/\ell)$ of the real subfield,
which gives small coefficients for cyclotomic numbers; the expressions
are verified exactly and combined by arithmetic in
$\mathbb{Q}(\zeta_\ell)$.

*Guessing the result.* Alternatively, an integer relation $g$ between
$1, w, \ldots, w^k$ for $k = d, d/2, d/4$ ($d = \max(m, n)$) is sought
with LLL from a numerical value of $w$, and the identity $g(w) = 0$ is
proved by evaluating $g(w)$ at a precision justified by the following
bound.

::: theorem Proving a guessed annihilating polynomial {#thm-guess}
Let $w = x \circ y$, let $g \in \Z[z]$ be nonzero of degree $k$ with
$\|g\|_1 = \sum |g_i|$, let $N = mn$, let $L = \operatorname{lc}(P)^n
\operatorname{lc}(Q)^m$ (with $Q$ replaced by its reversal for division)
and let $R \ge 1$ bound the absolute values of all numbers $a_i \circ
b_j$. If

$$
|g(w)| < |L|^{-k} \bigl(\|g\|_1 R^k\bigr)^{-(N-1)},
$$

then $g(w) = 0$.
:::

::: proof
Suppose $g(w) \ne 0$ and let $f \in \Z[z]$ be the minimal polynomial of
$w$, of degree $e \le N$. Then $g(v) \ne 0$ for every conjugate $v$ of
$w$ (the conjugates are the roots of $f$, and $f$ is irreducible), so
the resultant $\operatorname{Res}(f, g) = \operatorname{lc}(f)^k
\prod_{f(v) = 0} g(v)$ is a nonzero integer and has absolute value at
least 1. The polynomial $L \prod_{i,j} (z - a_i \circ b_j)$ has integer
coefficients (it is a resultant) and vanishes at $w$, so the primitive
polynomial $f$ divides it in $\Z[z]$ by Gauss's lemma, and comparing
leading coefficients, $\operatorname{lc}(f)$ divides $L$. Every
conjugate $v$ is one of the $a_i \circ b_j$, so $|g(v)| \le \|g\|_1
R^k$, and this bound is at least 1. Hence

$$
1 \le |\operatorname{lc}(f)|^k \, |g(w)| \prod_{v \ne w} |g(v)|
\le |L|^k \, |g(w)| \, \bigl(\|g\|_1 R^k\bigr)^{N-1},
$$

which contradicts the hypothesis.
:::

Here $R$ is computed from root bounds for $P$ and $Q$ (`fmpz_poly_bound_roots`):
their sum for addition and subtraction, their product for
multiplication and division (with the bound for the reversed
polynomial). The test needs a precision of about $k \log_2 |L| + (N -
1)(\log_2 \|g\|_1 + k \log_2 R)$ bits; before paying for it, $g(w)$ is
checked at three times the LLL precision, which eliminates most spurious
relations. If the test succeeds, $g$ is factored (it has degree $k$,
much smaller than $N$) and the factor is selected as above.

These LLL searches can fail, and their cost grows with the precision.
The implementation keeps a list of *streams* of attempts, one per
method and lattice dimension, each with a doubling sequence of
precisions, and always runs the cheapest next attempt according to a
cost model for LLL, until the total cost exceeds a budget proportional
to the estimated cost of the fallback (computing and factoring the
composed polynomial). Without detected cyclotomic structure, relations
with small coefficients are less likely and the budget is halved.

## The complete strategy {#strategy}

`qqbar_binary_op` tries, in order:

1. equal operands: $2x$, 0, $x^2$ or 1;
2. a quadratic operand: closed forms or certified norms
   ([](#quadratic));
3. for operands of equal degree at least 4: a rational result, guessed
   numerically with `qqbar_guess` at 128 bits and verified by the exact
   inverse operation (for example $x = t - y$ for a sum $t$), which is
   cheap since one operand is then rational;
4. structured methods within a cost budget ([](#structured));
5. distinct conjugates: the symmetric composed polynomial of half degree
   ([](#composed)), factored;
6. otherwise: the composed polynomial, factored, and the factor selected
   numerically ([](#select)).

Powers $x^n$ are computed by evaluating $t^n$ at $x$ in
$\mathbb{Q}(x)$, after checking for roots of unity, for principal roots
of rational numbers, and for minimal polynomials of the form $P(t) =
H(t^n)$, in which case $x^n$ is a root of $H$. Roots $x^{1/n}$ are
computed by factoring $P(t^n)$ and selecting the principal root;
Capelli's theorem [@Sch2000] often certifies that $P(t^n)$ is
irreducible, so that no factorisation is needed.

## Use in FLINT {#use}

- `qqbar_add`, `qqbar_sub`, `qqbar_mul` and `qqbar_div` handle rational
  operands with `qqbar_scalar_op` and call `qqbar_binary_op` otherwise;
  `qqbar_re`, `qqbar_im` and `qqbar_abs2` combine $x$ with its conjugate
  and benefit from the symmetric composed polynomials.
- `qqbar_fmpz_poly_composed_op` and `qqbar_fmpz_poly_symmetric_composed_op`
  compute composed polynomials from power sums.
- `qqbar_evaluate_fmpq_poly` computes minimal polynomials of elements of
  $\mathbb{Q}(x)$ through multiplication matrices (`nf_elem_rep_mat`
  and `fmpq_mat_minpoly`).
- `qqbar_express_in_field` and `qqbar_guess` expose the LLL searches of
  [](#structured) (see [](algorithms/exact/simplification)).
- `qqbar_evaluate_fmpz_mpoly` and `ca_get_qqbar` evaluate larger
  expressions with these operations, under limits on degrees and sizes
  that bound the cost of the composed polynomials.
