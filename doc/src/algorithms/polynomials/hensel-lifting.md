---
title: Hensel lifting
summary: Lifting a factorisation modulo a prime p into coprime monic factors to a factorisation modulo p^N, for two factors with the Bézout cofactors and for many factors along a factor tree, with quadratic convergence; the correctness proof and the fused tree nodes used by fmpz_poly.
complexity: One quadratic step from $p^e$ to $p^{2e}$ for a factorisation of a degree-$n$ polynomial into $r$ factors costs $O(M(n) \log r)$ operations modulo $p^{2e}$; lifting to $p^N$ costs a constant times the last step.
order: 65
---

Hensel lifting turns a factorisation of $f \in \Z[x]$ modulo a prime
$p$ into one modulo $p^N$ for any $N$. It is the bridge between the
cheap computations modulo a prime, such as factorisation over
$\mathbb{F}_p$, and the integer answer: once $p^N$ exceeds twice a bound
for the coefficients of a factor, the factor is determined by its image
modulo $p^N$ ([](algorithms/polynomials/factoring-integers#bounds)).
The method goes back to Hensel's $p$-adic numbers [@Hen1908];
Zassenhaus introduced quadratic lifting [@Zas1969], and the factor tree
for many factors is described in [@GG2003, Ch. 15].

## Lifting two factors {#two}

The basic step lifts a factorisation into two coprime factors together
with their Bézout cofactors. Below, $m \ge 2$ and $m_1 \ge 1$ are
integers with $m_1 \mid m$ (for example $m = p^e$ and $m_1 = p^{e'}$
with $e' \le e$), and congruences between polynomials hold
coefficientwise.

::: theorem Hensel step {#thm-step}
Let $f, g, h, a, b \in \Z[x]$ with $g$ and $h$ monic, $\deg f = \deg g
+ \deg h$, $\deg a < \deg h$, $\deg b < \deg g$, and

$$
f \equiv gh \pmod m, \qquad ag + bh \equiv 1 \pmod m, \qquad
\operatorname{lc}(f) \equiv 1 \pmod{m m_1}.
$$

Let $e = (f - gh)/m$, and let $\delta$, $\varepsilon$ be the remainders
of $be$ modulo $g$ and of $ae$ modulo $h$, computed modulo $m_1$. Then
$G = g + m\delta$ and $H = h + m\varepsilon$ are monic, congruent to
$g$ and $h$ modulo $m$, and

$$
f \equiv GH \pmod{m m_1}.
$$

With $e_2 = (1 - aG - bH)/m$, the polynomials $A = a + m (a e_2 \bmod
H)$ and $B = b + m (b e_2 \bmod G)$, computed modulo $m_1$, satisfy
$AG + BH \equiv 1 \pmod{m m_1}$.
:::

::: proof
The case $m_1=1$ is immediate, so assume $m_1>1$.
All computations of $\delta$ and $\varepsilon$ are modulo $m_1$, where
$ag + bh \equiv 1$ since $m_1 \mid m$. Write $be = q_1 g + \delta$ and
$ae = q_2 h + \varepsilon$. Then, modulo $m_1$,

$$
e = e(ag + bh) = g(ae) + h(be) = g \varepsilon + h \delta + gh(q_1 + q_2).
$$

The coefficient of $x^{\deg f}$ in $e$ is
$(\operatorname{lc}(f)-1)/m$, which vanishes modulo $m_1$ by the
leading-coefficient hypothesis. Thus $e \bmod m_1$ has degree less
than $\deg f$; also $\deg(g\varepsilon + h\delta) < \deg f$.
Consequently $gh(q_1 + q_2) \equiv e - g\varepsilon - h\delta$ has
degree less than $\deg(gh)$ modulo $m_1$. Since $gh$ is monic,
$q_1 + q_2 \equiv 0 \pmod{m_1}$. So $e \equiv g\varepsilon +
h\delta \pmod{m_1}$, and

$$
GH = gh + m(g\varepsilon + h\delta) + m^2 \delta\varepsilon
\equiv gh + m e = f \pmod{m m_1},
$$

because $m_1 \mid m$ gives $m m_1 \mid m^2$. For the cofactors,
$aG + bH \equiv ag + bh \equiv 1 \pmod{m}$, so $e_2$ is integral, and
$\deg(1 - aG - bH) < \deg f$. The same argument with $G, H$ in place
of $g, h$ gives $e_2 \equiv G (a e_2 \bmod H) + H (b e_2 \bmod G)
\pmod{m_1}$, and therefore $AG + BH = aG + bH + m(e_2 + m_1 \cdot
(\ldots)) \equiv 1 \pmod{m m_1}$.
:::

The leading-coefficient condition is necessary: a product of monic
factors of degrees $\deg g$ and $\deg h$ has leading coefficient 1.
It is automatic when $f$ is monic, and also allows other integer
representatives of a monic polynomial modulo $m m_1$.

::: proposition Uniqueness [@GG2003, Ch. 15] {#prop-unique}
Under the hypotheses of [](#thm-step), let $G',H'$ be monic of the
same degrees and residues modulo $m$ as $g,h$, with
$f\equiv G'H'\pmod{m m_1}$. Then $G'\equiv G$ and $H'\equiv H$
modulo $m m_1$.
:::

::: proof
The case $m_1=1$ follows from the prescribed residues modulo $m$.
Otherwise let $G',H'$ be another such pair, and write
$G'=G+m u$, $H'=H+m v$. Since the factors are monic of the same
degrees, $\deg u<\deg g$ and $\deg v<\deg h$. Subtracting the two
factorisations and dividing by $m$ gives
$$
gv+hu \equiv 0 \pmod{m_1},
$$
because $m_1 \mid m$. In $(\Z/m_1\Z)[x]$, reduction modulo $g$ and
$bh\equiv1\pmod g$ give $u\equiv0\pmod g$. The degree bound on $u$
and monicity of $g$ imply $u=0$ in this coefficient ring. Similarly,
reduction modulo $h$ and $ag\equiv1\pmod h$ give $v=0$. Hence
$G'\equiv G$ and $H'\equiv H$ modulo $m m_1$.
:::

Taking $m_1 = m$ doubles the exponent at each step
(quadratic lifting); taking $m_1 = p$ gains one power of $p$ (linear
lifting). Quadratic lifting reaches $p^N$ in $\lceil \log_2 N \rceil$
steps, and since the cost of a step grows with the size of the
modulus, the total cost is dominated by the last step.

### Normalising the target {#normalisation}

For an initial factorisation $f\equiv\ell gh\pmod m$, where
$\ell=\operatorname{lc}(f)$ is a unit modulo $m m_1$, choose a monic
representative $F\equiv\ell^{-1}f\pmod{m m_1}$. Here $\ell^{-1}$ is
the inverse in $\Z/(m m_1)\Z$. Now $F\equiv gh\pmod m$, so
[](#thm-step) applies to $F$ with the same Bézout cofactors. The
resulting factors satisfy $f\equiv\ell GH\pmod{m m_1}$.

For example, take $m=m_1=5$, $f=6x^2-x$, $g=x$, $h=x-1$, $a=1$ and
$b=-1$. Then $f\equiv gh\pmod5$ and $ag+bh=1$, but the leading
coefficient of $f$ is not 1 modulo 25. Here $e=x^2$, $\delta=0$ and
$\varepsilon=1$, so the correction formulas give $G=x$, $H=x+4$ and
$$
f-GH=5x(x-1)\not\equiv0\pmod{25}.
$$
The leading-coefficient hypothesis therefore cannot be omitted.

Instead, $6^{-1}\equiv21\pmod{25}$ gives the monic target
$F=x^2+4x\equiv21f\pmod{25}$. Now $e=(F-gh)/5=x$, and the same
corrections give $GH=F$ and $6GH\equiv f\pmod{25}$. One choice of
lifted cofactors is $A=6$, $B=-6$, for which
$AG+BH=-24\equiv1\pmod{25}$. The start and continuation routines
perform this modular normalisation before lifting the factor tree.

## Many factors: the factor tree {#tree}

For $f \equiv \ell g_1 \cdots g_r \pmod p$, with
$p\nmid\ell=\operatorname{lc}(f)$ and $r>2$ monic, pairwise coprime
$g_i$, the factors are arranged in a binary tree whose leaves are the
$g_i$ and whose inner nodes are the products of their two children.
At target precision $p^N$, its root represents the monic polynomial
$F\equiv\ell^{-1}f\pmod{p^N}$, with the inverse taken modulo $p^N$.
Each pair of siblings $(v, w)$ carries Bézout cofactors. One lifting
step lifts the root pair, then the pairs below, using the lifted
product at each node as the "$f$" of its children. By [](#thm-step),
the factorisation lifts node by node.
If the tree is balanced, the polynomials at each level have total
degree at most $n$, and a level costs $O(M(n))$ operations, so a step
costs $O(M(n) \log r)$ operations modulo the new modulus [@GG2003, Ch. 15].

::: algorithm Multifactor Hensel lifting {#alg-tree}
**Input:** $f$ with $p \nmid \ell = \operatorname{lc}(f)$, monic pairwise
coprime $g_1, \ldots, g_r$ with $f \equiv \ell \prod g_i \pmod p$, a
target $N \ge 1$.

**Output:** monic $G_i \equiv g_i \pmod p$ with $f \equiv \ell \prod G_i
\pmod{p^N}$.

1. Build the tree: repeatedly multiply the two polynomials of smallest
   degree among the leaves and the products made so far. Compute the
   Bézout cofactors of each pair of siblings modulo $p$ with the
   extended Euclidean algorithm.
2. Replace $f$ by $\ell^{-1} f \bmod p^N$, which is monic.
3. Let $N = e_0 > e_1 > \cdots > e_k = 1$ with $e_{i+1} = \lceil e_i/2
   \rceil$. For $i = k - 1, \ldots, 0$, lift the whole tree from
   $p^{e_{i+1}}$ to $p^{e_i}$ with [](#thm-step) at every node, using
   $m = p^{e_{i+1}}$ and $m_1 = p^{e_i - e_{i+1}}$; in the last step,
   do not lift the cofactors.
4. Return the lifted leaves.
:::

Choosing the two factors of smallest degree, as in Huffman coding,
balances the tree by degree rather than by the number of leaves, which
minimises the total degree of the products that have to be lifted.

## Use in FLINT {#use}

- `fmpz_poly_hensel_build_tree` builds the tree from an
  `nmod_poly_factor_t`, storing the nodes in arrays `v` (factors and
  products), `w` (cofactors) and `link` (the tree structure).
- `_fmpz_poly_hensel_start_lift` makes $f$ monic modulo $p^N$, builds
  the tree and runs the steps of [](#alg-tree) with
  `fmpz_poly_hensel_lift_tree`; it returns the precision of the
  cofactors so that `_fmpz_poly_hensel_continue_lift` can resume the
  lift later, which van Hoeij's algorithm uses to double the precision
  ([](algorithms/polynomials/factoring-integers#van-hoeij)).
  `fmpz_poly_hensel_lift_once` lifts once to a fixed precision, as the
  Zassenhaus algorithm needs.
- `fmpz_poly_hensel_lift`, `fmpz_poly_hensel_lift_without_inverse`
  and `fmpz_poly_hensel_lift_only_inverse` implement [](#thm-step)
  directly, with the reductions modulo $m_1$ done in `fmpz_mod_poly`.
- `fmpz_poly_hensel_lift_tree_recursive` lifts the tree. For nodes
  that are not tiny, it fuses the two halves of [](#thm-step): since
  $G \equiv g$ and $H \equiv h$ modulo $m_1$, the moduli $g \bmod m_1$
  and $h \bmod m_1$ and their Newton preinverses
  ([](algorithms/polynomials/newton)) serve all four divisions, and
  when $p \nmid f(0)$ it computes $H = f/G$ modulo $m m_1$ by exact
  division instead of lifting $h$. The arithmetic runs in `gr` rings
  chosen by the size of the modulus: `nmod` for one word, `mpn_mod`
  for a few words and `fmpz_mod` above; since quadratic lifting walks
  through the moduli $p, p^2, p^4, \ldots$, most levels are cheap even
  when the final precision is large.
- The same technique lifts factorisations in other settings:
  multivariate Hensel lifting in `fmpz_mpoly_factor` and
  `nmod_mpoly_factor` ([](algorithms/multivariate/overview)) and
  $p$-adic root lifting in `fmpz_poly_roots_fmpq`, which lifts roots
  instead of factors by Newton iteration.

::: functions fmpz_poly_hensel_lift_once _fmpz_poly_hensel_start_lift _fmpz_poly_hensel_continue_lift fmpz_poly_hensel_build_tree fmpz_poly_hensel_lift_tree fmpz_poly_hensel_lift_tree_recursive fmpz_poly_hensel_lift :::
