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
with their Bézout cofactors. Below, $m$ and $m_1$ are moduli with $m_1
\mid m$ (for example $m = p^e$ and $m_1 = p^{e'}$ with $e' \le e$), and
congruences between polynomials hold coefficientwise.

::: theorem Hensel step {#thm-step}
Let $f, g, h, a, b \in \Z[x]$ with $g$ and $h$ monic, $\deg f = \deg g
+ \deg h$, $\deg a < \deg h$, $\deg b < \deg g$, and

$$
f \equiv gh \pmod m, \qquad ag + bh \equiv 1 \pmod m.
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
All computations of $\delta$ and $\varepsilon$ are modulo $m_1$, where
$ag + bh \equiv 1$ since $m_1 \mid m$. Write $be = q_1 g + \delta$ and
$ae = q_2 h + \varepsilon$. Then, modulo $m_1$,

$$
e = e(ag + bh) = g(ae) + h(be) = g \varepsilon + h \delta + gh(q_1 + q_2).
$$

Since $f$ and $gh$ have the same degree and $f \equiv gh \pmod m$ with
$gh$ monic, $\deg e < \deg f$; also $\deg(g\varepsilon + h\delta) <
\deg g + \deg h = \deg f$. Hence $gh(q_1 + q_2) \equiv e - g\varepsilon
- h\delta$ has degree less than $\deg(gh)$ modulo $m_1$, and since $gh$
is monic, $q_1 + q_2 \equiv 0 \pmod{m_1}$. So $e \equiv g\varepsilon +
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

The lifted factors are unique: if $f \equiv G'H'$ modulo $m m_1$ with
monic $G' \equiv g$, $H' \equiv h$ modulo $m$, and $g, h$ are coprime
modulo $p$, then $G' = G$ and $H' = H$ modulo $m m_1$ [@GG2003, Ch. 15]. Taking $m_1 = m$ doubles the exponent at each step
(quadratic lifting); taking $m_1 = p$ gains one power of $p$ (linear
lifting). Quadratic lifting reaches $p^N$ in $\lceil \log_2 N \rceil$
steps, and since the cost of a step grows with the size of the
modulus, the total cost is dominated by the last step.

## Many factors: the factor tree {#tree}

For $f \equiv \ell g_1 \cdots g_r \pmod p$ with $r > 2$ monic, pairwise
coprime $g_i$, the factors are arranged in a binary tree whose leaves
are the $g_i$ and whose inner nodes are the products of their two
children; the root is $f / \ell$. Each pair of siblings $(v, w)$
carries Bézout cofactors. One lifting step lifts the root pair, then
the pairs below, using the lifted product at each node as the "$f$" of
its children. By [](#thm-step), the factorisation lifts node by node.
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
