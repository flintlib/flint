---
title: GCDs over small finite fields
summary: How the multivariate GCD algorithms over $\mathbb{F}_p$ and $\mathbb{F}_q$ cope with fields too small for random evaluation points, by computing images over extension fields $\mathbb{F}_{p^k}$ or modulo irreducible polynomials of increasing degree; Brown's algorithm over $\mathbb{F}_p$ with pairs of points and early termination; and how `nmod_mpoly_gcd` chooses among Brown's, Zippel's and Hensel-lifting algorithms.
complexity: Over a large field, the cost of the underlying algorithm (Brown, Zippel or Hensel lifting). Over a small field, the same algorithm over $\mathbb{F}_{p^d}$ for $d = O(\log_p D)$, where $D$ bounds the degrees, plus Chinese remaindering in $\mathbb{F}_p[x]$.
order: 45
---

The multivariate GCD algorithms of FLINT reduce the problem to fewer
variables by *evaluation*: a variable $x_j$ is replaced by a value
$\alpha \in K$, the GCD of the images is computed recursively, and the
GCD is recovered from enough images by interpolation (Brown's dense
algorithm [@Bro1971]), by sparse interpolation of the unknown
coefficients of a known support (Zippel's algorithm [@Zip1979]), or by
Hensel lifting of a factorisation $A = G\bar{A}$, $B = G\bar{B}$ of an
image. These algorithms, common to all coefficient rings, are described
in [](algorithms/multivariate/gcd),
[](algorithms/multivariate/sparse-interpolation) and
[](algorithms/multivariate/hensel-lifting). What is specific to finite
fields $K = \mathbb{F}_q$, $q = p^d$, with a word-size characteristic is
that the field may be too small to provide the evaluation points these
algorithms need. FLINT then computes the images over larger
fields, in two ways:

- *medium extensions*: the inputs are mapped into an extension
  $\mathbb{F}_{q^k}$ and the algorithm runs there with evaluation points in
  the extension;
- *reduction modulo irreducible polynomials*: the variable $x_j$ is not
  evaluated at a point but reduced modulo an irreducible polynomial
  $m(x_j)$, so that the image lies over the field $K[x_j]/(m) \cong
  \mathbb{F}_{q^{\deg m}}$, and images are combined by the Chinese
  remainder theorem in $K[x_j]$.

In the sources, versions over $K$ itself are named `smprime`, those over
extensions represented with Zech logarithms `medprime`, and those over
$K[x]/(m)$ with arbitrary $\deg m$ `lgprime` ("small", "medium" and
"large primes", the irreducible polynomials playing the role of primes).
The approach of working modulo prime ideals of higher degree when the
field is too small goes back to the analysis of the modular GCD algorithm
over small fields by Kaltofen and Monagan [@KM1999]. The implementation
is due to Daniel Schultz, whose notes [@Sch-mpoly] describe the
algorithms of the `mpoly` modules.

## Why small fields are a problem {#problem}

Every evaluation algorithm needs points at which the inputs keep their
shape. For Brown's algorithm in the variable $x_n$ with the other
variables $x_1, \ldots, x_{n-1}$, a point $\alpha$ is *bad* if it is a
root of $\gamma = \gcd(\mathrm{lc}(A), \mathrm{lc}(B))$ (the leading
coefficients with respect to $x_1, \ldots, x_{n-1}$, polynomials in
$x_n$), and *unlucky* if the GCD of the images $A(\alpha)$, $B(\alpha)$
has a larger degree than the image of the GCD. Both sets are contained in
the root set of a nonzero polynomial $\Delta \in K[x_n]$: the leading
coefficient for bad points, and a subresultant of the cofactors for
unlucky ones ([@GG2003, §6.7]). Sparse interpolation and Hensel lifting
choose a random point $(\alpha_2, \ldots, \alpha_n)$ for several
variables at once, and have similar conditions expressed by nonzero
polynomials in several variables.

::: lemma Schwartz–Zippel {#lem-sz}
Let $\Delta \in K[x_1, \ldots, x_n]$ be nonzero of total degree $D$ and
let $S \subseteq K$ be finite. For $\alpha$ chosen uniformly at random
from $S^n$, the probability that $\Delta(\alpha) = 0$ is at most
$D/|S|$ [@Sch1980] [@Zip1979].
:::

Over a large field the probability of failure is negligible. Over
$\mathbb{F}_p$ with $p$ small compared with the degrees it is not, and
there may be no good point at all: $\Delta = x^p - x$ vanishes on all of
$\mathbb{F}_p$. Taking points from an extension field restores the
bound of [](#lem-sz) with $|S| = q^k$.

## Extensions do not change the GCD {#extension}

::: lemma {#lem-ext}
Let $L/K$ be an extension of finite fields and $A, B \in K[x_1, \ldots,
x_n]$, not both zero. The monic GCD of $A$ and $B$ in $L[x_1, \ldots,
x_n]$ (with respect to any monomial ordering) has coefficients in $K$ and
is the monic GCD in $K[x_1, \ldots, x_n]$.
:::

::: proof
Let $G$ be the monic GCD over $K$ and write $A = G\bar{A}$, $B =
G\bar{B}$ with $\bar{A}, \bar{B}$ coprime over $K$. Let $H$ be the monic
GCD of $\bar{A}$ and $\bar{B}$ over $L$. The extension $L/K$ is Galois;
an automorphism $\sigma$ of $L$ over $K$, applied to the coefficients, is
a ring automorphism of $L[x_1, \ldots, x_n]$ fixing $\bar{A}$ and
$\bar{B}$, so $\sigma(H)$ is again a monic GCD of $\bar{A}$ and $\bar{B}$,
hence $\sigma(H) = H$. Thus $H \in K[x_1, \ldots, x_n]$. Division by the
monic $H$ with respect to the monomial ordering gives the same quotient
over $K$ and over $L$, so $H$ divides $\bar{A}$ and $\bar{B}$ over $K$ and
$H = 1$. Hence $\bar{A}$ and $\bar{B}$ are coprime over $L$ and $G$ is the
monic GCD over $L$.
:::

So a GCD computed over $\mathbb{F}_{q^k}$ with any algorithm, once made
monic, is the GCD over $\mathbb{F}_q$; the cofactors $A/G$ and $B/G$ also
have coefficients in $\mathbb{F}_q$. This justifies the `medprime`
Hensel-lifting GCD (`nmod_mpolyl_gcd_hensel_medprime`), which maps the
inputs to $\mathbb{F}_{p^k}$ represented by `fq_zech` and lifts there.

## Brown's algorithm over $\mathbb{F}_p$ {#brown-smprime}

Over a field with enough points, Brown's algorithm evaluates the last
variable $x = x_n$ at points $\alpha \in \mathbb{F}_p$, computes the GCDs
of the images in $\mathbb{F}_p[x_1, \ldots, x_{n-1}]$ recursively, and
interpolates the GCD $G$ and the cofactors $\bar A$, $\bar B$ together
[@Bro1971]. The version in FLINT, `nmod_mpolyn_gcd_brown_smprime`
(inputs in the `nmod_mpolyn_t` representation, with coefficients dense
in $x$), follows the notes [@Sch-mpoly] and adds two things to the
textbook algorithm:

- *pairs of points*: for odd $p$ the points are taken in pairs $\pm\alpha$,
  for $\alpha = (p-1)/2 - 1, (p-1)/2 - 2, \ldots, 2$. Writing a coefficient
  as $c(x) = c_0(x^2) + x\, c_1(x^2)$, both values $c(\pm \alpha) = c_0(\alpha^2)
  \pm \alpha c_1(\alpha^2)$ come from one pass over the coefficients, and
  the interpolants are updated modulo $x^2 - \alpha^2$ in one step. For
  $p = 2$ the function fails, and the caller uses extensions;
- *early termination*: an interpolant that does not change when two new
  images are added has probably reached its final value, and the
  algorithm tests it by division before the degree bound is reached.

As in [](#lgprime), let $\gamma = \gcd(\mathrm{lc}(A), \mathrm{lc}(B)) \in
\mathbb{F}_p[x]$ for the leading coefficients with respect to $x_1,
\ldots, x_{n-1}$ in lex order, let leading monomials be taken with
respect to these variables, and let $\operatorname{pp}$ denote the
primitive part over $\mathbb{F}_p[x]$.

::: algorithm Brown's GCD over $\mathbb{F}_p$ with early termination {#alg-brown-sm}
Input: nonzero $A, B \in \mathbb{F}_p[x_1, \ldots, x_{n-1}][x]$, $p$ odd;
a limit $\lambda$ for early termination.
Output: $G = \gcd(A, B)$, $\bar A = A/G$, $\bar B = B/G$, or failure.

1. Remove the contents of $A$ and $B$ in $\mathbb{F}_p[x]$ (their GCD and
   cofactors are multiplied back at the end). Set $\beta = 1 + \deg \gamma +
   \max(\deg_x A, \deg_x B)$ and $M = 1$.
2. Take the next point $\alpha$; fail if there is none. If
   $\gamma(\alpha)\gamma(-\alpha) = 0$, repeat this step.
3. Compute $g_\pm = \gcd(A(\pm\alpha), B(\pm\alpha))$ and the cofactors
   recursively (with the Euclidean algorithm in $\mathbb{F}_p[x_1]$ at the
   bottom); on failure go to step 2. If $g_+$ or $g_-$ is constant, return
   $G = 1$. If their leading monomials differ, go to step 2. Scale them to
   have leading coefficients $\gamma(\pm\alpha)$.
4. If $M \ne 1$: if the leading monomial of $g_\pm$ is larger than that of
   $G$, go to step 2; if smaller, set $M = 1$.
5. Combine $G$ with $g_\pm$ modulo $M (x^2 - \alpha^2)$. If $M \ne 1$, $G$
   did not change and $\deg M < \lambda$: let $T = \operatorname{pp}(G)$;
   if $T \mid A$ and $T \mid B$, return $T$ (made monic) and the quotients.
   Do the same for $\bar A$ (if $T = \operatorname{pp}(\bar A)$ divides $A$
   and $A/T$ divides $B$, return $G = A/T$) and for $\bar B$. Set $M \gets M
   (x^2 - \alpha^2)$.
6. If $\deg M < \beta$, go to step 2.
7. If $\deg \gamma + \deg_x A = \deg_x G + \deg_x \bar A$ and $\deg \gamma +
   \deg_x B = \deg_x G + \deg_x \bar B$, return $\operatorname{pp}(G)$ and
   the cofactors divided by its leading coefficient. Otherwise set $M = 1$
   and go to step 2.
:::

::: theorem {#thm-brown-sm}
A result returned by [](#alg-brown-sm) is correct.
:::

::: proof
For step 7 the proof of [](#thm-brown-lg) applies with $m = (x -
\alpha)(x + \alpha)$. For step 5, let $H = \gcd(A, B)$, primitive over
$\mathbb{F}_p[x]$ because the contents were removed. Every image has
leading monomial at least that of $H$, as in the proof of
[](#thm-brown-lg), so the same holds for $G$ and $T$, and the leading
monomial of an image of $\bar A$ is at most $\mathrm{lm}(A) -
\mathrm{lm}(H)$, the leading monomial of $A/H$. If $T$ divides $A$ and
$B$, then $T \mid H$, so $\mathrm{lm}(T) \le \mathrm{lm}(H)$; hence the
leading monomials are equal, $H/T$ lies in $\mathbb{F}_p[x]$, and it is a
unit because $H$ is primitive. If $T = \operatorname{pp}(\bar A)$ divides
$A$ and $G' = A/T$ divides $B$, then $G'$ is a common divisor with
$\mathrm{lm}(G') = \mathrm{lm}(A) - \mathrm{lm}(T) \ge \mathrm{lm}(H)$,
and $G' \mid H$; as before $G'$ and $H$ differ by a unit, using that
$G'$, a divisor of the primitive $A$, is primitive.
:::

The early tests cost a division each and are only made while $\deg M <
\lambda$; `mpoly_gcd_info_get_brown_upper_limit` sets $\lambda$ to at
most half of $\beta$, smaller for dense inputs, and to 0 (no early
termination) when no reliable degree bounds are known. In the bivariate
case, `nmod_mpolyn_gcd_brown_smprime_bivar` does the same with univariate
images, and `nmod_mpolyn_gcd_brown_smprime_threaded_pool` distributes the
points over threads. The version over $\Z$,
[](algorithms/multivariate/gcd#alg-brown), combines images modulo
word-size primes computed by this algorithm and replaces the
divisibility tests by a bound on the coefficients.

## Reduction modulo irreducible polynomials {#lgprime}

Instead of mapping all coefficients into an extension, Brown's and
Zippel's algorithms can treat one variable $x = x_n$ specially. For an
irreducible $m \in K[x]$ of degree $e$, the reduction
$$
\pi_m : K[x_1, \ldots, x_{n-1}][x] \to \mathbb{F}_m[x_1, \ldots, x_{n-1}],
\qquad \mathbb{F}_m = K[x]/(m) \cong \mathbb{F}_{q^e},
$$
replaces the evaluation $x = \alpha$, which is the case $m = x - \alpha$.
All the properties of evaluation used by Brown's algorithm hold for
$\pi_m$ with "$\Delta(\alpha) \ne 0$" replaced by "$m \nmid \Delta$", and
the Chinese remainder theorem in the principal ideal domain $K[x]$
replaces Lagrange interpolation: from images modulo pairwise coprime
$m_1, \ldots, m_r$ it recovers a polynomial with coefficients in $K[x]$
of degree less than $\sum_i \deg m_i$.

::: lemma Few unlucky moduli {#lem-unlucky}
Let $\Delta \in K[x]$ be nonzero of degree $D$. Distinct monic
irreducible polynomials $m_1, m_2, \ldots$ that divide $\Delta$ satisfy
$\sum_i \deg m_i \le D$. In particular at most $D/e$ of them have degree
$e$, while the number $I_e$ of monic irreducible polynomials of degree $e$
over $\mathbb{F}_q$ satisfies $e I_e \ge q^e - 2q^{e/2}$.
:::

::: proof
The $m_i$ are pairwise coprime divisors of $\Delta$, so their product
divides $\Delta$. For the count, Gauss's formula $e I_e = \sum_{c \mid e}
\mu(c) q^{e/c}$ gives $e I_e \ge q^e - \sum_{j=1}^{\lfloor e/2 \rfloor}
q^j \ge q^e - q^{e/2+1}/(q-1) \ge q^e - 2q^{e/2}$.
:::

So whatever the size of $q$, moduli of degree $e \approx \log_q D$
suffice: a random irreducible polynomial of degree $e$ divides $\Delta$
with probability at most $D/(q^e - 2q^{e/2})$. FLINT does not pick
the moduli at random: it uses the irreducible polynomial chosen by
`fq_nmod_ctx_init_ui` for $\mathbb{F}_{p^e}$ (a Conway polynomial or a
fixed polynomial of minimal weight) for $e = e_0, e_0 + 1, e_0 + 2,
\ldots$, starting from $e_0 \approx \max(2, \lfloor 20/\beta \rfloor)$
where $\beta$ is the bit length of $p$, so that $p^{e_0}$ has about 20
bits or more. Moduli of distinct degrees are automatically pairwise
coprime, and by [](#lem-unlucky) the moduli that are bad or unlucky have
total degree at most $\deg \Delta$, so the process terminates. Since the
degrees grow, after $r$ images the modulus has degree about $r e_0 +
r^2/2$, and $O(\sqrt{D})$ images suffice to reach a degree bound $D$.

The resulting algorithm for a GCD in $K[x_1, \ldots, x_{n-1}][x]$ is
Brown's algorithm [@Bro1971] with these "primes". Let $\gamma = \gcd(
\mathrm{lc}(A), \mathrm{lc}(B)) \in K[x]$, where the leading coefficients
are taken with respect to $x_1, \ldots, x_{n-1}$ in lexicographic order,
and assume that the contents of $A$ and $B$ with respect to $x_1, \ldots,
x_{n-1}$ have been removed.

::: algorithm Brown's GCD modulo irreducible polynomials {#alg-brown-lg}
Input: $A, B \in K[x_1, \ldots, x_{n-1}][x]$, primitive over $K[x]$.
Output: $G = \gcd(A, B)$ and the cofactors $\bar{A}, \bar{B}$.

1. Set $\beta = 1 + \deg \gamma + \max(\deg_x A, \deg_x B)$, $M = 1$ and
   $e = e_0$.
2. Let $m$ be the irreducible polynomial of degree $e$, and $e \gets e + 1$.
   If $m \mid \gamma$, repeat this step.
3. Compute $g_m = \gcd(\pi_m(A), \pi_m(B))$ over $\mathbb{F}_m$ (by
   evaluation of the next variable at points of $\mathbb{F}_m$, which is
   now large), scaled to have leading coefficient $\pi_m(\gamma)$, and
   the cofactors $\bar{a}_m$, $\bar{b}_m$.
4. If $g_m = 1$, return $G = 1$. If $M \ne 1$ and the leading monomial of
   $g_m$ is larger than that of the current $G$, the modulus is unlucky:
   go to step 2. If it is smaller, all previous moduli were unlucky: set
   $M = 1$.
5. Combine: $G \equiv g_m$, $\bar{A} \equiv \bar{a}_m$, $\bar{B} \equiv
   \bar{b}_m \pmod m$ and $\equiv$ the previous values modulo $M$ (Chinese
   remainder theorem in $K[x]$); $M \gets M m$.
6. If $\deg M < \beta$, go to step 2.
7. If $\deg \gamma + \deg_x A = \deg_x G + \deg_x \bar{A}$ and $\deg \gamma
   + \deg_x B = \deg_x G + \deg_x \bar{B}$, return the primitive part of $G$
   and the corresponding cofactors. Otherwise go to step 2.
:::

::: theorem {#thm-brown-lg}
When [](#alg-brown-lg) returns in step 7, its output is correct.
:::

::: proof
Every image satisfies $g_m \bar{a}_m = \pi_m(\gamma A)$ and $g_m \bar{b}_m
= \pi_m(\gamma B)$, so $G\bar{A} \equiv \gamma A$ and $G \bar{B} \equiv
\gamma B \pmod M$, the congruences holding coefficientwise in $K[x]$. By
the test of step 7, $\deg_x(G \bar{A}) \le \deg_x G + \deg_x \bar{A} =
\deg \gamma + \deg_x A < \beta \le \deg M$, and likewise for $\gamma A$,
so the congruence is an equality: $G \bar{A} = \gamma A$, and similarly
$G\bar{B} = \gamma B$. Thus $G$ is a common divisor of $\gamma A$ and
$\gamma B$, and its primitive part $G'$ divides $A$ and $B$ by Gauss's
lemma, so $G'$ divides the true GCD $H$. Leading monomials are taken with
respect to $x_1, \ldots, x_{n-1}$. Since $\mathrm{lc}(H)$ divides $\gamma$
and $m \nmid \gamma$, $\pi_m(H)$ has the leading monomial of $H$, and it
divides $g_m$; hence the leading monomial of every image, and so of $G$
and $G'$, is at least that of $H$. As $G'$ divides $H$, the quotient
$H/G'$ has leading monomial 1, that is, lies in $K[x]$, and is a unit
because $H$ is primitive.
:::

Termination follows from [](#lem-unlucky): only finitely many moduli
are bad or unlucky, and after them every image has the correct leading
monomial, so $G$ converges.

Zippel's algorithm is adapted in the same way. `nmod_mpolyu_gcdm_zippel`
first tries the algorithm over $\mathbb{F}_p$ (`nmod_mpolyu_gcdp_zippel`),
and if that fails, reduces the last variable modulo irreducible
polynomials of increasing degree, computes each image over
$\mathbb{F}_{p^e}$ with Zippel's algorithm (the first image fixes the
monomial structure, the later ones are found by solving for the
coefficients), and combines the images by Chinese remaindering in
$\mathbb{F}_p[x]$. The two-main-variable variant
`nmod_mpolyl_gcd_zippel_lgprime` (sparse interpolation with two main
variables and a known leading coefficient, as in [@HuMon2016]) reduces
the last variable modulo $m$ and evaluates the others at random points of
$\mathbb{F}_m$.

For `fq_nmod_mpoly` the base field is already an extension
$\mathbb{F}_q = \mathbb{F}_p[\theta]/(f)$, and the reduction is modulo an
irreducible $h \in \mathbb{F}_q[x]$ of degree $n$. FLINT presents
$\mathbb{F}_q[x]/(h)$ as $\mathbb{F}_p[\varphi]/(g)$ with $\deg g = n
\deg f$, and an embedding records the images of $\theta$ and $x$; since
$f$ has $\deg f$ roots in the large field, there are $\deg f$ embeddings
for each $g$, which are tried in turn before $n$ is increased
(`bad_fq_nmod_mpoly_embed_chooser_init` and
`bad_fq_nmod_mpoly_embed_chooser_next`).

## Choosing an algorithm {#choice}

`nmod_mpoly_gcd` and `nmod_mpoly_gcd_cofactors` (and, for a word-size
modulus, `fmpz_mod_mpoly_gcd` after converting to `nmod_mpoly`) proceed
as follows; `fq_nmod_mpoly_gcd` is analogous.

1. If an input is zero or a monomial, or if the inputs have the same
   shape after removing the term contents (the cofactors may be
   monomials), the GCD is found directly. Inputs whose exponents need more
   than one word are repacked, after deflation if necessary.
2. The exponents are deflated: the common powers of each variable are
   removed and the exponents divided by the GCD of their differences
   (`nmod_mpoly_deflation`). If the essential parts of $A$ and $B$ share
   no variable, the GCD is a monomial; if they share exactly one variable
   and depend on no other, the computation is univariate
   (`nmod_poly_gcd`); if a variable appears in only one of them, the GCD is
   the GCD of the other polynomial with the coefficients of that one with
   respect to the variable.
3. Degree bounds for the GCD are estimated by evaluating all variables but
   one at random points, first over $\mathbb{F}_p$, then (if the bounds
   are not reliable, or $p < 8 - m$ for $m$ variables) over
   $\mathbb{F}_{p^k}$ with Zech logarithms, then over larger extensions.
   If the bounds show that the GCD is 1, or that it equals $A$ or $B$ (a
   divisibility test), the computation stops.
4. With two variables, Brown's algorithm is tried first if the inputs are
   dense (density above 0.05, or above 0.2 when the degree in the second
   variable exceeds $p/2$), otherwise Hensel lifting first.
5. With more variables: for dense inputs (density above 0.08) Hensel
   lifting is tried first when some degree exceeds $p/2$, then Brown's
   algorithm; then Hensel lifting if the density (relative to the total
   degree) exceeds 0.05 (0.001 when a degree exceeds $p/2$); then the two
   sparse interpolation algorithms, Zippel's first if the inputs are short
   compared with the degree in the second main variable; finally Hensel
   lifting and Brown's algorithm.

A degree above $p/2$ in an evaluated variable means that $\mathbb{F}_p$
does not have enough points for dense interpolation; Brown's algorithm
would then work modulo irreducible polynomials, while Hensel lifting
needs only one good evaluation point, possibly in an extension. Each
algorithm returns failure rather than an incorrect result, and the next
one is tried. The returned GCD is monic, and the cofactors are scaled
accordingly.

## Use in FLINT {#use}

- `nmod_mpoly_gcd`, `nmod_mpoly_gcd_cofactors`: the selection of
  [](#choice), through `_nmod_mpoly_gcd_algo`.
- `nmod_mpoly_gcd_brown`, `nmod_mpoly_gcd_zippel`,
  `nmod_mpoly_gcd_zippel2`, `nmod_mpoly_gcd_hensel`: one algorithm each,
  falling back to extensions over small fields.
- Brown modulo irreducible polynomials: `nmod_mpolyn_gcd_brown_lgprime`
  ([](#alg-brown-lg)), with the interpolation helpers
  `nmod_mpolyn_interp_reduce_lg_mpolyn`, `nmod_mpolyn_interp_lift_lg_mpolyn`
  and `nmod_mpolyn_interp_crt_lg_mpolyn`; over $\mathbb{F}_p$ itself
  `nmod_mpolyn_gcd_brown_smprime` ([](#alg-brown-sm)), with the
  interpolation helpers `nmod_mpolyn_interp_reduce_2sm_mpolyn` and
  `nmod_mpolyn_interp_crt_2sm_mpolyn`.
- Zippel: `nmod_mpolyu_gcdm_zippel` (with fallback),
  `nmod_mpolyu_gcdp_zippel` (over $\mathbb{F}_p$),
  `nmod_mpolyl_gcd_zippel_smprime` and `nmod_mpolyl_gcd_zippel_lgprime`
  (two main variables).
- Hensel lifting: `nmod_mpolyl_gcd_hensel_smprime`, and
  `nmod_mpolyl_gcd_hensel_medprime` over $\mathbb{F}_{p^k}$ with $p^k \le
  10^6$ in `fq_zech` representation ([](#lem-ext)).
- Over $\mathbb{F}_q$: `fq_nmod_mpoly_gcd` with
  `fq_nmod_mpolyn_gcd_brown_lgprime`, `fq_nmod_mpolyu_gcdm_zippel` and
  `fq_nmod_mpolyl_gcd_zippel_lgprime` using the embeddings of
  `fq_nmod_mpoly.h`; `fq_zech_mpoly_gcd` converts to `fq_nmod_mpoly`.
- `fmpz_mod_mpoly_gcd` converts to `nmod_mpoly` when the modulus fits a
  word; otherwise $p > 2^{\bits}$ and only the algorithms over
  $\mathbb{F}_p$ are needed.

::: functions nmod_mpoly_gcd nmod_mpoly_gcd_cofactors nmod_mpoly_gcd_brown nmod_mpoly_gcd_zippel nmod_mpoly_gcd_zippel2 nmod_mpoly_gcd_hensel fq_nmod_mpoly_gcd fmpz_mod_mpoly_gcd :::
