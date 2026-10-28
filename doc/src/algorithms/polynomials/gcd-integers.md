---
title: GCDs of polynomials over Z
summary: The three GCD algorithms of fmpz_poly - the subresultant pseudo-remainder sequence, the heuristic GCD through one integer GCD, and the multimodular algorithm with adaptive termination and certification by multimodular division - with the proofs that make their results correct, and the multimodular extended GCD.
complexity: The subresultant algorithm uses $O(n^2)$ operations on integers of $O(n(b + \log n))$ bits for degree $n$ and $b$-bit coefficients; the multimodular algorithm costs about one modular GCD, $O(M(n) \log n)$ word operations, per prime, with a number of primes proportional to the size of the GCD (plus the verification).
order: 36
---

The GCD of $A, B \in \Z[x]$ is the product of the GCD of their contents
and the GCD of their primitive parts, so it suffices to consider
primitive polynomials. The Euclidean algorithm does not work directly
in $\Z[x]$, and in $\mathbb{Q}[x]$ it is slow, because the
coefficients of the remainders grow exponentially in size unless
contents are removed; even then, intermediate remainders are much
larger than the GCD, which is usually 1. There are three ways out, all
in `fmpz_poly`: pseudo-remainder sequences whose growth is controlled
exactly, evaluation at a large integer, and computation modulo many
primes. Generic versions over other rings are compared in
[](algorithms/polynomials/gcd); over a field the asymptotically fast
method is the [half-GCD](algorithms/polynomials/half-gcd), which
`fmpz_poly_gcd_modular` uses modulo each prime.

Throughout, $A$ and $B$ are primitive with $\deg A = m \ge \deg B = n$,
$G = \gcd(A, B)$ with positive leading coefficient, and $b$ bounds the
bit size of the coefficients.

## The subresultant algorithm {#subresultant}

Pseudo-division ([](fmpz_poly/divrem#pseudo)) gives $\ell^{\delta + 1}
A = QB + R$ with $\delta = \deg A - \deg B$ and $\ell =
\operatorname{lc}(B)$. Repeating it gives a pseudo-remainder sequence
whose coefficients grow exponentially. Collins [@Col1967] and Brown and
Traub [@BT1971] showed that each pseudo-remainder can be divided by a
known factor, and that the results are the *subresultants*,
determinants of submatrices of the Sylvester matrix, whose size is
bounded by Hadamard's inequality.

::: algorithm Subresultant GCD (Cohen's Algorithm 3.3.1) {#alg-subres}
**Input:** primitive $A$, $B$ with $\deg A \ge \deg B \ge 1$.

**Output:** $\gcd(A, B)$ up to the content.

1. $g \gets 1$, $h \gets 1$.
2. $\delta \gets \deg A - \deg B$; $R \gets$ the pseudo-remainder of
   $A$ by $B$ with exponent $\delta + 1$.
3. If $R = 0$, return $\operatorname{pp}(B)$. If $\deg R = 0$, return 1.
4. $A \gets B$, $B \gets R / (g h^\delta)$, $g \gets
   \operatorname{lc}(A)$, $h \gets g^\delta / h^{\delta - 1}$; go to 2.
:::

The divisions in step 4 are exact, and the coefficients of every
remainder are bounded by those of minors of the Sylvester matrix of
size at most $m + n$, that is, by $O((m + n)(b + \log(m + n)))$ bits
[@Coh1996, §3.3], [@Yap2000, Ch. 3]. The algorithm performs $O(mn)$
operations on such integers, which is efficient only for short
polynomials; `_fmpz_poly_gcd` uses it for length less than 6.
`fmpz_poly_resultant_euclidean` uses the same sequence to compute
resultants ([](algorithms/polynomials/resultant-integers)).

## The heuristic GCD {#heuristic}

The heuristic GCD of Char, Geddes and Gonnet [@CGG1989], analysed by
Liao and Fateman [@LF1995], evaluates at a large integer $\xi$ and
replaces the polynomial GCD by one integer GCD.

::: theorem {#thm-gcdheu}
Let $A, B \in \Z[x]$ be primitive and nonzero, and let $\xi \ge 2 + 2
\min(\|A\|_\infty, \|B\|_\infty)$ be an integer. Let $C \in \Z[x]$ be
the polynomial whose coefficients are the digits of $\gamma =
\gcd(A(\xi), B(\xi))$ in base $\xi$ in the symmetric range $[-\xi/2,
\xi/2)$. If $\operatorname{pp}(C)$ divides $A$ and $B$, then
$\operatorname{pp}(C) = \pm\gcd(A, B)$.
:::

::: proof
Let $D = \gcd(A, B)$. Since $\operatorname{pp}(C)$ divides $A$ and $B$,
it divides $D$; write $D = \operatorname{pp}(C) E$ with $E \in \Z[x]$.
By construction $C(\xi) = \gamma$, and $C = c \operatorname{pp}(C)$
with $|c| \le \|C\|_\infty \le \xi/2$. Since $D$ divides $A$ and $B$,
$D(\xi)$ divides $\gamma = c \operatorname{pp}(C)(\xi)$, so $E(\xi)$
divides $c$. Suppose $\deg E \ge 1$, and say $\|A\|_\infty \le
\|B\|_\infty$. Every root $\alpha$ of $E$ is a root of $A$, so by
Cauchy's bound $|\alpha| < 1 + \|A\|_\infty \le \xi/2$, and
$|E(\xi)| \ge \prod_\alpha |\xi - \alpha| > (\xi/2)^{\deg E} \ge
\xi/2 \ge |c|$. Since $c \ne 0$, this contradicts $E(\xi) \mid c$.
Hence $E$ is a constant, and since $D$ and $\operatorname{pp}(C)$ are
primitive, $E = \pm 1$.
:::

`_fmpz_poly_gcd_heuristic` takes $\xi = 2^k$ with $k$ at least 6 more
than the bit size of the smaller of the two coefficient bounds (3
would suffice for the theorem, including a bit for signs) and at least
one more than the larger, packs $A(\xi)$ and $B(\xi)$ into limb arrays
with `_fmpz_poly_bit_pack`, computes their GCD with
`flint_mpn_gcd`, and unpacks the result. The divisibility test is
done on the packed integers, $\gamma \mid A(\xi)$ and $\gamma \mid
B(\xi)$, followed by unpacking the quotients: if the coefficients of
the candidate and of the quotient are small enough that their product
cannot overflow the fields, the integer identity implies the
polynomial identity, and otherwise the product is computed. If the
candidate fails, the function returns 0. Failure is rare for small
coefficients, and the method then costs a few integer operations on
numbers of $(m + n) k$ bits, which is why `_fmpz_poly_gcd` tries it
first when the coefficients of both inputs together have at most 40
bits.

## The multimodular algorithm {#modular}

The modular GCD algorithm of Brown [@Bro1971] computes $\gcd(A \bmod p,
B \bmod p)$ in $\mathbb{F}_p[x]$ for many primes $p$ and reconstructs
$G$ by the Chinese remainder theorem.

::: lemma Good and unlucky primes {#lem-unlucky}
Let $p$ be a prime that divides neither $\operatorname{lc}(A)$ nor
$\operatorname{lc}(B)$, and let $g_p$ be the monic GCD of $A \bmod p$
and $B \bmod p$. Then $\deg g_p \ge \deg G$, with equality if and only
if $p$ does not divide $\operatorname{res}(A/G, B/G)$. In the case of
equality, $g_p = (G \bmod p) / \operatorname{lc}(G)$.
:::

::: proof
$G \bmod p$ divides both reductions and has degree $\deg G$ because
$\operatorname{lc}(G) \mid \operatorname{lc}(A)$, so it divides $g_p$.
The cofactors $A/G$ and $B/G$ keep their degrees modulo $p$, so their
reductions are coprime if and only if their resultant does not vanish
modulo $p$, and $g_p = (G/\operatorname{lc}(G)) \bmod p$ exactly in that
case.
:::

Only finitely many primes are *unlucky*, and they are detected because
they give images of too large a degree. To combine images, they must
be normalised consistently: with $\gamma = \gcd(\operatorname{lc}(A),
\operatorname{lc}(B))$, which is a multiple of $\operatorname{lc}(G)$,
$\gamma g_p$ is the image of $H = (\gamma / \operatorname{lc}(G))\, G
\in \Z[x]$ for every lucky prime. Mignotte's bound
([](algorithms/polynomials/factoring-integers#thm-mignotte)) bounds
the coefficients of $H$, so finitely many primes determine it, and $G =
\operatorname{pp}(H)$.

::: theorem Certification {#thm-certify}
Let $d$ be the least degree of $g_p$ over the primes used, and let
$C \in \Z[x]$ be primitive of degree $d$ with positive leading
coefficient. If $C$ divides $A$ and $B$, then $C = G$.
:::

::: proof
$C$ divides $G$, and $d \ge \deg G$ by [](#lem-unlucky). So $C$ and $G$
have the same degree, $G/C$ is a constant, and since both are
primitive with positive leading coefficient, $C = G$.
:::

The theorem makes the output correct regardless of how the primes are
chosen and when reconstruction is attempted: a wrong candidate is
rejected by the division test, and the algorithm continues with more
primes. `_fmpz_poly_gcd_modular` uses this freedom to adapt the work
to the size of $G$ rather than to an a priori bound:

- **Batches.** $A$ and $B$ are reduced modulo a batch of primes at once
  (with `fmpz_multi_mod_ui` when that is faster than separate
  reductions), and the modular GCDs of a batch, by `_nmod_poly_gcd`,
  are independent and run in parallel. The batch size grows
  geometrically, so at most a constant fraction of the work is wasted.
- **Probes.** Instead of reconstructing $H$ after every batch, only a
  random signed sum $s = \sum_i \pm h_i$ of its coefficients is
  reconstructed, incrementally, which costs one integer CRT. A full
  reconstruction, by a single fast multi-CRT of all stored images, is
  attempted when the symmetric residue of $s$ is 8 bits smaller than
  the modulus, which happens by accident with probability about
  $2^{-7}$; or when the product of the primes exceeds the Mignotte
  bound.
- **Primes.** The first prime is a 20-bit prime unless the
  coefficients are small: the modular GCD is cheapest for small primes,
  and the first prime usually proves that the GCD is 1 (an image of
  degree 0). The following primes are the largest primes below
  $2^{62}$, for which the Euclidean algorithm modulo $p$ can use
  Montgomery arithmetic, or, for $n \ge 2500$, 50-bit primes $p = c
  \cdot 2^{32} + 1$ for which `fft_small` multiplies with a single
  transform, which speeds up the half-GCD.
- **Normalisation.** If $\gcd(A(0), B(0))$ is smaller than $\gamma$,
  the images are normalised by their constant coefficients instead,
  which keeps $H$ small when $A$ and $B$ have a large common factor in
  their leading coefficients.
- **Shortcut.** If an image has degree $\deg B$, the algorithm tests
  directly whether $B$ divides $A$.
- **Verification.** Exact division over $\Z$ is expensive, so the
  quotients $A/C$ and $B/C$ are computed modulo the stored primes, and
  more primes if needed, until they stabilise; then $C \cdot (A/C) = A$
  follows from the congruences if the modulus exceeds twice a bound
  for both sides, and is otherwise checked with one multiplication.

## Extended GCD and cofactors {#xgcd}

Over $\Z$, the Bézout identity holds only up to a constant: if $A$
and $B$ are coprime, there are unique $s, t \in \mathbb{Q}[x]$ with
$\deg s < \deg B$, $\deg t < \deg A$ and $sA + tB = 1$, and by
Cramer's rule applied to the Sylvester system, $rs$ and $rt$ have
integer coefficients for $r = \operatorname{res}(A, B)$.
`fmpz_poly_xgcd_modular` computes $r$
([](algorithms/polynomials/resultant-integers)), then the extended GCDs
modulo primes $p \nmid r \operatorname{lc}(A) \operatorname{lc}(B)$,
scaled by $r \bmod p$, and combines the cofactors by incremental CRT.
When the reconstructed cofactors stop changing, they are checked modulo
each further prime, and the loop stops when the modulus exceeds $2^4$
times a bound for $|sA|$, $|tB|$ and $|r|$; then $sA + tB - r$ is
divisible by a modulus larger than its coefficients, so it is zero.

## Use in FLINT {#use}

- `fmpz_poly_gcd` removes common powers of $x$, then uses
  `fmpz_poly_gcd_subresultant` for length less than 6, tries
  `fmpz_poly_gcd_heuristic` when the coefficients of both inputs have
  at most 40 bits in total, and otherwise, or if the heuristic fails,
  calls `fmpz_poly_gcd_modular`.
- `fmpz_poly_lcm`, `fmpz_poly_squarefree_part`,
  `fmpz_poly_is_squarefree` and `fmpz_poly_factor_squarefree` are built
  on `fmpz_poly_gcd`; `fmpq_poly_gcd` uses it on the numerators.
- `fmpz_poly_xgcd` is `fmpz_poly_xgcd_modular`.
- The multivariate analogues are in `fmpz_mpoly`
  ([](algorithms/multivariate/overview)).

::: functions fmpz_poly_gcd fmpz_poly_gcd_modular fmpz_poly_gcd_heuristic fmpz_poly_gcd_subresultant fmpz_poly_xgcd_modular :::
