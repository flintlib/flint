---
title: Half-GCD
summary: The half-gcd algorithm computes the quotients of the first half of the Euclidean remainder sequence from the top halves of the inputs, recursively, giving GCDs, extended GCDs and resultants in $O(M(n) \log n)$ operations; with a correctness proof and FLINT's cutoffs for `nmod_poly`.
complexity: $O(M(n) \log n)$ field operations for the half-gcd, GCD, extended GCD and resultant of polynomials of degree $n$, against $O(n^2)$ for the Euclidean algorithm.
order: 40
---

The Euclidean algorithm computes $\gcd(a, b)$ for polynomials of degree
at most $d$ with $O(d^2)$ operations in the coefficient field. The
half-gcd algorithm reduces this to $O(M(d) \log d)$, where $M(d)$ is the
cost of a multiplication. It rests on the observation that the first
quotients of the Euclidean algorithm depend only on the leading
coefficients of the inputs, so that about half of the quotients can be
computed recursively from the top halves of $a$ and $b$. The idea is due
to Knuth [@Knu1970] and Schönhage [@Sch1971] for integers and continued
fractions, and to Moenck [@Moe1973] for polynomials; a correct and
complete version for polynomials was given by Brent, Gustavson and Yun
[@BGY1980], and Thull and Yap [@ThullYap1990] gave a unified treatment
for polynomials and integers. A textbook account is
[@GG2003, Ch. 11]. Thomé [@Tho2002] generalised
the recursion to matrix polynomials.

Throughout, the coefficients lie in a field $K$; for `nmod_poly` this
means that the modulus must be prime (more precisely, every leading
coefficient that occurs must be invertible).

## Generic version {#generic}

The half-gcd algorithm is implemented once, in `gr_poly`, for any field
in which the leading coefficients that occur can be inverted; a ring
$\Z/N\Z$ pretending to be a field (`gr_ctx_set_is_pretend_field`) runs
it as well, and a non-invertible leading coefficient then reveals a
factor of $N$. The generic GCD, extended GCD and resultant use it over
finite fields from length 200, with inner cutoff 100, and the
Euclidean algorithm over other fields, where the half-gcd algorithm
suffers from coefficient growth (exact infinite fields such as $\mathbb{Q}$) or
numerical instability (floating-point and ball arithmetic).
The type-specific modules call the same functions with their own
cutoffs ([](algorithms/polynomials/gcd#specialisations)); those of
`nmod_poly` are listed in [](#use). The $2 \times 2$ matrix products are the ring method
`_gr_poly_hgcd_mat_mul`, so a ring can replace the products described
below. [](algorithms/polynomials/gcd) compares the half-gcd algorithm
with the Euclidean and subresultant algorithms over all rings.

## The remainder sequence and its matrices {#sequence}

For $a, b \in K[x]$ with $\deg a > \deg b$, the Euclidean algorithm
produces remainders $r_0 = a$, $r_1 = b$ and

$$
r_{i-1} = q_i r_i + r_{i+1}, \qquad \deg r_{i+1} < \deg r_i,
$$

until some $r_{k+1} = 0$; then $r_k$ is a GCD. With the quotient matrices
$Q_i = \begin{pmatrix} q_i & 1 \\ 1 & 0 \end{pmatrix}$ we have
$(r_{i-1}, r_i)^T = Q_i (r_i, r_{i+1})^T$, so

$$
\begin{pmatrix} a \\ b \end{pmatrix} = M_j \begin{pmatrix} r_j \\ r_{j+1} \end{pmatrix},
\qquad M_j = Q_1 Q_2 \cdots Q_j, \qquad \det M_j = (-1)^j.
$$

Since $\deg r_{i-1} = \deg q_i + \deg r_i$, every $q_i$ has degree at
least 1, and the entries of $M_j$ have degree at most
$\sum_{i \le j} \deg q_i = \deg a - \deg r_j$.

::: lemma Uniqueness {#lem-unique}
Suppose that $(a, b)^T = Q_1 \cdots Q_j (A, B)^T$ for polynomials $q_1, \ldots, q_j$ of
degree at least 1, with $\deg A > \deg B$. Then $q_1, \ldots, q_j$ are the
first $j$ Euclidean quotients of $(a, b)$, and $A = r_j$, $B = r_{j+1}$.
:::

::: proof
Let $(u_{i-1}, u_i)^T = Q_i \cdots Q_j (A, B)^T$, so that $u_j = A$,
$u_{j+1} = B$ and $u_{i-1} = q_i u_i + u_{i+1}$. By downward induction,
$\deg u_i > \deg u_{i+1}$ for all $i$: this holds for $i = j$, and then
$\deg u_{i-1} = \deg q_i + \deg u_i > \deg u_i$. Hence each $u_{i+1}$ is
the remainder and $q_i$ the quotient of the division of $u_{i-1}$ by
$u_i$, which are unique. By induction from $i = 1$, the $u_i$ are the
Euclidean remainders.
:::

## The half-gcd {#hgcd}

::: definition {#def-hgcd}
Let $\deg a = d > \deg b$. The *half-gcd* of $(a, b)$ is the triple
$(M_j, r_j, r_{j+1})$ for the unique $j$ with

$$
\deg r_j \ge \lceil d/2 \rceil > \deg r_{j+1},
$$

that is $2 \deg r_j \ge d > 2 \deg r_{j+1}$.
:::

The index $j$ exists and is unique because the degrees of the remainders
decrease strictly from $\deg r_0 = d$ (and $\deg 0 = -\infty$). The
entries of $M_j$ have degree at most $d - \deg r_j \le \lfloor d/2 \rfloor$.

The key fact is that the half-gcd of the top halves of $a$ and $b$ is
an initial part of the half-gcd of $a$ and $b$.

::: lemma Lifting {#lem-lift}
Let $\deg a = d > \deg b \ge m$, write $a = x^m a_1 + a_0$ and
$b = x^m b_1 + b_0$ with $\deg a_0, \deg b_0 < m$, and let $d' = d - m$
and $(M, A^*, B^*)$ be the half-gcd of $(a_1, b_1)$. Then
$(A, B)^T = M^{-1} (a, b)^T$ satisfies

$$
(A, B) = x^m (A^*, B^*) + M^{-1} (a_0, b_0), \qquad
\deg A = m + \deg A^* \ge m + \lceil d'/2 \rceil > \deg B,
$$

and $A$, $B$ are consecutive remainders of $(a, b)$, with
$(a, b)^T = M (A, B)^T$.
:::

::: proof
The first identity is linearity. $M^{-1} = \pm \operatorname{adj}(M)$
has entries of degree at most $d' - \deg A^* \le \lfloor d'/2 \rfloor \le
\deg A^*$, so $M^{-1}(a_0, b_0)$ has degree less than $m + \deg A^*$.
Hence $\deg A = m + \deg A^*$. Also $\deg B^* < \lceil d'/2 \rceil \le
\deg A^*$, so $\deg B < m + \deg A^* = \deg A$. Since $M$ is a product of
quotient matrices with quotients of degree at least 1 (those of
$(a_1, b_1)$, where $\deg a_1 = d' > \deg b_1$), [](#lem-unique) applies.
:::

::: algorithm Half-gcd {#alg-hgcd}
**Input:** $a$, $b$ with $\deg a = d > \deg b$.

**Output:** the half-gcd $(M, A, B)$ of $(a, b)$.

1. $m \gets \lceil d/2 \rceil$. If $\deg b < m$, return $(I, a, b)$.
2. If $d$ is below a cutoff, run the Euclidean algorithm, multiplying up
   the quotient matrices, until the remainder degree drops below $m$, and
   return.
3. $(R, \cdot, \cdot) \gets$ half-gcd of $(a \operatorname{div} x^m,
   b \operatorname{div} x^m)$; $(A, B) \gets R^{-1}(a, b)$, computed as in
   [](#lem-lift) from the outputs of the recursive call and the low parts
   $a_0$, $b_0$.
4. If $\deg B < m$, return $(R, A, B)$.
5. Divide: $A = qB + D$. Let $k = 2m - \deg B$.
6. $(S, \cdot, \cdot) \gets$ half-gcd of $(B \operatorname{div} x^k,
   D \operatorname{div} x^k)$; $(A', B') \gets S^{-1}(B, D)$, again via
   [](#lem-lift).
7. Return $(R \, Q \, S, A', B')$ with $Q = \begin{pmatrix} q & 1 \\ 1 & 0 \end{pmatrix}$.
:::

::: theorem {#thm-hgcd}
[](#alg-hgcd) returns the half-gcd of $(a, b)$.
:::

::: proof
By induction on $d$. Step 2 is the definition. In step 3, with
$d' = d - m = \lfloor d/2 \rfloor$, [](#lem-lift) shows that $A, B$ are
consecutive remainders with $(a, b)^T = R (A, B)^T$ and
$\deg A \ge m + \lceil d'/2 \rceil \ge m$. If $\deg B < m$, they straddle
$m = \lceil d/2 \rceil$ as required.

Otherwise $m \le \deg B < m + \lceil d'/2 \rceil$, so $0 < k \le m$, and
$B, D$ are the next pair of remainders. The polynomial $B
\operatorname{div} x^k$ has degree $d'' = \deg B - k = 2 \deg B - 2m$,
and $\lceil d''/2 \rceil = \deg B - m$. Applying [](#lem-lift) to
$(B, D)$ with $k$ in place of $m$ (if $\deg D < k$ the recursive call
returns the identity and $(A', B') = (B, D)$), $A'$ and $B'$ are
consecutive remainders of $(B, D)$, hence of $(a, b)$, with

$$
\deg A' \ge k + \lceil d''/2 \rceil = m > \deg B'.
$$

The matrix is $RQS$ because $(a, b)^T = R(A, B)^T = RQ(B, D)^T =
RQS(A', B')^T$.
:::

::: proposition Complexity {#prop-hgcd-cost}
[](#alg-hgcd) uses $O(M(d) \log d)$ operations in $K$.
:::

::: proof
Both recursive calls receive inputs of degree at most about $d/2$
($d' = \lfloor d/2 \rfloor$, and $d'' < 2 \lceil d'/2 \rceil$). The other
steps are a constant number of products of polynomials of degree
$O(d)$ (applying $R^{-1}$ and $S^{-1}$, whose entries have degree at
most $d/4$ or so, to vectors of degree $O(d)$, one division, and one
$2 \times 2$ matrix product), which cost $O(M(d))$. Hence
$H(d) \le 2H(d/2) + O(M(d))$, which gives $O(M(d) \log d)$ since
$M(d)/d$ is nondecreasing.
:::

## GCD, extended GCD and resultant {#applications}

**GCD.** Starting from $(a, b)$, one half-gcd step produces consecutive
remainders $(A, B)$ with $\deg B < d/2$. After one division the degree
has halved, and repeating gives the GCD. The total cost is dominated by
the first call, so it is also $O(M(d) \log d)$. Only the remainders are
needed, not the matrix, so the top-level call can skip the final matrix
product.

**Extended GCD.** The cofactors $s$, $t$ with $sa + tb = \gcd(a, b)$ are
read off from the product of all the matrices: if
$(a, b)^T = M (r_k, 0)^T$, then $M^{-1} = \pm \operatorname{adj}(M)$ gives
$r_k$ as a combination of $a$ and $b$.

**Resultant.** For $a = qb + r$ with $\deg a = \alpha$, $\deg b = \beta$,
$\deg r = \rho$ and leading coefficient $c$ of $b$,

::: lemma {#lem-res}
$\operatorname{Res}(a, b) = (-1)^{\alpha\beta} c^{\alpha - \rho}
\operatorname{Res}(b, r)$ if $r \ne 0$, $\operatorname{Res}(a, b) = c^\alpha$
if $\beta = 0$, and $\operatorname{Res}(a, b) = 0$ if $r = 0$ and $\beta > 0$.
:::

::: proof
$\operatorname{Res}(b, a) = c^\alpha \prod_{b(\xi) = 0} a(\xi)$ over the
roots $\xi$ of $b$ in an algebraic closure, with multiplicity. Since
$a(\xi) = r(\xi)$, this equals $c^{\alpha - \rho} \operatorname{Res}(b, r)$,
and $\operatorname{Res}(a, b) = (-1)^{\alpha\beta} \operatorname{Res}(b, a)$.
The other cases follow from the same product formula.
:::

So the resultant is the product, over the steps of the Euclidean
algorithm, of factors that depend only on the degrees and leading
coefficients of the remainders [@GG2003, Ch. 6 and §11.2]. The half-gcd algorithm
performs the same divisions as the Euclidean algorithm, except that in
the recursive calls the remainders are truncated: they have the correct
leading coefficients, and their degrees are offset by the number of
dropped coefficients. FLINT accumulates the factors of [](#lem-res) during
the half-gcd, tracking that offset. The last remainder produced inside a
recursive call is not yet a true remainder (its low coefficients are
corrected only by the lifting step of [](#lem-lift)), so its factor is
deferred until its true degree is known. Fast computation of
resultants and subresultants along these lines is due to Reischert
[@Rei1997]; see also [@GG2003, §11.2].

## Use in FLINT {#use}

The algorithm is implemented once, in `gr_poly` (`_gr_poly_hgcd`,
`_gr_poly_gcd_hgcd`, `_gr_poly_xgcd_hgcd`, `_gr_poly_resultant_hgcd`),
and `nmod_poly` calls it with the `gr` context for $\Z/n\Z$ and cutoffs
tabulated by the bit length $\beta$ of the modulus:

| Cutoff | Function | Range of values | Meaning |
|---|---|---|---|
| inner | `nmod_poly_hgcd_iter_recursive_cutoff` | 60–244 | step 2 of [](#alg-hgcd): Euclid inside the recursion |
| outer | `nmod_poly_hgcd_outer_cutoff` | 310–1585 | in the GCD and resultant loops, finish with Euclid below this length |
| GCD | `nmod_poly_gcd_hgcd_cutoff` | 470–2199 | `nmod_poly_gcd` and `nmod_poly_resultant` use the half-gcd from this length of the shorter input |
| extended GCD | `nmod_poly_xgcd_hgcd_cutoff` | 72–388 | `nmod_poly_xgcd` uses the half-gcd from this length |

The extended GCD switches earlier than the GCD, since the Euclidean
algorithm must then also update the cofactors, while the half-gcd
computes the transformation matrices in any case.

Below the GCD cutoff, `nmod_poly_gcd` and `nmod_poly_resultant` use the
Euclidean algorithm, in Montgomery representation
(`nmod_poly_gcd_euclidean_redc_fast`) when $\ell \ge 32$ and
$\bits/2 \le \beta \le \bits - 2$; see
[](algorithms/arithmetic/montgomery). `nmod_poly_resultant` first tries
special methods for very short inputs.

The $2 \times 2$ polynomial matrix products of step 7 use eight
multiplications when an entry is shorter than 20, and Strassen's
seven-multiplication scheme otherwise. With `fft_small`, when all
entries have length at least 128, each of the eight entries is
transformed once and the four products are accumulated in the
transformed domain, so that 8 forward and 4 inverse transforms replace
the 24 transforms of eight separate FFT products; from length 1024 the
transforms are distributed over the thread pool.

The functions:

- `nmod_poly_hgcd` returns the half-gcd as in [](#def-hgcd), with
  $M = (m_{11}, m_{12}; m_{21}, m_{22})$ and the return value
  $\det M = \pm 1$; `_nmod_poly_hgcd` is the underlying function.
  `nmod_poly_hgcd_ref` is a simple quadratic reference implementation.
- `nmod_poly_gcd_hgcd`, `nmod_poly_xgcd_hgcd` and
  `nmod_poly_resultant_hgcd` force the half-gcd method;
  `nmod_poly_gcd`, `nmod_poly_xgcd`, `nmod_poly_resultant`,
  `nmod_poly_gcdinv`, `nmod_poly_invmod` and `nmod_poly_discriminant`
  use it above the cutoffs.

The functions are described, with examples, in [](nmod_poly/gcd) and
[](nmod_poly/resultant).

::: functions nmod_poly_hgcd nmod_poly_gcd_hgcd nmod_poly_xgcd_hgcd nmod_poly_resultant_hgcd :::
