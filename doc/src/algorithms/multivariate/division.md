---
title: Division, exact division and square roots
summary: Multivariate division with remainder with respect to a monomial ordering, and its heap implementation by Monagan and Pearce, in which a heap of size at most the number of terms of the divisor supplies the terms of the current remainder; early termination of exact division, quasi-division over the integers, division by several polynomials, parallel exact division, and the analogous heap algorithm for square roots with its non-square certificates.
complexity: Heap division $O(\ell_B \ell_Q \log \ell_B)$ monomial comparisons and $\ell_B \ell_Q$ coefficient multiplications for a quotient with $\ell_Q$ terms, $O(\ell_B)$ extra memory; heap square root $O(\ell_Q^2 \log \ell_Q)$.
order: 30
---

Division in $R[x_0, \ldots, x_{n-1}]$ is defined with respect to the
monomial ordering of the context ([](algorithms/multivariate/representation#orderings)).
Its main uses in FLINT are *exact division*, to compute $A/B$ when $B$
divides $A$ or to decide whether it does, which is needed throughout the
GCD and factoring code; and *division with remainder*, which is the
reduction step of Gröbner basis computations [@Buc1965; @CLO2015].
The heap algorithms of Monagan and Pearce [@MonPea2007; @MonPea2011]
compute the quotient term by term, in the time of the multiplication
$QB$, without forming intermediate remainders.

## Division with remainder {#divrem}

Over a field $K$, division of $A$ by $B \ne 0$ produces $Q$ and $R$ with
$A = QB + R$ and no term of $R$ divisible by $\operatorname{lm}(B)$:

::: algorithm Division with remainder {#alg-divrem}
Input: $A, B \in K[x]$, $B \ne 0$.
Output: $Q, R$ with $A = QB + R$, no monomial of $R$ divisible by
$\operatorname{lm}(B)$.

1. Set $Q = R = 0$ and $P = A$.
2. While $P \ne 0$: let $c x^\mu = \operatorname{lt}(P)$. If
   $\operatorname{lm}(B)$ divides $x^\mu$, set $t = (c/\operatorname{lc}(B))
   \, x^\mu / \operatorname{lm}(B)$, $Q \gets Q + t$ and $P \gets P - tB$;
   otherwise $R \gets R + c x^\mu$ and $P \gets P - c x^\mu$.
:::

Each step removes the leading term of $P$ and adds only smaller ones, so
the algorithm terminates because the ordering is a well-order, and $A =
QB + R + P$ is invariant. Quotient and remainder are unique
[@CLO2015, §2.3]. With several divisors $B_0, \ldots, B_{k-1}$
(`*_divrem_ideal`) the first $B_i$ whose leading monomial divides
$x^\mu$ is used; the remainder then depends on the order of the $B_i$
unless they form a Gröbner basis.

Over $\Z$ the quotient $c/\operatorname{lc}(B)$ is generally not an
integer. FLINT offers two variants, described with examples on
[](fmpz_mpoly/division):

- `fmpz_mpoly_divrem` uses the floor division $c = q \operatorname{lc}(B)
  + r$, takes $q$ as quotient coefficient and moves $r x^\mu$ to the
  remainder;
- `fmpz_mpoly_quasidivrem` keeps a scale $s$ and, when $\operatorname{lc}(B)
  \nmid c$, multiplies $P$, $Q$, $R$ and $s$ by $u = |\operatorname{lc}(B)|
  / \gcd(c, \operatorname{lc}(B))$ before continuing.

::: lemma Quasi-division {#lem-quasi}
`fmpz_mpoly_quasidivrem` returns $s > 0$, $Q$ and $R$ with $sA = QB + R$
and no monomial of $R$ divisible by $\operatorname{lm}(B)$; $s$ divides
$|\operatorname{lc}(B)|^{m}$, where $m$ is the number of terms of $Q$. Over
$\Q$, $Q/s$ and $R/s$ are the quotient and remainder of
[](#alg-divrem).
:::

::: proof
The invariant $sA = QB + R + P$ holds initially and is preserved by the
scaling (all four quantities are multiplied by $u$) and by the steps of
[](#alg-divrem), which after scaling divide exactly. A scaling happens at
most once per quotient term and multiplies $s$ by a divisor of
$|\operatorname{lc}(B)|$. Dividing $sA = QB + R$ by $s$ gives a
decomposition with the property of the remainder, which is unique.
:::

## Heap division {#heap}

In [](#alg-divrem), $P = A - Q_k B$ after $k$ quotient terms $q_0, \ldots,
q_{k-1}$. Its leading term is all that is needed, so instead of storing
$P$, the algorithm merges the terms of $A$ with the products $-q_j b_i$,
$i \ge 1$, in a heap ordered by monomial, exactly like the heap
multiplication of [](algorithms/multivariate/multiplication#alg-heap)
applied to $(B - \operatorname{lt}(B)) \times Q$ while $Q$ is still being
computed. (The products $q_j b_0$ cancel the leading terms by
construction and are never formed.)

::: algorithm Heap division {#alg-heap-div}
Input: $A$, $B$ with $\ell_A$, $\ell_B$ terms.
Output: $Q$ and $R$ as in [](#alg-divrem).

1. Put the first term of $A$ into the heap.
2. While the heap is not empty: pop all entries with the largest
   monomial $x^\mu$ and add their coefficients ($a$ for a term of $A$,
   $-b_i q_j$ for a product) into an accumulator $c$; reinsert successors
   (the next term of $A$; the products $(i + 1, j)$ and $(i, j + 1)$ by the
   rule of the multiplication heap).
3. If $c \ne 0$: if $\operatorname{lm}(B) \mid x^\mu$, the next quotient
   term is $q_k = (c / \operatorname{lc}(B)) \, x^{\mu} /
   \operatorname{lm}(B)$, and the product $(1, k)$ becomes available;
   otherwise $c x^\mu$ is a term of $R$.
:::

The product $(1, k)$ for a new quotient term must be inserted only when
the row of $q_k$ is allowed to start; Monagan and Pearce keep a counter
$s$ of rows waiting for the newest quotient term and insert it lazily,
which keeps the invariant of
[](algorithms/multivariate/multiplication#lem-heap) and the bound
on the heap size.

::: theorem Heap division {#thm-heap-div}
[](#alg-heap-div) computes $Q$ and $R$ using at most $\ell_A + (\ell_B -
1)\ell_Q$ heap insertions on a heap with at most $\ell_B$ entries, hence
$O((\ell_A + \ell_B \ell_Q) \log \ell_B)$ monomial comparisons, and
$O(\ell_B)$ memory besides $A$, $Q$ and $R$.
:::

::: proof
The terms entering the heap are the terms of $A$, one at a time, and the
products $b_i q_j$, $i \ge 1$; the latter are in rows indexed by $i$, each
with at most one entry in the heap as in the multiplication. The
monomial popped is the leading monomial of $P = A - Q_kB - R_k$, so the
steps coincide with those of [](#alg-divrem).
:::

This *divisor heap* is the variant of [@MonPea2011] suited to divisors
with few terms, the common case in exact division ($\ell_B \le \ell_Q$
for most GCD and factoring applications). Coefficients are accumulated in
three words when the sizes permit, as in multiplication.

### Exact division {#exact}

When $B$ is expected to divide $A$, the remainder is not wanted and the
division can stop at the first sign of inexactness.

::: lemma Early termination {#lem-early}
If $B$ divides $A$ with quotient $Q$, then in [](#alg-heap-div):

1. every leading monomial $x^\mu$ of $P$ is divisible by $\operatorname{lm}(B)$,
   and every accumulated coefficient is divisible by $\operatorname{lc}(B)$;
2. every such $x^\mu$ is at least the last monomial of $A$;
3. the monomials of $Q$ fit in the exponent fields of $A$;
4. $\deg_{x_i} B \le \deg_{x_i} A$ for each $i$.
:::

::: proof
If $B \mid A$, the remainder is zero, so no term goes to $R$, which is (1).
If $x^\mu = \operatorname{lm}(B) x^\nu$ with $x^\nu$ a monomial of $Q$, then
$x^\mu \ge \operatorname{tm}(B) x^\nu \ge \operatorname{tm}(B)
\operatorname{tm}(Q) = \operatorname{tm}(A)$, where tm is the last
(trailing) monomial, which is (2). Every monomial of $Q$ divides a
monomial of $A = QB$ in the sense that $\deg_{x_i} Q \le \deg_{x_i} A$
(degrees add in an integral domain), which gives (3) and (4).
:::

`fmpz_mpoly_divides` tests (4) before starting, computes in the packing
width of $A$ and returns 0 as soon as one of (1)–(3) fails, using
[the packed divisibility test](algorithms/multivariate/representation#lem-divides)
for (1) and the overflow test for (3). Since in typical uses most trial
divisions that fail do so on the first few terms, this is much faster
than computing a remainder. `fmpz_mpoly_divexact` assumes divisibility.

### Parallel exact division {#parallel}

Monagan and Pearce [@MonPea2010] parallelise exact division although the
quotient terms are produced sequentially: worker threads compute the
products $b_i q_j$ for ranges of the terms of $A$ (cut at monomials chosen
with `mpoly_divides_select_exps`) as soon as the quotient terms they need
exist, and subtract them from their range; the thread owning the current
leading range divides and publishes new quotient terms. FLINT uses this
in `fmpz_mpoly_divides_heap_threaded` (and the `nmod_mpoly` analogue)
when $A$ has more than 1024 terms, threads are available, and the
platform has a strongly ordered memory model (`FLINT_KNOW_STRONG_ORDER`),
which the lock-free publication of quotient terms requires.

### Dense array division {#array}

`fmpz_mpoly_divides_array` and `fmpz_mpoly_divrem_array` run
[](#alg-divrem) on a dense array of coefficients indexed by the
exponents, in lex or degree slices as for
[array multiplication](algorithms/multivariate/multiplication#array).
They are not selected automatically by `fmpz_mpoly_divides` and
`fmpz_mpoly_divrem`.

## Square roots {#sqrt}

If $A = Q^2$ with $Q = \sum_{j} q_j x^{\nu_j}$ sorted, then
$\operatorname{lt}(A) = q_0^2 x^{2\nu_0}$, and for $k \ge 1$ the
coefficient of $x^{\nu_0 + \nu_k}$ in $A - (q_0 x^{\nu_0} + \cdots + q_{k-1}
x^{\nu_{k-1}})^2$ is $2 q_0 q_k$, its leading term. So the terms of $Q$
can be computed one by one like the quotient terms of a division by $2
q_0 x^{\nu_0}$, with a heap supplying the products $2 q_i q_j$, $0 < i <
j$, of the partial square root. FLINT's `fmpz_mpoly_sqrt_heap`
and `nmod_mpoly_sqrt_heap` do this, accumulating the cross products once
and doubling.

A polynomial that is not a square is usually recognised quickly:

::: lemma Certificates of non-squareness {#lem-nonsquare}
Let $A \in \Z[x]$ be nonzero. $A$ is not a square in $\Z[x]$ if one of the
following holds:

1. some $\deg_{x_i} A$, or for a degree ordering the total degree of $A$,
   is odd;
2. $\operatorname{lc}(A)$ is not a square, or $\operatorname{lm}(A)$ is
   not a square monomial;
3. in the computation above, a monomial $x^\mu$ of the remainder is not
   divisible by $x^{\nu_0}$, the coefficient is not divisible by $2q_0$, or
   $x^\mu$ is smaller than $x^{\nu_0} \cdot \sqrt{\operatorname{tm}(A)}$;
4. there are a prime $p$ and $\alpha \in \F_p^n$ with $A(\alpha) \bmod p$ a
   quadratic non-residue, that is, with Jacobi symbol $\bigl(\frac{A(\alpha)}{p}\bigr) = -1$.
:::

::: proof
Degrees, leading terms and trailing terms of $Q^2$ are twice, squares of,
and squares of those of $Q$, in an integral domain. (3) is (1) and (2) of
[](#lem-early) for the division by $2 q_0 x^{\nu_0}$. For (4), $A = Q^2$
implies $A(\alpha) = Q(\alpha)^2$ in $\F_p$.
:::

`mpoly_is_proved_not_square` implements (1); `fmpz_mpoly_sqrt_heap` checks
(2) and (3) as it goes and, when the partial square root becomes longer
than $A$, tries (4) at a few random points modulo word-size primes before
continuing. `fmpz_mpoly_is_square` returns as soon as the square root is
determined.

## Generic and specialised implementations {#use}

The generic versions in `gr_mpoly` work for any coefficient ring:
`gr_mpoly_divides` (with a direct method for monomial divisors, and
`gr_mpoly_divides_heap` or `gr_mpoly_divides_heap_threaded` in integral
domains), `gr_mpoly_divrem_heap` and `gr_mpoly_divrem_ideal` for division
with remainder (typically over fields),
`gr_mpoly_quasidivrem_heap` and `gr_mpoly_quasidivrem_ideal_heap` for
pseudo-division, and the *weak* variants `gr_mpoly_divrem_weak` and
`gr_mpoly_divrem_ideal_weak`, which divide coefficients with remainder
when an exact division is impossible, like `fmpz_mpoly_divrem` over
$\Z$. The specialised modules
implement the same algorithms with packed coefficient accumulation:

- `fmpz_mpoly`: `fmpz_mpoly_divides` ([](#lem-early); heap:
  `fmpz_mpoly_divides_monagan_pearce`, threaded:
  `fmpz_mpoly_divides_heap_threaded`), `fmpz_mpoly_divrem` and
  `fmpz_mpoly_div` (floor division of coefficients, heap:
  `fmpz_mpoly_div_monagan_pearce`), `fmpz_mpoly_quasidivrem` and
  `fmpz_mpoly_quasidiv` ([](#lem-quasi), `fmpz_mpoly_quasidivrem_heap`),
  `fmpz_mpoly_divrem_ideal` and `fmpz_mpoly_quasidivrem_ideal`; array
  versions `fmpz_mpoly_divides_array` and `fmpz_mpoly_divrem_array`;
  square roots `fmpz_mpoly_sqrt_heap`, `fmpz_mpoly_sqrt`,
  `fmpz_mpoly_is_square` ([](#lem-nonsquare)).
- `fmpq_mpoly`: division by the integer parts, with the quasi-division
  providing the remainder over $\Q$.
- `nmod_mpoly`: `nmod_mpoly_divides` (heap, threaded heap, and a dense
  method `nmod_mpoly_divides_dense`), `nmod_mpoly_divrem` and
  `nmod_mpoly_divrem_ideal` (Monagan–Pearce heap), `nmod_mpoly_sqrt_heap`.
- `fmpz_mod_mpoly`, `fq_nmod_mpoly`, `fq_zech_mpoly`: heap division and
  square roots with the same structure.

::: functions fmpz_mpoly_divides fmpz_mpoly_divrem fmpz_mpoly_quasidivrem fmpz_mpoly_divrem_ideal fmpz_mpoly_divides_heap_threaded fmpz_mpoly_sqrt_heap mpoly_is_proved_not_square gr_mpoly_divides :::
