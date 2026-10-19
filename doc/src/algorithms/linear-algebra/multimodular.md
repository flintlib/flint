---
title: Multimodular and p-adic linear algebra over the integers
summary: How FLINT computes determinants, solutions of linear systems, ranks, echelon forms, characteristic and minimal polynomials and products of integer matrices by working modulo word-size primes or prime powers, with Hadamard bounds, rational reconstruction, output-sensitive termination and exact certificates for every answer.
complexity: For an $n \times n$ matrix with $b$-bit entries and $\beta$-bit words, $O(n^{\omega + 1}(b + \log n)/\beta)$ word operations for the determinant and for multimodular solving, and $O(n^3 (b + \log n)^2/\beta^2 + n^\omega)$ for $p$-adic solving with one right-hand side; the output-sensitive solvers replace $n(b + \log n)$ by the actual size of the solution.
order: 50
---

Exact linear algebra over $\Z$ suffers from the growth of intermediate
results: Gaussian elimination over $\Q$ creates fractions whose sizes
grow quickly, and even fraction-free elimination
([](algorithms/linear-algebra/fraction-free)) works with entries as large
as the final result in each of its $O(n^3)$ operations. Modular methods
avoid this. They solve the problem modulo word-size primes, where
`nmod_mat` runs at the speed of floating-point arithmetic
([](algorithms/linear-algebra/matrix-multiplication)), and recover the
integer or rational result from its residues. FLINT does this for
determinants, linear systems, ranks and echelon forms, characteristic
polynomials and products, and makes every answer rigorous.

The methods are specific to $\Z$ and $\Q$. There is no generic
multimodular framework in `gr_mat`; instead, the generic functions
dispatch to these ones for the ring of integers: `gr_mat_mul`,
`gr_mat_det`, `gr_mat_charpoly` and `gr_mat_permanent` call the
corresponding `fmpz_mat` functions.

Throughout, $A$ is an $n \times n$ (or $m \times n$) integer matrix with
entries bounded by $2^b$ in absolute value, $\lVert \cdot \rVert$ is the
maximum norm of a matrix, $\lVert \cdot \rVert_2$ the Euclidean norm of a
row or column, and $\beta$ is the word size.

## Bounds {#bounds}

::: theorem Hadamard's inequality {#thm-hadamard}
For a square matrix $A$ with rows $a_1, \ldots, a_n$ and columns
$c_1, \ldots, c_n$,
$$
|\det A| \le \min\Big(\prod_i \lVert a_i \rVert_2, \; \prod_j \lVert c_j \rVert_2\Big).
$$
:::

::: proof
The Gram–Schmidt orthogonalisation $a_1^*, \ldots, a_n^*$ of the rows
satisfies $|\det A| = \prod_i \lVert a_i^* \rVert_2$ and
$\lVert a_i^* \rVert_2 \le \lVert a_i \rVert_2$ [@Had1893]. Apply this to
$A$ and to $A^T$.
:::

With entries below $2^b$, this gives $|\det A| \le n^{n/2} 2^{nb}$, a
number of about $n(b + \frac12 \log_2 n)$ bits. `fmpz_mat_det_bound`
evaluates the bound with upward-rounded `mag_t` arithmetic.
`fmpz_mat_det_bound_submatrix` omits zero rows and columns from the
products, which bounds every minor of a (possibly rectangular) matrix: the
rows of a $k \times k$ submatrix are no longer than the corresponding rows
of $A$, and the same holds for columns.

::: lemma Solution bound {#lem-cramer}
Let $A$ be nonsingular with columns $c_1, \ldots, c_n$ and let $B$ have
columns $b_1, \ldots, b_k$. Every entry of $\det(A) A^{-1} B$ is at most
$$
N = \prod_i \lVert c_i \rVert_2 \cdot \max_l \lVert b_l \rVert_2
$$
in absolute value.
:::

::: proof
By Cramer's rule entry $(j, l)$ is $\det A_j(b_l)$, where $A_j(b_l)$ is
$A$ with column $j$ replaced by $b_l$. By [](#thm-hadamard) applied to
columns it is at most $\lVert b_l \rVert_2 \prod_{i \ne j}
\lVert c_i \rVert_2 \le N$, since the nonzero integer column $c_j$ has
$\lVert c_j \rVert_2 \ge 1$.
:::

`fmpz_mat_solve_bound` returns $N$ and the bound $D$ of
[](#thm-hadamard) for $|\det A|$. The row product cannot replace the
column product in $N$: it bounds $\det A$ but not the cofactors along a
column.

## Chinese remaindering {#crt}

If $x \in \Z$ satisfies $|x| \le H$ and $p_1, \ldots, p_k$ are pairwise
coprime with $M = p_1 \cdots p_k > 2H$, then $x$ is the unique integer in
$(-M/2, M/2]$ congruent to its residues. For matrices,
`fmpz_mat_multi_mod_ui` reduces all entries modulo $k$ primes and
`fmpz_mat_multi_CRT_ui` reconstructs them, both through the subproduct
tree of `fmpz_comb_t`: reducing a $kc$-bit integer modulo $k$ primes of
$c$ bits by a tree of remainders costs $O(M(kc) \log k)$ word operations
instead of the $O(k^2)$ of reducing modulo each prime separately, and the
same holds for reconstruction [@GG2003, §10.3]. Both functions distribute
the entries over threads. FLINT uses primes of 59 to 63 bits, the size
for which `nmod_mat` has its fastest kernels
(`NMOD_MAT_OPTIMAL_MODULUS_BITS` = 59).

::: definition {#def-bad}
A prime $p$ is *bad* for a computation if the reduction modulo $p$ of a
quantity the algorithm relies on behaves differently from the integer
quantity: for example if $p$ divides $\det A \ne 0$, so that $A$ is
singular modulo $p$, or if the rank of $A$ modulo $p$ is less than its
rank over $\Q$.
:::

The rank of $A$ modulo $p$ is at most its rank $r$ over $\Q$, with
equality unless $p$ divides all $r \times r$ minors of $A$, in particular
some fixed nonzero minor $\Delta$. Since $|\Delta|$ is at most the bound
$H$ of `fmpz_mat_det_bound_submatrix`, the product of the distinct bad
primes is at most $H$, so a set of distinct primes whose product exceeds
$H$ contains a good prime. FLINT uses this both ways: random large primes
are good with overwhelming probability, and when all primes of such a set
show a rank drop, the rank drop is real.

## Determinants {#det}

::: algorithm Determinant with a known divisor {#alg-det}
**Input:** $A$, $n \times n$; a positive integer $d$ dividing $\det A$.

**Output:** $\det A$.

1. $H \gets \lceil 2 \cdot \text{Hadamard}(A) / d \rceil$, $x \gets 0$,
   $M \gets 1$.
2. While $M \le H$: take the next prime $p > 2^{59}$ not dividing $d$ and
   compute an LU decomposition of $A$ modulo $p$.
   1. If the rank modulo $p$ is less than $n$ and $\det A = 0$ can be
      certified by [](#prop-kernel), return 0.
   2. $x_p \gets \pm \prod_i U_{ii} \cdot d^{-1} \bmod p$, the sign being
      that of the row permutation.
   3. $x \gets$ the integer in $(-Mp/2, Mp/2]$ congruent to $x$ modulo $M$
      and to $x_p$ modulo $p$; $M \gets Mp$.
3. Return $dx$.
:::

This is `fmpz_mat_det_modular_given_divisor`; `fmpz_mat_det_modular`
uses $d = 1$. Since $|\det(A)/d| \le H/2 < M/2$ at the end, the result is
correct. The number of primes is about
$\log_2(H)/59 \approx n(b + \frac12 \log_2 n)/59$, and each costs a
reduction of $A$ and an LU decomposition modulo $p$, so the total is
$O(n^{\omega + 1}(b + \log n)/\beta)$ word operations, against
$O(n^5 b^2)$ with classical integer arithmetic for fraction-free
elimination.

### Early termination {#early}

With `proved` = 0, the loop stops as soon as the reconstructed value has
not changed while primes with more than 100 bits in total were added. If
$\det(A)/d$ is much smaller than the Hadamard bound, which is typical for
structured input, this saves most of the primes. A value that stays the
same while $k$ further primes are added is wrong only if all $k$ primes
divide the difference between the true and the reconstructed value. For
random primes this is very unlikely [@Kal2002], but it is not a proof,
and `fmpz_mat_det` does not use it.

### A large divisor {#divisor}

The number of primes depends on $\det(A)/d$, so a large known divisor
$d$ saves work.

::: lemma {#lem-sn}
Let $A$ be nonsingular with Smith normal form
$\operatorname{diag}(s_1, \ldots, s_n)$. Then $s_n A^{-1}$ is integral.
Hence for every integer vector $v$ the least common denominator of
$A^{-1} v$ divides $s_n$, which divides $\det A$.
:::

::: proof
With $A = U^{-1} S V^{-1}$ for unimodular $U$ and $V$,
$A^{-1} = V S^{-1} U$, and $s_n S^{-1}$ is integral since $s_i \mid s_n$
for all $i$.
:::

`fmpz_mat_det_divisor` solves $Ax = v$ for $v = (-1, 1, -1, \ldots)^T$ by
$p$-adic lifting ([](#dixon)) and returns the least common multiple $d$
of the denominators of $x$, a divisor of $s_n$. For random matrices,
$s_n$ is $\det(A)$ up to a small factor and $d = s_n$ with high
probability [@AbbottBronsteinMulders1999; @Pan1988; @EGV2000], so that
$\det(A)/d$ needs only one or two primes and the $p$-adic solve, which
costs about as much as $O(n)$ primes, dominates.
`fmpz_mat_det_modular_accelerated` uses this, and `fmpz_mat_det` calls it
for $n \ge 60$ when $n$ is at least the bit size of the entries. For
matrices whose determinant is spread over many invariant factors, $d$ is
small and nothing is gained. Asymptotically faster determinant algorithms
exist [@EGV2000; @KV2004; @Sto2005] but are not implemented.

## Rank and singularity {#rank}

Full rank modulo a prime proves full rank over $\Q$. A rank deficiency
modulo $p$ is either real or caused by a bad prime, and is certified as
follows.

::: proposition Kernel vector certificate {#prop-kernel}
Let $A$ be $m \times n$, let $p$ be a prime and let $r < n$ be the rank
of $A$ modulo $p$, with pivot rows $i_1, \ldots, i_r$ and pivot columns
$j_1, \ldots, j_r$, so that the submatrix $P = A[\{i_k\}, \{j_k\}]$ is
nonsingular modulo $p$. Let $f$ be a non-pivot column and let
$z \in \Q^r$ solve $Pz = -A[\{i_k\}, f]$. Define $v \in \Q^n$ by
$v_{j_k} = z_k$, $v_f = 1$ and $v_j = 0$ otherwise. If $Av = 0$, then $A$
has rank less than $n$ over $\Q$. If the rank of $A$ over $\Q$ is $r$,
then $Av = 0$.
:::

::: proof
$P$ is nonsingular modulo $p$, hence over $\Q$, so $z$ exists, and
$v \ne 0$ since $v_f = 1$; this proves the first claim. If $A$ has rank
$r$, its rows are combinations of the rows $i_1, \ldots, i_r$, which
annihilate $v$ by construction.
:::

`_fmpz_mat_certify_singular_lu_mod_p` computes $z$ exactly (by Cramer's
rule, fraction-free elimination or $p$-adic solving with the same prime)
and checks the $m - r$ other rows. `fmpz_mat_is_singular` (with primes of
17 bits), the determinant and the nonsingular solvers use it. If the
certificate fails, another prime is tried, and once the product of the
primes with a rank drop exceeds the Hadamard bound, $A$ is singular by
the argument of [](#crt).

For the rank itself and the reduced row echelon form, the analogous
certificate checks that all non-pivot rows lie in the span of the pivot
rows.

::: proposition Echelon form certificate {#prop-rref}
In the situation of [](#prop-kernel), let $E$ be the $r \times n$ matrix
$P^{-1} A[\{i_k\}, \cdot]$, which has an identity matrix in the pivot
columns, and let $N$ be the $n \times (n - r)$ matrix whose columns are
the nullspace vectors defined by $E$, one per non-pivot column. If $E$
is in reduced row echelon form and $FN = 0$, where $F$ consists of the
non-pivot rows of $A$, then $A$ has rank $r$ and its reduced row echelon
form is $E$ followed by zero rows.
:::

::: proof
The rows of $E$ span the same space $W$ as the pivot rows of $A$, and the
columns of $N$ span its orthogonal complement. $FN = 0$ puts the rows of
$F$ in the orthogonal complement of that complement, which is $W$. So
the row space of $A$ is $W$, of dimension $r$, and $E$ is its reduced
echelon basis.
:::

`fmpz_mat_rref_mul` and `fmpz_mat_rank_certify_lu_mod_p` implement this
with primes of 17 bits, following [@Stein2007]; $E$ is computed with
`fmpz_mat_solve`, so that $E$ is obtained as an integer matrix and a
denominator.

## Solving linear systems {#solve}

### p-adic lifting {#dixon}

::: algorithm Dixon's algorithm {#alg-dixon}
**Input:** nonsingular $A$, $n \times n$; $B$, $n \times k$; a prime $p$
not dividing $\det A$; a number of steps $s$.

**Output:** $X \equiv A^{-1} B \pmod{p^s}$.

1. Compute $C = A^{-1} \bmod p$, or an LU decomposition of $A$ modulo
   $p$; $R \gets B$; $X \gets 0$.
2. For $i = 0, \ldots, s - 1$:
   1. $Y \gets C R \bmod p$, with entries in $[0, p)$;
   2. $X \gets X + p^i Y$;
   3. $R \gets (R - AY)/p$, an exact division.
3. Return $X$.
:::

::: lemma {#lem-dixon}
After step $i$, $AX = B - p^{i+1} R$; in particular
$AX \equiv B \pmod{p^{i+1}}$. The entries of $R$ stay bounded by
$n \lVert A \rVert + \lVert B \rVert / p^{i+1}$.
:::

::: proof
By induction on $i$: $R - AY \equiv R - ACR \equiv 0 \pmod p$, so the
division is exact, and $A(X + p^i Y) = B - p^i R + p^i AY =
B - p^{i+1} (R - AY)/p$. For the bound,
$\lVert (R - AY)/p \rVert \le (\lVert R \rVert + n \lVert A \rVert
(p - 1))/p$.
:::

The rational solution $X^* = A^{-1} B$ has entries $u/\det(A)$ with
$|u| \le N$ by [](#lem-cramer), so after $s$ steps with $p^s > 2ND$ it is
determined by rational reconstruction ([](#reconstruction)). Each step
costs a product with $C$ (or two triangular solves) modulo $p$ and a
product of $A$ with an $n \times k$ matrix of word-size entries, about
$n^2 k \lceil b/\beta \rceil$ word operations. With
$s \approx n(b + \log n)/\log_2 p$ steps, the total for $k = 1$ is
$O(n^3 (b + \log n)^2/\beta^2)$ word operations plus one LU
decomposition modulo $p$ [@Dix1982; @MC1979]. For one right-hand side
this beats the $O(n^{\omega + 1} b/\beta)$ of the multimodular method.
Storjohann's high-order lifting reaches $O^\sim(n^\omega b)$ bit
operations [@Sto2003; @Sto2005] but is not implemented.

`fmpz_mat_solve_dixon` is the plain algorithm: it lifts to the worst-case
bound $2\max(N, D)^2$ and returns $X$ modulo $p^s$ without reconstruction.
`fmpz_mat_solve_dixon_den` refines it:

- It starts with a prime near $2^{30}$, for which the LU decomposition
  and the triangular solves are about twice as fast as with a 59-bit
  prime, at the same cost per lifted bit.
- The products $AY$ are computed with `fmpz_mat_mul`, with the classical
  algorithm (one multiplication of a multiword integer by a word per
  entry), or multimodularly from precomputed residues of $A$, depending on
  $n$, $k$ and $b$.
- When the solution turns out to be long, the digit is increased from $p$
  to $q = p^e$: the inverse $C = A^{-1} \bmod q$ is obtained by Newton
  iteration, $C \gets C(2I - AC)$ modulo increasing powers of $p$, and the
  products $CR$ and $AY$ are computed multimodularly with tables of
  residues of $C$ and $A$. The digit size is doubled while the lifting
  still to be done is estimated to cost at least twice as much as the
  change, which bounds the overhead by a factor 1.5 whatever the size of
  the solution.
- The digits are accumulated in a binary counter of partial sums, so that
  forming $X$ modulo $p^i$ costs $O(M(i \log p) \log i)$ per entry rather
  than $O(i^2)$.
- Termination is output-sensitive ([](#termination)).

### Chinese remaindering {#multi-mod}

`fmpz_mat_solve_multi_mod_den` solves $AX = B$ modulo primes
$p > 2^{59}$ with `nmod_mat_solve`, in batches processed in parallel,
and combines the residues of $X^* = A^{-1} B$ modulo the product $M$ of
the primes. The residues are combined only when a reconstruction is
attempted, so the total cost of the Chinese remaindering is that of a few
full reconstructions. Each prime costs $O(n^\omega + n^2 k)$ word
operations and a reduction of $A$ and $B$, so the total is
$O(n^{\omega + 1} (b + \log n)/\beta)$ for $k \le n$. When $k \ge n/2$
this beats $p$-adic lifting, whose steps cost $O(n^2 k b/\beta)$ each, and
`fmpz_mat_solve` uses it then.

### Rational reconstruction {#reconstruction}

::: theorem {#thm-ratrec}
Let $m \ge 2$, $N, D \ge 1$ with $2ND < m$, and $x \in \Z$. There is at
most one fraction $u/v$ with $|u| \le N$, $0 < v \le D$, $\gcd(v, m) = 1$
and $u \equiv vx \pmod m$. If it exists, it is found by the extended
Euclidean algorithm applied to $m$ and $x$, stopped at the first
remainder at most $N$.
:::

::: proof
If $u_1/v_1$ and $u_2/v_2$ both qualify, then $u_1 v_2 \equiv x v_1 v_2
\equiv u_2 v_1 \pmod m$ and $|u_1 v_2 - u_2 v_1| \le 2ND < m$, so
$u_1 v_2 = u_2 v_1$. For the statement about the Euclidean algorithm see
[@WGD1982] and [@GG2003, §5.10].
:::

Reconstructing every entry separately would need $m > 2ND$ for each. FLINT
reconstructs a whole matrix with a running common denominator
(`_fmpz_mat_reconstruct_matwise`): each entry $x$ is first multiplied by
the denominator $d$ found so far, and if $dx \bmod m$ in the symmetric
range is at most $N$, it is accepted as the numerator; otherwise the entry
is reconstructed with denominator bound $D/d$, and the new factor of the
denominator is multiplied into $d$ and into the entries done so far. Most
entries then cost one multiplication and one reduction [@CheSto2005].

### Output-sensitive termination {#termination}

The bounds $N$ and $D$ are often far larger than the actual solution. The
solvers therefore stop as soon as the solution is determined, using a
certificate that needs no a priori bound.

::: theorem Size certificate {#thm-certificate}
Let $Z$ be an integer $n \times k$ matrix and $d$ a nonzero integer with
$AZ \equiv dB \pmod m$. If
$$
n \lVert A \rVert \lVert Z \rVert + |d| \lVert B \rVert < m,
$$
then $AZ = dB$, so that $Z/d$ is the solution of $AX = B$.
:::

::: proof
Each entry of $AZ - dB$ is a sum of $n$ products bounded by
$\lVert A \rVert \lVert Z \rVert$ and one term bounded by
$|d| \lVert B \rVert$, so it is less than $m$ in absolute value. It is
divisible by $m$, hence zero.
:::

This criterion goes back to Cabay [@Cab1971] and is used in [@CheSto2005;
@Stef2010]. To apply it, the solvers must find a candidate $(Z, d)$ early,
without knowing the size of the solution:

1. *Probe.* A random combination $s = \sum_{i,j} u_{ij} x_{ij} \bmod m$ of
   the entries of the current solution modulo $m$, with weights
   $u_{ij} \in [1, 16]$, is maintained incrementally at negligible cost.
   It is the residue of the rational number
   $s^* = \sum u_{ij} X^*_{ij}$, whose denominator is the common
   denominator of the solution with high probability, and whose numerator
   has about the size of the largest numerator.
2. *Probe reconstruction* (`_fmpz_mat_solve_probe_reconstruct`). First $s$
   and $s d'$, for the denominator $d'$ found by the previous probe, are
   tested for being small in the symmetric range, with a margin of 30
   bits. Otherwise the continued fraction expansion of $s/m$ is computed.
   The true value $u/v$ is the convergent preceding a large partial
   quotient as soon as $m > 2^{30} |u| v$, whatever the balance between
   $u$ and $v$, which halves the precision that a fixed split
   $N = D = \sqrt{m/2}$ would need [@Mon2004]. A partial quotient of at
   least 30 bits is required, so a spurious result has probability about
   $2^{-30}$.
3. *Attempt.* When the probe has succeeded, and either the certificate can
   hold at the current modulus or an explicit check $AZ = dB$ is estimated
   to cost less than the lifting steps still missing, the whole solution
   is reconstructed with a denominator bound $2^{16}$ times the probe
   denominator (the margin doubles after failures) and checked with
   [](#thm-certificate) or by multiplication. During the first eight steps
   attempts are made without a probe, since they are cheap and the probe
   is not yet informative.
4. *Fallback.* When $m > 2ND$, reconstruction with the global bounds is
   unique and cannot fail.

A spurious probe costs time but never produces a wrong result, since
every accepted solution satisfies [](#thm-certificate) or has been
verified by a multiplication. Finally the gcd of the denominator and the
entries is divided out, so that the denominator is the least common
denominator.

## Consistency of general systems {#can-solve}

For an $m \times n$ matrix $A$ of rank $r$, the system $AX = B$ is solved
on the pivot rows and columns found modulo a prime, and the remaining
$m - r$ equations are checked exactly. If they hold, the solution is
correct. If one fails, either the system is inconsistent or the prime is
bad.

::: proposition Inconsistency certificate {#prop-inconsistent}
If there is a vector $y$ with $y^T A = 0$ and $y^T B \ne 0$, then
$AX = B$ has no solution.
:::

::: proof
$y^T B = y^T A X = 0$ for any solution $X$.
:::

Conversely, by the Fredholm alternative such a $y$ exists whenever the
system is inconsistent. `fmpz_mat_can_solve_multi_mod_den` constructs a
candidate from a failing equation $i_0$: $y$ has a denominator in
position $i_0$ and minus the solution of the transposed pivot system in
the pivot rows, so that $y^T A$ vanishes on the pivot columns by
construction, and $y^T A = 0$ is then checked exactly on the other columns
[@GieLobSau1998; @MulSto2004]. If the check fails, the prime was bad and
another one is tried; exact work is redone only when the rank increases.
Once the product of the primes exceeds the bound on all minors of $A$,
some prime was good ([](#crt)) and its verdict, "inconsistent", is final.

## Characteristic and minimal polynomials {#charpoly}

`fmpz_mat_charpoly_modular` uses a fixed number of primes determined by
the coefficient bound of `fmpz_mat_charpoly_bound`, computes the
characteristic polynomials modulo the primes in parallel
([](algorithms/linear-algebra/charpoly)), and reconstructs them. There are
no bad primes, since the characteristic polynomial commutes with
reduction.

The minimal polynomial does not: modulo a bad prime its degree drops.
`fmpz_mat_minpoly_modular` keeps the primes of largest degree and, when
the reconstruction stabilises, checks the candidate $f$ exactly:
$f(A) e_i = 0$ for the unit vectors $e_i$ used as Krylov generators
modulo the primes. Their Krylov spaces span $\mathbb{F}_p^n$, hence
$\Q^n$, so $f(A) = 0$ and $\mu_A \mid f$; and $\deg f \le \deg \mu_A$
because the minimal polynomial modulo $p$ divides the reduction of
$\mu_A$. Hence $f = \mu_A$. The number of primes is capped by a bound of
Dumas [@Dum2007] on the coefficients of $\mu_A$.

## Matrix multiplication {#mul}

The product $C = AB$ of integer matrices is determined by its residues
modulo numbers whose product exceeds $2\lVert C \rVert$, and there are no
bad primes. `fmpz_mat_mul_multi_mod` and `fmpz_mat_mul_blas` use this;
see [](algorithms/linear-algebra/integer-matrix-multiplication).

## Use in FLINT {#use}

| Task | Functions |
|---|---|
| Determinant | `fmpz_mat_det_modular`, `fmpz_mat_det_modular_accelerated`, `fmpz_mat_det_modular_given_divisor`, `fmpz_mat_det_divisor` |
| Singularity, rank | `fmpz_mat_is_singular`, `fmpz_mat_rank`, `fmpz_mat_rank_certify_lu_mod_p`, `_fmpz_mat_certify_singular_lu_mod_p` |
| Echelon form | `fmpz_mat_rref_mul`, `fmpz_mat_rref_upper_certify_lu_mod_p` |
| Nonsingular systems | `fmpz_mat_solve_dixon`, `fmpz_mat_solve_dixon_den`, `fmpz_mat_solve_multi_mod_den` |
| General systems | `fmpz_mat_can_solve_dixon_den`, `fmpz_mat_can_solve_multi_mod_den`, `fmpz_mat_can_solve` |
| Rational matrices | `fmpq_mat_solve_dixon`, `fmpq_mat_solve_multi_mod`, `fmpq_mat_can_solve` and variants |
| Polynomials | `fmpz_mat_charpoly_modular`, `fmpz_mat_minpoly_modular` |
| Products | `fmpz_mat_mul_multi_mod`, `fmpz_mat_mul_blas` |
| Building blocks | `fmpz_mat_multi_mod_ui`, `fmpz_mat_multi_CRT_ui`, `fmpz_mat_CRT_ui`, `fmpq_mat_set_fmpz_mat_mod_fmpz` |
