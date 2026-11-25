---
title: Rational reconstruction
summary: Recovers a fraction $n/d$ with $|n| \le N$ and $0 < d \le D$ from its residue modulo $m > 2ND$ by stopping the extended Euclidean algorithm at the first remainder below $N$. Used to return rational results from multimodular and p-adic algorithms.
complexity: $O(\log^2 m)$ bit operations with the Euclidean algorithm, $O(M(\log m) \log\log m)$ with the half-GCD algorithm.
order: 70
---

Multimodular and $p$-adic algorithms compute in $\Z/m\Z$, where a
rational number $n/d$ with $\gcd(d, m) = 1$ is represented by the residue
$a \equiv n d^{-1} \pmod m$. When the answer of a computation over $\Q$ —
the solution of a linear system, the coefficients of a polynomial — has
been found modulo a large enough $m$, *rational reconstruction* recovers
it from $a$. The method goes back to Wang [@WGD1982]; see also
[@CollinsEncarnacion1995], [@GG2003, §5.10] and [@Sho2009, §4.6].

Throughout, $m > 2$, $0 \le a < m$, and $N, D \ge 1$ are bounds with $2ND
< m$.

## Uniqueness {#uniqueness}

::: lemma {#lem-unique}
There is at most one fraction $n/d$ with $|n| \le N$, $0 < d \le D$ and
$n \equiv ad \pmod m$, up to the representation of the fraction.
:::

::: proof
If $n_1/d_1$ and $n_2/d_2$ both satisfy the conditions, then $n_1 d_2
\equiv a d_1 d_2 \equiv n_2 d_1 \pmod m$, and $|n_1 d_2 - n_2 d_1| \le 2ND
< m$, so $n_1 d_2 = n_2 d_1$.
:::

## Reconstructible residues {#density}

A reduced fraction $n/d$ has a residue modulo $m$ only when
$\gcd(d,m)=1$ [@Mon2004]. This condition also follows from
$n \equiv ad \pmod m$ and $\gcd(n,d)=1$: any prime dividing both $d$ and
$m$ would divide $n$.

Let $F_m(N,D)$ count the reduced fractions with $|n| \le N$,
$1 \le d \le D$ and $\gcd(d,m)=1$, for integer bounds $N,D \ge 1$.
Under $2ND < m$, [](#lem-unique) shows that these fractions have distinct
residues, so $F_m(N,D)$ is exactly the number of reconstructible residues.
Writing $\mu$ for the Möbius function, we have

$$
F_m(N,D) = 1 + 2
\sum_{\substack{1 \le d \le D \\ \gcd(d,m)=1}}
\sum_{e \mid d} \mu(e) \left\lfloor \frac{N}{e} \right\rfloor.
$$

The term 1 counts zero. For positive numerators, the identity
$\sum_{e \mid \gcd(n,d)} \mu(e) = 1$ if $\gcd(n,d)=1$, and 0 otherwise,
gives the inner sum; negative numerators give the factor 2.

If $m$ is prime and $D<m$, every denominator in the sum is a unit.
Interchanging the sums, with $E=\min(N,D)$, gives

$$
\begin{aligned}
F_m(N,D)
&= 1 + 2\sum_{e \le E} \mu(e)
  \left\lfloor \frac{N}{e} \right\rfloor
  \left\lfloor \frac{D}{e} \right\rfloor \\
&= \frac{12}{\pi^2} ND
   + O\bigl((N+D)\log(2E)\bigr).
\end{aligned}
$$

For $m=p^k$ with a fixed prime $p$ and balanced bounds $N=D=B$, the
denominator condition is $p \nmid d$. Writing $d=ej$ in the exact count
therefore gives

$$
\begin{aligned}
F_{p^k}(B,B)
&= 1 + 2\sum_{\substack{e \le B \\ p \nmid e}} \mu(e)
  \left\lfloor \frac{B}{e} \right\rfloor
  \left(\left\lfloor \frac{B}{e} \right\rfloor
       -\left\lfloor \frac{B}{pe} \right\rfloor\right) \\
&= \frac{12}{\pi^2}\frac{p}{p+1} B^2 + O_p(B\log B).
\end{aligned}
$$

To obtain the estimates, replace each floor by its argument plus $O(1)$.
The rounding errors sum to $O((N+D)\log(2E))$ and $O_p(B\log B)$,
respectively, since $\sum_{e \le E}1/e \le 1+\log E$.
Extending the main sums to infinity adds $O(ND/E)$ and $O_p(B)$, using
$\sum_{e>E}1/e^2 \le 1/E$. The constants follow from the Euler products

$$
\sum_{e \ge 1}\frac{\mu(e)}{e^2} = \frac{6}{\pi^2},
\qquad
\sum_{\substack{e \ge 1 \\ p \nmid e}}\frac{\mu(e)}{e^2}
= \frac{6}{\pi^2(1-p^{-2})}.
$$

Since $ND/E=\max(N,D)$, the tails are absorbed by the stated errors.
In the prime-power case, multiplying the second Euler product by
$2(1-1/p)B^2$ gives the factor $p/(p+1)$ in the leading term.

The prime-modulus count also explains the scale of the uniqueness bound:
if the bounds are enlarged so that $N,D \to \infty$ and $ND/m \to \infty$,
while $D<m$, there are more bounded fractions than residues.

For the balanced choice
$B=\lfloor\sqrt{\lfloor(m-1)/2\rfloor}\rfloor$, we have $B^2/m \to 1/2$.
The proportion of reconstructible residues therefore tends to
$6/\pi^2 \approx 61\%$ as $m$ tends to infinity through primes.
For $m=p^k$ with $p$ fixed and $k \to \infty$, it instead tends to
$(6/\pi^2)p/(p+1)$, giving $4/\pi^2 \approx 40.5\%$ for powers of 2.

For example, $m=64$ has balanced bound $B=5$. The admissible reduced
fractions have the following denominators and numerators:

| Denominator | Numerators | Count |
| --- | --- | --- |
| $1$ | $-5,-4,\ldots,5$ | 11 |
| $3$ | $\pm1,\pm2,\pm4,\pm5$ | 8 |
| $5$ | $\pm1,\pm2,\pm3,\pm4$ | 8 |

Since $2B^2=50<64$, these 27 fractions give distinct residues by
[](#lem-unique). Thus exactly 27 of the 64 residues are reconstructible.

## The extended Euclidean algorithm {#euclid}

Apply the extended Euclidean algorithm to $r_0 = m$ and $r_1 = a$: with
the quotients $q_i = \lfloor r_{i-1} / r_i \rfloor$,

$$
r_{i+1} = r_{i-1} - q_i r_i, \qquad s_{i+1} = s_{i-1} - q_i s_i, \qquad
t_{i+1} = t_{i-1} - q_i t_i,
$$

starting from $(s_0, t_0) = (1, 0)$ and $(s_1, t_1) = (0, 1)$, until
$r_{\ell + 1} = 0$. Then $r_i = s_i m + t_i a$, so $r_i \equiv t_i a
\pmod m$ for every $i$. We use the following standard facts
[@GG2003, §3.2, §4.6].

1. The $t_i$ alternate in sign and $|t_1| < |t_2| < \cdots < |t_{\ell+1}|$
   (with $|t_2| = q_1 \ge 1 = |t_1|$ the only possible equality).
2. $s_i t_{i+1} - s_{i+1} t_i = (-1)^i$, so $\gcd(s_i, t_i) = 1$.
3. $r_i |t_{i+1}| + r_{i+1} |t_i| = m$ for $0 \le i \le \ell$.
4. The convergents of the continued fraction of $a/m$ are the fractions
   $-s_i/t_i$, $1 \le i \le \ell + 1$.

Fact 3 follows from $r_i t_{i+1} - r_{i+1} t_i = (s_i t_{i+1} - s_{i+1}
t_i) m = \pm m$ and the alternating signs of the $t_i$.

::: algorithm Rational reconstruction {#alg-rr}
**Input:** $m > 2$, $0 \le a < m$, $N, D \ge 1$ with $2ND < m$.

**Output:** the fraction $n/d$ in lowest terms with $|n| \le N$, $0 < d
\le D$ and $n \equiv ad \pmod m$, or failure if none exists.

1. If $a \le N$, return $a/1$. If $m - a \le N$, return $(a - m)/1$.
2. Run the extended Euclidean algorithm on $(m, a)$, computing only the
   $r_i$ and $t_i$, up to the first index $j$ with $r_j \le N$.
3. If $|t_j| \le D$ and $\gcd(r_j, t_j) = 1$, return
   $\operatorname{sign}(t_j)\, r_j / |t_j|$. Otherwise fail.
:::

::: theorem {#thm-rr}
[](#alg-rr) is correct.
:::

::: proof
A returned fraction satisfies the conditions: in step 1 trivially, and in
step 3 because $r_j \equiv t_j a \pmod m$, $r_j \le N$, $|t_j| \le D$ and
$\gcd(r_j, t_j) = 1$. By [](#lem-unique) it is the only one.

Conversely, suppose that $n/d$ in lowest terms satisfies the conditions
and that step 1 does not apply, so $a > N$ and $j \ge 2$. Write $n = km +
da$. Then $\gcd(k, d)$ divides $n$ and $d$, so $\gcd(k, d) = 1$, and

$$
\left| \frac{a}{m} - \frac{-k}{d} \right| = \frac{|n|}{dm} \le
\frac{N}{dm} < \frac{1}{2d D} \le \frac{1}{2d^2}.
$$

By Legendre's theorem [@HarWri2008, Theorem 184], $-k/d$ is a convergent
of $a/m$, so by fact 4 and fact 2, $(k, d) = \pm (s_i, t_i)$ for some $i
\ge 1$, and then $n = km + da = \pm r_i$. Since $r_i = |n| \le N$, we have
$i \ge j$. If $i > j$, then by facts 3 and 1,

$$
m = r_j |t_{j+1}| + r_{j+1} |t_j| \le 2 r_j |t_{j+1}| \le 2N |t_i| = 2Nd
\le 2ND < m,
$$

a contradiction. Hence $i = j$, $(r_j, t_j) = \pm(n, d)$, and step 3
returns $n/d$.
:::

The remainders $r_i$ shrink and the cofactors $t_i$ grow as the
algorithm proceeds, with $r_{i-1} |t_i| \le m$; the bound $N$ decides
where the algorithm stops, trading the size of the numerator against the
size of the denominator. Monagan's *maximal quotient* variant
[@Mon2004] does not fix $N$ in advance but stops after the largest
quotient $q_i$, which succeeds with a modulus only slightly larger than
$2|n|d$ with high probability; FLINT implements the bounded version,
which is a proof when the bounds are known.

## Complexity {#complexity}

The Euclidean algorithm on $b$-bit numbers costs $O(b^2)$ bit operations.
Lehmer's method [@Leh1938] simulates several steps with word-size
leading parts of $r_{i-1}$ and $r_i$, so that most of the work is done
with single-word $2 \times 2$ matrices. The half-GCD algorithm
[@Sch1971], [@ThullYap1990] computes the quotient sequence up to a given
size recursively from the leading halves and costs $O(M(b) \log b)$. The
stopping rule $r_j \le N$ is a *partial* GCD, which the half-GCD algorithm
handles naturally: it reduces $(m, a)$ to a pair straddling a given size,
then continues with Lehmer steps.

## Several entries with a common denominator {#vectors}

When a vector or matrix of rationals is reconstructed, the entries often
share most of their denominator. Multiplying each residue by the
denominator $\delta$ of the entries found so far, the next entry is
usually an integer $\delta x \bmod m$ in the symmetric range $|\cdot| \le
N$, recognised without any GCD; only otherwise is the entry
reconstructed, with the denominator bound $D/\delta$, and $\delta$
multiplied by the new denominator.

## Use in FLINT {#use}

`_fmpq_reconstruct_fmpz_2` implements [](#alg-rr). After the quick
integer tests of step 1, it chooses by the size of $m$:

- for $m$ of one or two limbs, a Euclidean algorithm in one or two words,
  with a binary GCD for the final coprimality test;
- up to 12 limbs, Lehmer's method on stack arrays of limbs;
- for larger $m$, Lehmer's method on `fmpz`, switching to the half-GCD
  algorithm `_fmpq_hgcd` of the continued fraction code
  ([](algorithms/number-theory/continued-fractions)) when $m$ has at least
  `FMPQ_RECONSTRUCT_HGCD_CUTOFF` (500) more limbs than $N$, and to plain
  Euclidean steps when they differ by fewer than 3 limbs.

`_fmpq_reconstruct_fmpz_2_naive` is the plain algorithm on `fmpz`, used
to check the fast version in tests. `fmpq_reconstruct_fmpz` uses the
balanced bounds $N = D = \lfloor \sqrt{\lfloor (m-1)/2 \rfloor}
\rfloor$.

The reconstruction returns rational results in
`fmpq_mat_set_fmpz_mat_mod_fmpz`, `fmpz_mat_det_divisor` and the Dixon
solver `fmpz_mat_solve_dixon_den` (through
`_fmpz_mat_reconstruct_matwise`), which all use the common denominator
technique above, and in `fmpq_poly_interpolate_multi_mod` and
`fmpz_poly_roots_fmpq`. How the multimodular and $p$-adic linear algebra
algorithms choose the bounds and decide when to stop is described in
[multimodular and $p$-adic linear algebra](algorithms/linear-algebra/multimodular#reconstruction).

::: functions fmpq_reconstruct_fmpz fmpq_reconstruct_fmpz_2 fmpq_mod_fmpz :::
