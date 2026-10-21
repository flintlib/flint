---
title: Computing the partition function
summary: Computes p(n) exactly by the Hardy–Ramanujan–Rademacher series, truncated with Rademacher's error bound, with the exponential sums factored into products of cosines and each term evaluated in ball arithmetic at the precision it needs; small values come from Euler's pentagonal recurrence.
complexity: $O(n^{1/2} \log^{4+o(1)} n)$ bit operations for $p(n)$, quasi-optimal since $p(n)$ has $\Theta(n^{1/2})$ bits; $O(n^{3/2 + o(1)})$ for all of $p(0), \ldots, p(n)$.
order: 60
---

The partition function $p(n)$ counts the partitions of $n$, the ways of
writing $n$ as a sum of positive integers without regard to order. Its
generating function is Euler's product

$$
\sum_{n=0}^{\infty} p(n) x^n = \prod_{k=1}^{\infty} \frac{1}{1 - x^k}.
$$

The values grow subexponentially: $\log p(n) \sim \pi \sqrt{2n/3}$, so
$p(n)$ has about $3.7 \sqrt{n}$ bits. FLINT's methods follow
[@Joh2012]; for the theory see [@Andrews1976, Chapter 5] and
[@Apostol1997, Chapter 5].

## All values up to $n$ {#vector}

By Euler's pentagonal number theorem,

$$
\prod_{k=1}^{\infty} (1 - x^k) = \sum_{j=-\infty}^{\infty} (-1)^j
x^{j(3j-1)/2} = 1 - x - x^2 + x^5 + x^7 - x^{12} - x^{15} + \cdots,
$$

a series with $O(\sqrt n)$ nonzero terms below $x^n$. Inverting it gives
the recurrence

$$
p(n) = \sum_{j \ge 1} (-1)^{j+1} \bigl( p(n - j(3j-1)/2) + p(n -
j(3j+1)/2) \bigr),
$$

which computes $p(0), \ldots, p(n)$ with $O(n^{3/2})$ additions of numbers
of $O(\sqrt n)$ bits. FLINT instead inverts the pentagonal series with
fast power series inversion (`arith_number_of_partitions_vec`, and
modulo a word-size integer `arith_number_of_partitions_nmod_vec`), which
costs $O(M(n))$ operations on coefficients of $O(\sqrt n)$ bits. The
recurrence is used for single values up to $n = 416$, in single words:
since it involves only additions and subtractions, it is exact modulo
$2^{64}$, hence exact while $p(n) < 2^{64}$.

## The Hardy–Ramanujan–Rademacher formula {#hrr}

Hardy and Ramanujan [@HardyRamanujan1918] found an asymptotic series for
$p(n)$, and Rademacher [@Rademacher1937] modified it into a convergent
series. In the form used by FLINT, with

$$
C = \frac{\pi}{6} \sqrt{24n - 1}, \qquad U(x) = \cosh x - \frac{\sinh
x}{x},
$$

it reads

$$
p(n) = \sum_{k=1}^{N} A_k(n) \sqrt{\frac{3}{k}} \frac{4}{24n - 1}\,
U\!\left(\frac{C}{k}\right) + R(n, N),
$$

where $A_k(n)$ is the exponential sum

$$
A_k(n) = \sum_{\substack{0 \le h < k \\ \gcd(h, k) = 1}} e^{\pi i
[s(h,k) - 2nh/k]},
$$

$s(h, k)$ is a Dedekind sum (see `fmpq_dedekind_sum`), and the remainder
satisfies [@Rademacher1937]

$$
|R(n, N)| < M(n, N) = \frac{44 \pi^2}{225 \sqrt 3} N^{-1/2} +
\frac{\pi \sqrt 2}{75} \left( \frac{N}{n-1} \right)^{1/2} \sinh\left(
\frac{\pi}{N} \sqrt{\frac{2n}{3}} \right).
$$

The $k$-th term is of the order of $e^{C/k}$, so the first term gives
almost all the bits of $p(n)$ and the terms decrease rapidly.
$|A_k(n)| \le k$, and the bound $M(n, N)$ falls below $1/2$ for $N$ of
the order of $\sqrt n$: for $n = 10^9$, $M(n, 10^4) \approx 0.32$. Since
$p(n)$ is an integer, an approximation with total error below $1/2$
determines it.

## Factoring the exponential sums {#expsum}

Evaluating $A_k(n)$ from its definition costs $O(k)$ Dedekind sums,
$O(k^2)$ operations for each $k$ and $O(n^{3/2})$ in total, which would
dominate. Lehmer and Whiteman [@Whiteman1956] showed that $A_k(n)$ is
multiplicative in a twisted sense and that its prime power factors are
single cosines. FLINT uses the following form.

::: proposition {#prop-factor}
If $k = k_1 k_2$ with $\gcd(k_1, k_2) = 1$, then $A_k(n) = A_{k_1}(n_1)
A_{k_2}(n_2)$ for explicit $n_1, n_2$ determined by congruences modulo
$k_1$ and $k_2$. For a prime $p > 3$ and $k = p^\lambda$, let $r = 1 -
24n \bmod k$. Then

- if $r$ is not a square modulo $p$, $A_k(n) = 0$;
- if $p \mid r$, $A_k(n) = (3/p) \sqrt p$ for $\lambda = 1$ and $A_k(n)
  = 0$ for $\lambda > 1$;
- otherwise, with $m$ such that $(24m)^2 \equiv r \pmod k$,
  $$
  A_k(n) = 2 \left(\frac{3}{k}\right) \sqrt k \cos\left( \frac{4\pi
  m}{k} \right).
  $$

Similar formulas, involving square roots modulo $8k$ and $3k$, hold for
powers of 2 and 3.
:::

So after factoring $k$, $A_k(n)$ is a product of $\omega(k)$ cosines of
rational multiples of $\pi$, an algebraic prefactor $\pm 2^a
\sqrt{u/v}$, and a sign, computed with $O(\log k)$ modular square roots
and inverses (`arith_hrr_expsum_factored`, which returns the
symbolic product as a `trig_prod_t`). The cost of $A_k(n)$ becomes
polylogarithmic in $k$, and many $A_k(n)$ vanish.

## Precision and error control {#precision}

Each term is evaluated at a precision tied to its magnitude: the term $k$
has absolute value at most $\sqrt{12k}\, e^{C/k}/(24n - 1)$, and it must
be known to an absolute error of about $2^{-g}/N$ with a few guard bits
$g$, so the working precision is

$$
\mathit{prec}_k \approx \frac{C/k}{\log 2} + \log_2 \frac{\sqrt{12 k}}{24n
- 1} + \log_2(8N) + O(\log \log n),
$$

where the last term accounts for the number of operations, which depends
on the number of prime factors of $k$. The precision thus decreases like
$1/k$: the first term needs the full $3.7 \sqrt n$ bits, the term $k$
about $1/k$ of that, and most of the $N = O(\sqrt n)$ terms need only a
few dozen bits. The cost of all terms is dominated by
the first few, and Johansson [@Joh2012] shows a total cost of $O(n^{1/2}
\log^{4 + o(1)} n)$ for the evaluation with $N = O(\sqrt n)$ terms,
within polylogarithmic factors of the size of the output.

All operations are done in ball arithmetic, so the computed sum $S$
comes with a radius $\rho$ that bounds every rounding error. Adding
$M(n, N)$ to the radius gives a ball $[S \pm (\rho + M(n,N))]$ that is
guaranteed to contain $p(n)$; if it contains exactly one integer, that
integer is $p(n)$. This makes the result a proof; the implementation of
[@Joh2012] used MPFR with error estimates that were not fully
rigorous.

The angles of the cosines are reduced to $[0, \pi/4]$ (using sines near
$\pi/2$), so that every factor is known with full relative precision, and
$\exp(C/k)$ is computed from $\exp(C)$ by a $k$-th root when that is
cheaper than a new exponential.

## Use in FLINT {#use}

- `partitions_fmpz_ui` and `partitions_fmpz_fmpz` use the table
  `partitions_lookup` for $n < 128$, the pentagonal recurrence in words
  for $n < 417$, and `mp_real_partitions_hrr` above.
- `mp_real_partitions_hrr` chooses the number of terms with
  $M(n, N) < 0.4$ (by bisection on a cheap estimate, then exactly),
  factors each $A_k(n)$ with `arith_hrr_expsum_factored`, and sorts the
  terms by the precision they need: terms of up to about 200 bits are
  evaluated in batches as `dfloat` balls of one to four doubles with the
  vectorised `gr_vec` functions, and the others as `mp_real` balls of the
  needed number of limbs. Above $n = 10^8$, the constants $\pi$, $C$ and
  $e^C$ are computed with all threads, and the terms are distributed over
  the threads of `flint_set_num_threads` in ranges of equal estimated
  cost. The remainder bound comes from `partitions_rademacher_bound`, and
  the result from `mp_real_unique_integer`.
- `arith_number_of_partitions` is a wrapper of `partitions_fmpz_ui`;
  `arith_number_of_partitions_vec` and
  `arith_number_of_partitions_nmod_vec` invert the pentagonal series.
- `partitions_leading_fmpz` evaluates the dominant part of the first
  term, and `arb_partitions_fmpz` approximations of $p(n)$ to a given
  precision.

The generic functions `gr_partitions_ui`, `gr_partitions_fmpz` and
`gr_partitions_vec` call `partitions_fmpz_ui`, `partitions_fmpz_fmpz` and
the vector functions over $\Z$ and $\Z/n\Z$, and `arb_partitions_fmpz`
in real and complex rings.

::: functions partitions_fmpz_fmpz partitions_rademacher_bound arith_hrr_expsum_factored arith_number_of_partitions_vec :::
