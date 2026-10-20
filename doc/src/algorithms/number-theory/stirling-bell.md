---
title: Stirling and Bell numbers
summary: Isolated Stirling numbers by binary splitting of rising factorials, power sums or power series, chosen by k/n; rows by convolution; Bell numbers by Dobiński's formula, a multimodular sum, the Bell triangle or generating functions modulo primes.
complexity: An isolated S(n,k) or B_n costs about $O(n^2 \log^{O(1)} n)$ bit operations; a row of Stirling numbers or a vector of Bell numbers modulo a word-size prime costs $O(M(n) \log n)$ word operations.
order: 80
---

The Stirling numbers $s(n,k)$ and $S(n,k)$ and the Bell numbers $B_n$
are defined on [](arith/stirling) and [](arith/bell). All of these numbers have
$\Theta(n \log n)$ bits for typical $k$, so the best possible cost of an
isolated value is quasi-linear in $n \log n$; the methods below are
simpler and, apart from the binary splitting method for $s(n,k)$, cost
about $n^2$ word operations, which is adequate up to $n$ of the order of
$10^5$. The references for the combinatorics are [@Comtet1974] and
[@GKP1994, §6.1].

## Stirling numbers of the first kind {#first-kind}

The unsigned numbers are the coefficients of the rising factorial:

$$
\sum_{k=0}^n |s(n,k)|\, x^k = x(x+1)\cdots(x+n-1).
$$

A row is computed by binary splitting: the product of the $n - 1$ linear
factors $x + j$, $1 \le j < n$, split recursively in halves, with fast
polynomial multiplication at the top levels. The coefficients of the
product of the factors with $j \in [a, b)$ have at most $(b-a)
\log_2 \max(b, b - a)$ bits, so the base cases with small products are
done in single words, and the total cost is $O(M(n^2 \log n) \log n)$
bit operations, quasi-linear in the size of the row.

For one coefficient $|s(n,k)|$, only a truncated product is needed: for
$k < n/2$ the coefficients up to $x^{k-1}$ of $\prod_{j=1}^{n-1} (x + j)$,
and otherwise the coefficients up to $x^{n-k}$ of the reversed product
$\prod_{j=1}^{n-1} (1 + jx)$. The last multiplication of the binary
splitting tree computes only the wanted coefficient, as a dot product.
For $k$ very close to $n$ the exponential generating function
$$
|s(n,k)| = \frac{n!}{k!} [x^{n-k}] \left(\frac{-\log(1-x)}{x}\right)^k
$$
is cheaper, since the power series has length only $n - k + 1$; FLINT
uses it for $n > 140$ and $k > 0.87 n$. The signed numbers are $s(n,k) =
(-1)^{n-k} |s(n,k)|$.

## Stirling numbers of the second kind {#second-kind}

By inclusion–exclusion over the empty blocks,

$$
S(n,k) = \frac{1}{k!} \sum_{j=0}^{k} (-1)^{k-j} \binom{k}{j} j^n,
$$

a sum of $k + 1$ terms of about $n \log k$ bits. Two observations reduce
the work. The powers $j^n$ are needed only for odd $j$: for $j = 2^e j'$
with $j'$ odd, $j^n = j'^n \cdot 2^{en}$, a shift. And the binomial
coefficients are symmetric, so only half of them are computed. Modulo a
word-size prime the powers of composite $j$ are products of the powers of
two factors, read from a divisor table, so that only $\pi(k)$
exponentiations are needed; the sum is then a dot product.

Isolated values use:

- the triangular recurrence $S(n,k) = k S(n-1, k) + S(n-1, k-1)$ in one
  or two words while $S(n,k)$ fits, that is for $n \le 43$ on 64-bit
  machines;
- the exact power sum for $k \le c_1 n$, where the terms are small;
- the exponential generating function
  $$
  S(n,k) = \frac{n!}{k!} [x^{n-k}] \left(\frac{e^x - 1}{x}\right)^k
  $$
  for $k \ge c_2 n$, a power of a series of length $n - k + 1$;
- the power sum modulo $\bits - 1$ bit primes followed by Chinese
  remaindering in between.

The cutoffs are $c_1 = c_2 = 0.9$ for $n < 200$; for larger $n$, $c_1$ is
$0.95 e^{-0.00022 n}$ (for $n < 3000$) or $1500/n$, kept in $[0.0002,
0.8]$, and $c_2 = \min(0.92 + 0.005 \log n, 0.98)$. The number of primes
is determined by the bound
$$
S(n,k) \le \frac{1}{2} \binom{n}{k} k^{n-k} \qquad (1 \le k < n)
$$
of Rennie and Dobson [@RenDob1969].

A row $S(n, 0), \ldots, S(n, \mathit{klen} - 1)$ is a convolution of the
sequences $j^n/j!$ and $(-1)^i/i!$:
$$
S(n,k) = \sum_{j=0}^{k} \frac{j^n}{j!} \cdot \frac{(-1)^{k-j}}{(k-j)!}.
$$
FLINT uses the triangular recurrence for $n \le 80$, the convolution over
$\Z$ with the common denominator $\ell!^2$, $\ell = \min(\mathit{klen},
n) - 1$, for $\mathit{klen} < n/2$, and otherwise the convolution modulo word-size
primes with Chinese remaindering, where each entry uses only as many
primes as its own bound requires.

## Bell numbers {#bell}

Summing the power sum formula over $k$ and exchanging the sums gives
$$
B_n = \sum_{k=0}^{n} S(n, k) = \sum_{j=0}^n \frac{j^n}{j!}
\sum_{i=0}^{n-j} \frac{(-1)^i}{i!},
$$
and letting the inner sum run to infinity gives Dobiński's formula
[@Dob1877]
$$
B_n = \frac{1}{e} \sum_{k=0}^{\infty} \frac{k^n}{k!}.
$$

**Dobiński's formula, exactly.** With $P = N! \sum_{k=0}^N k^n/k!$ and $Q =
N! \sum_{k=0}^N 1/k!$, both integers,
$$
B_n - \frac{P}{Q} = \frac{T_P - B_n T_Q}{e - T_Q}, \qquad
T_P = \sum_{k > N} \frac{k^n}{k!},\ T_Q = \sum_{k > N} \frac{1}{k!}.
$$
For $N \ge n$ each $k > N$ has $k^n > n^n > B_n$, so $P/Q \le B_n$; and $N
\approx n(1 + 1.2/\log n)$ makes the tail $T_P$ small enough that $P/Q >
B_n - 1$. Hence $B_n = \lceil P/Q \rceil$, computed with integer
arithmetic only (`arith_bell_number_dobinski`). The terms are
accumulated by Horner's rule in $k$, so the cost is $O(N)$ multiplications
of a large number by a word plus the $N/2$ powers $k^n$ for odd $k$, about
$O(n^2 \log^2 n)$ bit operations.

**Multimodular sum.** Modulo a prime $p > n$, the double sum above is
evaluated with $O(n)$ multiplications: the factorials are precomputed as
$n!/i!$, so that the inner sums are prefix sums and only one inversion
(of $n!^2$) is needed at the end, and the powers $j^n$ of composite $j$
come from a divisor table. FLINT accumulates the products in three words
and reduces once, uses Montgomery arithmetic
([](algorithms/arithmetic/montgomery)) with primes of $\bits - 3$ bits,
and processes the primes in parallel (`arith_bell_number_multi_mod`).
The number of primes comes from de Bruijn's asymptotic expansion
[@deB1958]
$$
\log B_n = n \left( \log n - \log \log n - 1 + \frac{\log\log n}{\log n}
+ \frac{1}{\log n} + \frac{1}{2} \left(\frac{\log\log n}{\log n}\right)^2
+ O\left(\frac{\log\log n}{\log^2 n}\right) \right).
$$

**Vectors.** The Bell triangle (Aitken's array) $a_{0,0} = 1$, $a_{i,0} =
a_{i-1,i-1}$, $a_{i,j} = a_{i,j-1} + a_{i-1,j-1}$ has $a_{i,0} = B_i$
and needs $n^2/2$ additions for $B_0, \ldots, B_{n-1}$. Modulo $m$ there
are two power series methods: the exponential generating function
$\exp(e^x - 1)$, computed with one power series exponential in $O(M(n))$
operations if $(n-1)!$ is invertible modulo $m$, and the ordinary
generating function
$$
\sum_{n \ge 0} B_n x^n = \sum_{k \ge 0} \frac{x^k}{(1 - x)(1 - 2x)
\cdots (1 - kx)},
$$
which follows from $\sum_n S(n,k) x^n = x^k / \prod_{j=1}^k (1 - jx)$.
Truncated at $k = n - 1$, the sum is a fraction $R/Q$ that binary
splitting computes as polynomials of degree $n$, followed by one power
series division: $O(M(n) \log n)$ operations, for any modulus.

## Use in FLINT {#use}

| Function | Method |
|---|---|
| `arith_stirling_number_1u`, `arith_stirling_number_1` | special cases; truncated binary splitting; EGF for $n > 140$, $k > 0.87 n$ |
| `arith_stirling_number_2` | special cases; triangle for $n \le 43$; power sum, EGF or multimodular power sum by $k/n$ |
| `arith_stirling_number_1u_vec`, `arith_stirling_number_1_vec` | binary splitting of the rising factorial |
| `arith_stirling_number_2_vec` | triangle for $n \le 80$; convolution over $\Z$ for $\mathit{klen} < n/2$; multimodular convolution |
| `arith_stirling_matrix_1u`, `arith_stirling_matrix_2`, ... | `gr_mat_stirling`: triangular recurrences row by row |
| `arith_bell_number` | table ($n \le 25$); triangle ($n \le 58$); Dobiński ($n \le 3400$); multimodular |
| `arith_bell_number_vec` | triangle for $n < 5000$; multimodular vector |
| `arith_bell_number_nmod` | table; periodicity for $m = 2, 3$; the $O(n)$ sum if $m > n$ and $n!$ is invertible; else the vector |
| `arith_bell_number_nmod_vec` | triangle for $n < 300$; EGF series if $m \ge n$ and invertible; triangle for $n < 500 + \mathit{bits}(m)^2$; OGF |

The generic ring functions `gr_stirling_s1u_uiui`, `gr_stirling_s2_uiui`,
`gr_stirling_s1_ui_vec` (and the other kinds), `gr_bellnum_ui` and
`gr_bellnum_vec` call these functions over $\Z$ (and
`arith_bell_number_nmod` and `arith_bell_number_nmod_vec` over $\Z/n\Z$)
and convert the result for other rings.

Bell numbers modulo 2 and 3 are periodic with periods 3 and 13, by
Touchard's congruence $B_{n+p} \equiv B_n + B_{n+1} \pmod p$ for prime
$p$ [@Comtet1974].

::: functions arith_stirling_number_2 arith_bell_number arith_bell_number_nmod_vec :::
