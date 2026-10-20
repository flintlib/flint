---
title: arith.h
header: arith.h
description: Number-theoretic and combinatorial sequences — Bernoulli, Euler, Stirling and Bell numbers, partitions, sums of squares, Ramanujan's tau function, divisors and Landau's function.
category: Number theory
order: 1
---

This module computes integer and rational sequences of number theory and
combinatorics: Bernoulli and Euler numbers and polynomials, Stirling
numbers of both kinds, Bell numbers (also modulo a word-size integer),
the partition function, the number of representations as a sum of
squares, Ramanujan's tau function, divisor lists and Landau's function.
Most functions come in a version for an isolated value and a version
for a whole vector, which use different algorithms: a single large value
is typically computed from an analytic formula or a multimodular sum,
and a vector by a recurrence or by power series arithmetic, often modulo
several word-size primes followed by Chinese remaindering.

Several functions are thin wrappers around specialised modules:
Bernoulli numbers come from `bernoulli`, the partition function from
`partitions`, Euler numbers from `arb`, the Stirling matrices from
`gr_mat`. The algorithms are described in
[](algorithms/number-theory/bernoulli),
[](algorithms/number-theory/partitions) and
[](algorithms/number-theory/stirling-bell), and for the combinatorial
numbers in [](arith/stirling) and [](arith/bell). Functions for primorials, divisor sums, the Möbius
and Euler functions and the classical orthogonal and cyclotomic
polynomials, formerly in this module, are now in `fmpz`, `fmpz_poly` and
`fmpq_poly`; the old names remain as macros (see
[Deprecated aliases](#deprecated)).

## Harmonic numbers {#harmonic}

::: function
void arith_harmonic_number(fmpq_t x, slong n)
void _arith_harmonic_number(fmpz_t num, fmpz_t den, slong n)

Set the result to the harmonic number $H_n = \sum_{k=1}^n 1/k$, or to 0
for $n \le 0$. Wrappers of `fmpq_harmonic_ui`.
:::

## Bernoulli numbers and polynomials {#bernoulli}

The Bernoulli numbers are defined by $x/(e^x - 1) = \sum_{n \ge 0} B_n
x^n/n!$, so that $B_1 = -1/2$ and $B_n = 0$ for odd $n \ge 3$. The
`bernoulli` module provides the underlying algorithms and a cache.

::: summary algorithms/number-theory/bernoulli :::

::: function
void arith_bernoulli_number(fmpq_t x, ulong n)
void _arith_bernoulli_number(fmpz_t num, fmpz_t den, ulong n)

Set the result to the Bernoulli number $B_n$, in lowest terms with a
positive denominator. Wrappers of `bernoulli_fmpq_ui`: the value is read
from the Bernoulli cache if present and otherwise computed from the zeta
function or by the multimodular algorithm.
:::

::: function
void arith_bernoulli_number_vec(fmpq * x, slong n)
void _arith_bernoulli_number_vec(fmpz * num, fmpz * den, slong n)

Set the entries to $B_0, B_1, \ldots, B_{n-1}$. For $n < 700$ the function
calls `_arith_bernoulli_number_vec_recursive`, otherwise
`_arith_bernoulli_number_vec_multi_mod`. For long vectors,
`bernoulli_fmpq_vec_no_cache` (which uses the zeta function and several
threads) is usually faster.
:::

::: function
void _arith_bernoulli_number_vec_recursive(fmpz * num, fmpz * den, slong n)

Sets the numerators and denominators of $B_0, \ldots, B_{n-1}$ using
Ramanujan's recurrence, which expresses $B_m$ in terms of the $B_k$ with
$k \equiv m \pmod 6$. Intended for $n$ up to a few thousand.

::: details
The first values are taken from a table. The recurrence is then run with
all numerators over the common denominator $(n+1)\#$ (the primorial), so
that only integer additions and multiplications by binomial coefficients
are needed; the binomial coefficients are updated with products of up to
six small factors per word. The cost is $O(n^2)$ operations on numbers of
$O(n \log n)$ bits. The function throws for $n$ beyond about $3 \cdot
10^9$ (64-bit), where the word-size products would overflow.
:::
:::

::: function
void _arith_bernoulli_number_vec_multi_mod(fmpz * num, fmpz * den, slong n)

Sets the numerators and denominators of $B_0, \ldots, B_{n-1}$ by power
series inversion modulo word-size primes and Chinese remaindering.

::: details
The function inverts
$$
\frac{\cosh(x) - 1}{x^2} = \sum_{k \ge 0} \frac{x^{2k}}{(2k+2)!},
\qquad
\frac{x^2}{\cosh(x) - 1} = \sum_{k=0}^{\infty} \frac{(2-4k) B_{2k}}{(2k)!}
x^{2k}
$$
modulo enough primes of $\bits - 1$ bits to determine the numerators
(bounded by `arith_bernoulli_number_size` plus the size of the
denominator), which saves half of the work compared with $x/(e^x-1)$ since
the odd terms vanish [@BuhlerCrandallSompolski1992]. The numerators are
recovered by Chinese remaindering with a comb that uses only as many
primes as each entry needs. The cost is $O(n^2 \log^{O(1)} n)$.
:::
:::

::: function
void arith_bernoulli_number_denom(fmpz_t den, ulong n)

Sets `den` to the denominator of $B_n$ in lowest terms: 1 for odd $n \ge
3$ and $n = 0$, 2 for $n = 1$, and for even $n \ge 2$ the product of the
primes $p$ with $(p - 1) \mid n$, by the von Staudt–Clausen theorem. Values
for $n \le 178$ are read from a table.
:::

::: function
double arith_bernoulli_number_size(ulong n)

Returns $b$ such that $|B_n| < 2^{\lfloor b \rfloor}$, from $|B_n| < 4
\cdot n!/(2\pi)^n$ and $n! \le (n+1)^{n+1} e^{-n}$. Odd $n$ are not
treated specially. The bound is evaluated in double precision and is not
guaranteed for $n > 10^{14}$.
:::

::: function
void arith_bernoulli_polynomial(fmpq_poly_t poly, ulong n)

Sets `poly` to the Bernoulli polynomial $B_n(x) = \sum_{k=0}^n
\binom{n}{k} B_k x^{n-k}$, computing $B_0, \ldots, B_n$ with
`_arith_bernoulli_number_vec` and scaling by binomial coefficients over a
common denominator.
:::

## Euler numbers and polynomials {#euler}

The Euler numbers are the integers $E_n$ defined by $1/\cosh(t) =
\sum_{n \ge 0} E_n t^n/n!$. With this convention the odd-indexed numbers
are zero and the even ones alternate in sign: $E_0, E_1, E_2, \ldots = 1,
0, -1, 0, 5, 0, -61, 0, 1385, \ldots$ The Euler polynomials are defined by
$2e^{xt}/(e^t + 1) = \sum_{n \ge 0} E_n(x) t^n/n!$, and $E_n =
2^n E_n(1/2)$.

::: function
void arith_euler_number(fmpz_t res, ulong n)

Sets `res` to $E_n$. A wrapper of `arb_fmpz_euler_number_ui`, which uses
a table for small $n$, the Dirichlet beta function $\beta(n+1)$ for $n <
2000$ and a combination of the beta function with a multimodular method
(as for Bernoulli numbers) for larger $n$.
:::

::: function
void arith_euler_number_vec(fmpz * res, slong n)

Sets `res` to $E_0, E_1, \ldots, E_{n-1}$, for $n \ge 0$.

::: details
The even-indexed $|E_{2k}|$ are computed by inverting the power series of
$\cosh$ modulo primes of $\bits - 1$ bits, as many as
`arith_euler_number_size` requires, and recovered by Chinese
remaindering, using for each entry only as many primes as it needs.
:::
:::

::: function
double arith_euler_number_size(ulong n)

Returns $b$ such that $|E_n| < 2^{\lfloor b \rfloor}$, from $|E_n| <
2^{n+2} n!/\pi^{n+1}$ and $n! \le (n+1)^{n+1} e^{-n}$. Odd $n$ are not
treated specially. Not guaranteed for $n > 10^{14}$.
:::

::: function
void arith_euler_polynomial(fmpq_poly_t poly, ulong n)

Sets `poly` to the Euler polynomial $E_n(x)$, using
$$
E_n(x) = \frac{2}{n+1}\left(B_{n+1}(x) -
2^{n+1}B_{n+1}\left(\frac{x}{2}\right)\right)
$$
with $B_{n+1}(x)$ computed once by `arith_bernoulli_polynomial` and
rescaled.
:::

## Stirling numbers {#stirling}

The Stirling numbers of the first kind $s(n,k)$ (signed) and $|s(n,k)|$
(unsigned) and of the second kind $S(n,k)$ are defined by

$$
x(x-1)\cdots(x-n+1) = \sum_{k=0}^n s(n,k)\, x^k, \qquad
x(x+1)\cdots(x+n-1) = \sum_{k=0}^n |s(n,k)|\, x^k, \qquad
x^n = \sum_{k=0}^n S(n,k)\, x(x-1)\cdots(x-k+1),
$$

and are zero for $k > n$. $|s(n,k)|$ counts the permutations of $n$
elements with $k$ cycles and $S(n,k)$ the partitions of a set of $n$
elements into $k$ blocks. See [](arith/stirling).

::: function
void arith_stirling_number_1u(fmpz_t s, ulong n, ulong k)
void arith_stirling_number_1(fmpz_t s, ulong n, ulong k)
void arith_stirling_number_2(fmpz_t s, ulong n, ulong k)

Set $s$ to $|s(n,k)|$, $s(n,k)$ or $S(n,k)$, for any $n, k \ge 0$. These
functions compute an isolated value, with algorithms chosen according to
$n$ and $k/n$ (binary splitting, power sums, power series and
multimodular methods).
:::

::: function
void arith_stirling_number_1u_vec(fmpz * row, ulong n, slong klen)
void arith_stirling_number_1_vec(fmpz * row, ulong n, slong klen)
void arith_stirling_number_2_vec(fmpz * row, ulong n, slong klen)

Set `row` to the values for $k = 0, 1, \ldots, \mathit{klen} - 1$ of
$|s(n,k)|$, $s(n,k)$ or $S(n,k)$. A full row has $\mathit{klen} = n + 1$;
entries with $k > n$ are set to zero.
:::

::: function
void arith_stirling_matrix_1u(fmpz_mat_t mat)
void arith_stirling_matrix_1(fmpz_mat_t mat)
void arith_stirling_matrix_2(fmpz_mat_t mat)

Set the entry $(n, k)$ of the $r \times c$ matrix `mat` to $|s(n,k)|$,
$s(n,k)$ or $S(n,k)$, for $0 \le n < r$ and $0 \le k < c$; the entries
above the diagonal are zero. The square matrices of $s$ and $S$ are
inverses of each other. Wrappers of `gr_mat_stirling`, which fills the
rows with the triangular recurrences.
:::

## Bell numbers {#bell}

The Bell number $B_n$ is the number of partitions of a set of $n$
elements, $B_n = \sum_{k=0}^n S(n,k)$, with exponential generating
function $\sum B_n x^n/n! = e^{e^x - 1}$. See [](arith/bell).

::: function
void arith_bell_number(fmpz_t b, ulong n)
void arith_bell_number_dobinski(fmpz_t res, ulong n)
void arith_bell_number_multi_mod(fmpz_t res, ulong n)

Set the result to the Bell number $B_n$. The default version uses a
table for $n \le 25$ (64-bit), the Bell triangle in fixed-size multiword
arithmetic for $n \le 58$, Dobiński's formula for $n \le 3400$ and the
multimodular algorithm above.

::: details
The `dobinski` version computes $P = N! \sum_{k=0}^N k^n/k!$ and $Q = N!
\sum_{k=0}^N 1/k! \approx N!\, e$ for $N \approx n (1 + 1.2/\log n)$ and
returns $\lceil P/Q \rceil$, an exact form of Dobiński's formula $B_n =
e^{-1} \sum_{k \ge 0} k^n/k!$ [@Dob1877] without floating-point
arithmetic. Powers of even $k$ are obtained from those of odd $k$ by
shifts.

The `multi_mod` version computes $B_n$ modulo primes of $\bits - 3$ bits
with the sum used by `arith_bell_number_nmod`, in Montgomery arithmetic
and with powers $j^n$ for composite $j$ obtained as products from a
divisor table, and combines the residues by Chinese remaindering. The
number of primes is determined by `arith_bell_number_size`. The primes
are processed in parallel with the threads set by
`flint_set_num_threads`.
:::
:::

::: function
void arith_bell_number_vec(fmpz * b, slong n)
void arith_bell_number_vec_recursive(fmpz * b, slong n)
void arith_bell_number_vec_multi_mod(fmpz * b, slong n)

Set $b$ to $B_0, B_1, \ldots, B_{n-1}$. The `recursive` version uses the
Bell triangle (Aitken's array), with $O(n^2)$ additions of numbers of
$O(n \log n)$ bits; the `multi_mod` version computes the vector modulo
word-size primes with `arith_bell_number_nmod_vec` and combines them by
Chinese remaindering, in $O(n^2 \log^{O(1)} n)$. The default version uses
the recursive algorithm for $n < 5000$.
:::

::: function
ulong arith_bell_number_nmod(ulong n, nmod_t mod)

Returns $B_n$ modulo the modulus of `mod`, which can be any integer $\ge
1$.

::: details
Small $n$ are read from the table, and the moduli 2 and 3 use the periods
3 and 13 of the sequence. If $n$ is less than the modulus and $n!$ is
invertible, the function evaluates
$$
B_n = \sum_{k=0}^n \frac{(n-k)^n}{(n-k)!} \sum_{j=0}^k
\frac{(-1)^j}{j!}
$$
with the factorials scaled by $n!$ so that the main loop only multiplies,
the powers of composite numbers sieved from those of their factors, and
the sum accumulated in three words; this costs $O(n)$ multiplications and
$O(n)$ memory. Otherwise it computes $B_0, \ldots, B_n$ with
`arith_bell_number_nmod_vec` (and throws if $n > 2^{\bits - 3}$).
:::
:::

::: function
void arith_bell_number_nmod_vec(nn_ptr b, slong n, nmod_t mod)
void arith_bell_number_nmod_vec_recursive(nn_ptr b, slong n, nmod_t mod)
void arith_bell_number_nmod_vec_ogf(nn_ptr b, slong n, nmod_t mod)
int arith_bell_number_nmod_vec_series(nn_ptr b, slong n, nmod_t mod)

Set $b$ to $B_0, \ldots, B_{n-1}$ modulo the modulus of `mod`.

::: details
The `recursive` version uses the Bell triangle, with $O(n^2)$ additions.
The `ogf` version expands the ordinary generating function
$$
\sum_{n \ge 0} B_n x^n = \sum_{k \ge 0} \frac{x^k}{(1 - x)(1 - 2x) \cdots
(1 - kx)}
$$
truncated to length $n$, with the partial sums computed by binary
splitting as a fraction $R/Q$ and one power series division, in $O(M(n)
\log n)$. The `series` version computes $\exp(e^x - 1)$ as a power series
and multiplies by factorials, in $O(M(n))$; it requires $(n-1)!$ to be
invertible modulo the modulus and returns 0 (with $b$ unspecified) if it
is not, 1 on success. The `recursive` and `ogf` versions work for any
modulus. The default version uses the recursive algorithm for $n < 300$,
the series version if the modulus is at least $n$ and it succeeds, the
recursive algorithm for $n < 500 + \mathit{bits}^2$ where $\mathit{bits}$
is the bit size of the modulus, and the ogf version otherwise.
:::
:::

::: function
double arith_bell_number_size(ulong n)

Returns $b$ such that $B_n < 2^{\lfloor b \rfloor}$, from de Bruijn's
asymptotic expansion of $\log B_n$ [@deB1958] with a margin. The bound is
verified for $n < 2^{64}$ but not proved in general; an earlier version
used the weaker bound $B_n < (0.792n/\log(n+1))^n$ of [@BerTas2010].
:::

## The partition function {#partitions}

The partition function $p(n)$ counts the ways of writing $n$ as a sum of
positive integers without regard to order. Its generating function is
$\sum p(n) x^n = \prod_{k \ge 1} (1 - x^k)^{-1}$.

::: summary algorithms/number-theory/partitions :::

::: function
void arith_number_of_partitions(fmpz_t x, ulong n)

Sets $x$ to $p(n)$. A wrapper of `partitions_fmpz_ui`, which uses the
Hardy–Ramanujan–Rademacher formula for large $n$.
:::

::: function
void arith_number_of_partitions_vec(fmpz * res, slong len)
void arith_number_of_partitions_nmod_vec(nn_ptr res, slong len, nmod_t mod)

Set `res` to $p(0), p(1), \ldots, p(\mathit{len} - 1)$, exactly or modulo
the modulus of `mod`, by inverting Euler's pentagonal series
$\prod_{k \ge 1} (1 - x^k) = \sum_{j \in \Z} (-1)^j x^{j(3j-1)/2}$ with
`_fmpz_poly_inv_series` or `_nmod_poly_inv_series`. The cost is
$O(M(n))$ operations on numbers of $O(\sqrt n)$ bits, or modulo the
modulus.
:::

::: type
trig_prod_t

A product of cosines of rational multiples of $\pi$ times an algebraic
prefactor, $\mathit{prefactor} \cdot \sqrt{\mathit{sqrt\_p} /
\mathit{sqrt\_q}} \cdot \prod_{i < \mathit{n}} \cos(\pi\,
\mathit{cos\_p}[i] / \mathit{cos\_q}[i])$, with the fields named as in
the formula (at most `FLINT_BITS` factors). An array of length one of
`trig_prod_struct`.
:::

::: function
void trig_prod_init(trig_prod_t prod)

Initialises `prod` to the empty product with prefactor 1. Inline; no
memory is allocated.
:::

::: function
void arith_hrr_expsum_factored(trig_prod_t prod, ulong k, ulong n)

Multiplies `prod` (which must be initialised with `trig_prod_init`) by a
symbolic expression for the exponential sum
$$
A_k(n) = \sum_{\substack{0 \le h < k \\ \gcd(h, k) = 1}}
\exp\left(\pi i \left[ s(h,k) - \frac{2hn}{k}\right]\right)
$$
of the Hardy–Ramanujan–Rademacher formula, where $s(h,k)$ is a Dedekind
sum, as a product of cosines.

::: details
The sum is factored into sums for the prime powers dividing $k$ using
the multiplicativity relations of Lehmer and Whiteman [@Whiteman1956],
and each prime-power sum is a single cosine (or zero) determined by a
square root modulo the prime power. The prefactor is set to zero when
$A_k(n) = 0$. Since $A_k(n)$ depends only on $n \bmod k$, $n$ should be
reduced modulo $k$; the computation needs $24k$ and $24n$ to fit in a
word. See [](algorithms/number-theory/partitions).
:::
:::

## Sums of squares {#squares}

::: function
void arith_sum_of_squares(fmpz_t r, ulong k, const fmpz_t n)

Sets $r$ to the number $r_k(n)$ of representations of $n$ as an ordered
sum of $k$ squares of integers (signs and order count), with $r_k(0) = 1$
and $r_k(n) = 0$ for $n < 0$.

::: details
For $k = 1$ the result is 2 or 0. For $k = 2$ and $k = 4$ the function
factors $n$ and uses Jacobi's formulas $r_2(n) = 4 \sum_{d \mid n,\,
d \text{ odd}} (-1)^{(d-1)/2}$ and $r_4(n) = 8 \sum_{d \mid n,\, 4
\nmid d} d$ [@Grosswald1985]. For $k = 3$ and $k = 5$ it uses the
recursion $r_k(n) = \sum_{j \in \Z} r_{k-1}(n - j^2)$, with $O(\sqrt n)$
calls for $k - 1 \in \{2, 4\}$; $n$ is then assumed to fit in a `ulong`.
For other $k$ it computes the coefficients of $\vartheta_3(q)^k$ up to
$q^n$ with `arith_sum_of_squares_vec`, which requires $n$ to fit in a
`slong` and throws otherwise.
:::
:::

::: function
void arith_sum_of_squares_vec(fmpz * r, ulong k, slong n)

Sets $r_i$ to $r_k(i)$ for $0 \le i < n$, the first $n$ coefficients of
$$
\vartheta_3(q)^k = \Bigl(\sum_{j \in \Z} q^{j^2}\Bigr)^k,
$$
computed from the sparse expansions of $\vartheta_3$ and
$\vartheta_3^2$ and a truncated power.
:::

## Ramanujan's tau function {#tau}

::: function
void arith_ramanujan_tau(fmpz_t res, const fmpz_t n)

Sets `res` to $\tau(n)$, the coefficient of $q^n$ in $\Delta(q) = q
\prod_{k \ge 1} (1 - q^k)^{24}$, or to 0 for $n \le 0$.

::: details
The function factors $n$, computes $\tau(p)$ for the prime factors with
`arith_ramanujan_tau_series`, and uses the multiplicativity of $\tau$ and
the recurrence $\tau(p^{r+1}) = \tau(p) \tau(p^r) - p^{11}
\tau(p^{r-1})$. The cost is dominated by the series of length $P + 1$,
where $P$ is the largest prime factor of $n$, which must fit in a word.
:::
:::

::: function
void arith_ramanujan_tau_series(fmpz_poly_t res, slong n)

Sets `res` to the polynomial with coefficients $\tau(0), \tau(1),
\ldots, \tau(n-1)$, using Jacobi's identity
$$
\Delta(q) = q \Biggl( \sum_{k \ge 0} (-1)^k (2k+1) q^{k(k+1)/2}
\Biggr)^8,
$$
evaluated with three squarings of which the first, of a very sparse
series, is done directly.
:::

## Divisors {#divisors}

::: function
void arith_divisors(fmpz_poly_t res, const fmpz_t n)

Sets the coefficients of `res` to the positive divisors of $|n|$ in
increasing order (so `res` has length $\sigma_0(n)$), and to the zero
polynomial if $n = 0$. Uses a table for $|n| < \bits$ and otherwise
factors $n$ with `fmpz_factor` and generates the products of prime
powers.
:::

## Landau's function {#landau}

::: function
void arith_landau_function_vec(fmpz * res, slong len)

Sets `res` to the first `len` values $g(0), g(1), \ldots$ of Landau's
function, the largest order of an element of the symmetric group $S_n$,
equivalently the largest $\operatorname{lcm}(n_1, \ldots, n_r)$ with $n_1
+ \cdots + n_r \le n$.

::: details
Implements the basic algorithm of [@DelegliseNicolasZimmermann2009], a
dynamic programme over the prime powers $p^e \le n$ for primes $p$ up to
about $1.33\sqrt{n \log n}$. The running time is $O(n^{3/2} / \sqrt{\log
n})$ operations on integers of $O(\sqrt{n \log n})$ bits.
:::
:::

## Tables {#tables}

::: variable
const ulong bell_number_tab[]

The Bell numbers $B_n$ for $n <$ `BELL_NUMBER_TAB_SIZE` (26 on 64-bit
machines, 16 on 32-bit machines).
:::

::: variable
const unsigned int partitions_lookup[128]

The values $p(n)$ for $n < 128$, all below $2^{32}$.
:::

::: variable
const ulong euler_number_small[]
const slong _bernoulli_numer_small[]

The absolute values $|E_{2k}|$ for $2k <$ `SMALL_EULER_LIMIT` (25 on
64-bit machines) and the numerators of $B_{2k}$ for $2k <$
`BERNOULLI_SMALL_NUMER_LIMIT` (35 on 64-bit machines). Static arrays
defined in the header.
:::

::: macro
BELL_NUMBER_TAB_SIZE
SMALL_EULER_LIMIT
BERNOULLI_SMALL_NUMER_LIMIT

The sizes of the tables above.
:::

## Deprecated aliases {#deprecated}

These macros map old names of functions that have moved to other
modules. New code should use the targets.

::: macro
arith_primorial
arith_divisor_sigma
arith_moebius_mu
arith_euler_phi

Aliases of `fmpz_primorial`, `fmpz_divisor_sigma`, `fmpz_moebius_mu` and
`fmpz_euler_phi`.
:::

::: macro
arith_cyclotomic_polynomial
_arith_cyclotomic_polynomial
arith_chebyshev_t_polynomial
_arith_chebyshev_t_polynomial
arith_chebyshev_u_polynomial
_arith_chebyshev_u_polynomial
arith_swinnerton_dyer_polynomial
_arith_swinnerton_dyer_polynomial
arith_legendre_polynomial

Aliases of `fmpz_poly_cyclotomic`, `fmpz_poly_chebyshev_t`,
`fmpz_poly_chebyshev_u`, `fmpz_poly_swinnerton_dyer` (and their
underscore versions) and `fmpq_poly_legendre_p`.
:::

::: macro
arith_dedekind_sum
arith_dedekind_sum_naive

Aliases of `fmpq_dedekind_sum` and `fmpq_dedekind_sum_naive`.
:::
