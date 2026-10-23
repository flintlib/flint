---
title: The p-adic exponential and logarithm
summary: Convergence and truncation bounds for the p-adic exponential and logarithm series, their evaluation by rectangular splitting with deferred denominators, the balanced (bit-burst) method with binary splitting, and the Satoh–Skjernaa–Taguchi reduction of the logarithm, as implemented in padic, qadic and padic_radix.
complexity: For an argument of valuation $v$ to precision $N$: rectangular splitting $O(\sqrt{N/v})$ multiplications of $O(N \log p)$-bit integers plus $O(N/v)$ cheap operations; balanced method $O(M(N(\log p + \log N)) \log^2 N)$ bit operations.
order: 30
---

The exponential and the logarithm of $p$-adic numbers are defined by the
same power series as over the reals,

$$
\exp(x) = \sum_{i \ge 0} \frac{x^i}{i!}, \qquad
\log(1 - y) = -\sum_{i \ge 1} \frac{y^i}{i},
$$

but the series converge only on small discs, and the denominators
$i!$ and $i$ are divisible by $p$ for large $i$. The same algorithms are used in
$\mathbf{Q}_p$ (`padic`, `padic_radix`) and in the unramified
extensions $\mathbf{Q}_q$ (`qadic`), where the arguments are polynomials
reduced modulo the defining polynomial but the coefficients of the
series are the same rational numbers.

## Convergence and truncation {#bounds}

Legendre's formula gives
$\operatorname{ord}_p(i!) = (i - s_p(i))/(p - 1)$, where $s_p(i) \ge 1$
is the sum of the base-$p$ digits of $i \ge 1$. Hence
$\operatorname{ord}_p(i!) \le (i - 1)/(p - 1)$, with equality when $i$
is a power of $p$.

::: lemma Terms of the exponential {#lem-exp-terms}
If $\operatorname{ord}_p(x) = v$ then
$\operatorname{ord}_p(x^i / i!) \ge iv - (i - 1)/(p - 1)$, with equality
for $i = p^k$. The series $\exp(x)$ converges if and only if
$v > 1/(p - 1)$, that is, $v \ge 1$ for odd $p$ and $v \ge 2$ for
$p = 2$. In that case every term with $i \ge 1$ has valuation at least
$v$, so $\exp(x) \in 1 + x\mathbf{Z}_p$ is a unit.
:::

::: proof
The bound follows from Legendre's formula. If $v > 1/(p - 1)$, the
right-hand side $i(v - 1/(p - 1)) + 1/(p - 1)$ tends to infinity, and it
is at least $v$ for $i \ge 1$ because it is increasing in $i$. If
$v \le 1/(p - 1)$, the terms with $i = p^k$ have valuation
$p^k(v - 1/(p - 1)) + 1/(p - 1) \le 1/(p - 1) < 1$, so infinitely many
terms are not divisible by $p$ and the series diverges.
:::

::: theorem Truncation of the exponential {#thm-exp-bound}
Let $1/(p - 1) < v < N$ and

$$
n = \left\lceil \frac{(p - 1)N - 1}{(p - 1)v - 1} \right\rceil .
$$

Then $\operatorname{ord}_p(x^i/i!) \ge N$ for all $i \ge n$, so
$\exp(x) \equiv \sum_{i < n} x^i / i! \pmod{p^N}$.
:::

::: proof
For $i \ge n$ we have $i((p - 1)v - 1) \ge (p - 1)N - 1$, which is
equivalent to $iv - (i - 1)/(p - 1) \ge N$; apply
[](#lem-exp-terms).
:::

`_padic_exp_bound` returns this $n$ when $p$ fits in a `slong`. For
larger $p$ it returns $\lceil N/v \rceil$, which is correct because then
$n \le N < p$ and $i!$ is a unit for all $i < p$.

::: theorem Truncation of the logarithm {#thm-log-bound}
Let $\operatorname{ord}_p(y) = v \ge 1$. Then
$\operatorname{ord}_p(y^i / i) \ge iv - \lfloor \log_p i \rfloor$, the
series $\log(1 - y)$ converges, and $\operatorname{ord}_p(y^i/i) \ge N$
for every $i \ge b$ whenever $bv - \lceil \log_p b \rceil \ge N$ and
$b \ge 2$.
:::

::: proof
The first statement holds because $p^{\operatorname{ord}_p(i)} \le i$.
The function $g(i) = iv - \lceil \log_p i \rceil$ satisfies
$g(i + 1) - g(i) \ge v - 1 \ge 0$, since $\lceil \log_p i \rceil$
increases by at most one from $i$ to $i + 1$; hence $g(i) \ge g(b) \ge
N$ for $i \ge b$, and $iv - \operatorname{ord}_p(i) \ge g(i)$.
:::

`_padic_log_bound` starts from the estimate
$b = \lceil (c + \lceil \log_p c \rceil + 1)/v \rceil$ with
$c = N - \lfloor \log_p v \rfloor$ and decreases $b$ while
$bv - \lceil \log_p b \rceil \ge N$ still holds, returning the least
such $b \ge 2$. The logarithm converges for $v \ge 1$ for every $p$,
including $p = 2$; then $\log(1 - y)$ has valuation at least $v$.

## Rectangular splitting {#rectangular}

A polynomial $\sum_{i < n} c_i x^i$ can be evaluated with
$m + n/m$ multiplications of full-size elements, by computing the powers
$x^2, \ldots, x^m$ once and writing the sum as
$\sum_j \bigl(\sum_{i < m} c_{jm + i} x^i\bigr) x^{jm}$, evaluated by
Horner's rule in $x^m$; the inner sums need only multiplications by the
scalars $c_i$ [@PS1973]. With $m \approx \sqrt n$ this is about
$2\sqrt n$ multiplications.

For the exponential, the coefficients $1/i!$ are not $p$-adic integers,
so FLINT multiplies through by $(n - 1)!$: the coefficients
$(n - 1)!/i!$ are integers, built block by block as products of small
integers, and the sum $S = \sum_{i < n} ((n-1)!/i!) x^i$ is computed
modulo $p^{N + k}$ with $k = \lfloor (n - 2)/(p - 1) \rfloor \ge
\operatorname{ord}_p((n - 1)!)$. Since $\exp(x)$ is a unit,
$\operatorname{ord}_p(S) = \operatorname{ord}_p((n - 1)!)$; the factors
of $p$ are removed from $S$ and $(n-1)!$, and the unit part of
$(n - 1)!$ is inverted modulo $p^N$ ([](algorithms/p-adic/hensel-lifting)).
The guard digits $k$ guarantee that $S / (n-1)!$ is correct modulo
$p^N$. The logarithm is handled in the same way, with the block
denominators $\prod_{i} (jb + i)$ in place of $(n - 1)!$ and
$k = \lfloor \log_p n \rfloor$ guard digits; each block is scaled by its
own inverted denominator times $p^k$, and the total is divided by $p^k$
at the end.

The cost is $O(\sqrt{n})$ multiplications of integers of
$O((N + k) \log p)$ bits, where $n \approx N/v$, plus $O(n)$
multiplications by single-word integers. Small arguments ($n \le 3$)
are handled by explicit formulas.

## The balanced method {#balanced}

Brent's bit-burst method [@Bre1976b] evaluates $\exp(x)$ for a real $x$
by splitting $x$ into pieces with exponentially growing numbers of bits
and handling each piece by binary splitting. FLINT uses the $p$-adic
analogue. Write the integer representative of $x$ modulo $p^N$ in base
$p$ and cut it at the digit positions $2, 4, 8, \ldots$:

$$
x = x_0 + x_1 + x_2 + \cdots, \qquad
x_j \in p^{w_j}\Z, \quad 0 \le x_j < p^{2 w_j},
\quad w_j = 2^j,
$$

where $x_0$ consists of the digits below $p^2$. Then
$\exp(x) = \prod_j \exp(x_j)$, and each factor is a series with
argument of valuation at least $w_j$, so by [](#thm-exp-bound) about
$N/w_j$ terms are needed, while the argument has only $2w_j$ digits.

::: algorithm Binary splitting for a short exponential series {#alg-bsplit}
**Input:** an integer $x$ and $1 \le a < b$.

**Output:** integers $P = x^{b - a}$, $Q = (b-1)!/(a-1)!$ and
$T$ with $T/Q = \sum_{i=a}^{b-1} x^{i - a + 1}(a-1)!/i!$.

1. If $b - a \le 2$, return the values directly.
2. Otherwise let $m = \lfloor (a + b)/2 \rfloor$, compute
   $(P_1, Q_1, T_1)$ for $[a, m)$ and $(P_2, Q_2, T_2)$ for $[m, b)$
   recursively, and return $P_1 P_2$, $Q_1 Q_2$ and
   $T_1 Q_2 + P_1 T_2$.
:::

The values are exact integers. At the top level, with $a = 1$,
$\exp(x_j) = 1 + T/Q$, and since $\exp(x_j)$ is a unit,
$\operatorname{ord}_p(T + Q) = \operatorname{ord}_p(Q)$: the common
power of $p$ is removed and the unit part of $Q$ inverted modulo $p^N$.
For the piece $x_j$ with $n_j \approx N/w_j$ terms the integers at the
root have $O(n_j \cdot w_j \log p + n_j \log n_j) = O(N \log p +
n_j \log n_j)$ bits, so the binary splitting costs
$O(M(N(\log p + \log N)) \log N)$, and the $O(\log N)$ pieces together
cost $O(M(N(\log p + \log N)) \log^2 N)$. For $p = 2$ the first piece
is zero, since $x$ has valuation at least 2.

For the logarithm a product decomposition is needed instead. With
$x = 1 - y$, write $y = r_0 + s$ where $r_0$ holds the digits of $y$
below $p^2$. Then

$$
1 - y = (1 - r_0)\Bigl(1 - \frac{s}{1 - r_0}\Bigr),
$$

and $t = s/(1 - r_0)$ has valuation at least 2. Repeating the step on
$t$ with the cut at $p^4$, then $p^8$, and so on, gives
$\log(1 - y) = \sum_j \log(1 - r_j)$ with $r_j$ of valuation at least
$2^j$ and fewer than $2^{j+1}$ digits. Each logarithm is summed by
binary splitting of $\sum_{i \ge 1} r_j^i / i$ (with $Q$ the product
of the indices), and the divisions by $1 - r_j$ are inversions modulo
$p^N$ that reuse one set of precomputed powers (`_padic_inv_precomp`).
The cost is of the same order as for the exponential.

## The Satoh–Skjernaa–Taguchi reduction {#satoh}

::: lemma {#lem-sst}
If $\operatorname{ord}_p(a - 1) \ge 1$ then
$\operatorname{ord}_p(a^p - 1) \ge \operatorname{ord}_p(a - 1) + 1$;
hence $\operatorname{ord}_p(a^{p^k} - 1) \ge k + 1$.
:::

::: proof
Write $a = 1 + h$ with $\operatorname{ord}_p(h) = w \ge 1$. Then
$a^p - 1 = ph + \binom{p}{2}h^2 + \cdots + h^p$. The middle binomial
coefficients are divisible by $p$, so those terms have valuation at
least $1 + 2w$, and $h^p$ has valuation $pw \ge w + 1$.
:::

Since $\log(a^{p^k}) = p^k \log(a)$, the logarithm of $a$ modulo $p^N$
is $p^{-k}$ times the logarithm of $a^{p^k}$ modulo $p^{N + k}$, and by
[](#lem-sst) the series for the latter has argument of valuation at
least $k + 1$, so only about $(N + k)/(k + 1)$ terms. Satoh, Skjernaa
and Taguchi used this to compute logarithms in point counting
[@SST2003]. `_padic_log_satoh` takes $k = \lfloor \sqrt N \rfloor$,
computes $a^{p^k}$ modulo $p^{N + k}$ by `fmpz_powm`, evaluates its
logarithm by rectangular splitting and divides by $p^k$ exactly; for
$N < 16$ it calls `_padic_log_rectangular` directly.

## Use in FLINT {#use}

| Function | Method | Selected by default |
|---|---|---|
| `_padic_exp_rectangular`, `padic_exp_rectangular` | rectangular splitting | `_padic_exp` for $N < 1024$ |
| `_padic_exp_balanced`, `padic_exp_balanced` | bit-burst with binary splitting | `_padic_exp` for $N \ge 1024$ |
| `_padic_log_rectangular`, `padic_log_rectangular` | rectangular splitting | `_padic_log` for $N < 512/\lceil\log_2 (p+1)\rceil$ |
| `_padic_log_balanced`, `padic_log_balanced` | product decomposition with binary splitting | `_padic_log` otherwise |
| `_padic_log_satoh`, `padic_log_satoh` | power $a^{p^k}$, then rectangular splitting | never |
| `_qadic_exp_rectangular`, `_qadic_exp_balanced` | as for `padic`, in $\mathbf{Z}_q$ | `_qadic_exp` for $N < 8192/\lceil\log_2 (p+1)\rceil$, resp. otherwise |
| `_qadic_log_rectangular`, `_qadic_log_balanced` | as for `padic`, in $\mathbf{Z}_q$ | `_qadic_log` for $N < 1024/\lceil\log_2 (p+1)\rceil$, resp. otherwise |
| `_padic_radix_exp_rectangular`, `_padic_radix_exp_balanced` | as above, in radix arithmetic | rectangular for $N$ below 80 to 200, depending on the size of $p$ |
| `_padic_radix_log_rectangular`, `_padic_radix_log_balanced` | as above, in radix arithmetic | rectangular for $N$ below 80 to 350, depending on the size of $p$ |

`qadic_norm_analytic` computes norms as $\exp(\operatorname{Tr}
\log(x))$ with these functions
([](algorithms/p-adic/unramified-extensions#norm)).

The radix versions in `padic_radix` refine the methods: the rectangular
splitting works in fixed point on limbs of radix $p^e$ and defers all
denominators to a single inversion; the balanced exponential tracks the
factor $p^{w_j}$ of each piece as a digit shift, so that the series is
evaluated with a working precision of about $N - w_j$ digits, and uses a
precomputed table of powers in the binary splitting tree (the scheme of
`_arb_exp_sum_bs_powtab`).
