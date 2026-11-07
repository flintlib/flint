---
title: Multiple zeta values and multiple polylogarithms
summary: Multiple zeta values by Akhilesh's double tails in fixed-point arithmetic with exact error counting, and in quasi-linear time by Hölder convolution and binary splitting; star values and Yamamoto's interpolation from shared subwords; multiple polylogarithms and colored multiple zeta values through Goncharov's iterated integrals, as power series at 0 or split at 1/2, by binary splitting with majorant truncation bounds.
order: 51
---

The multiple zeta values
$$
\zeta(s_1, \ldots, s_r) = \sum_{n_1 > n_2 > \cdots > n_r \ge 1} \frac{1}{n_1^{s_1} \cdots n_r^{s_r}}
\qquad (s_1 \ge 2)
$$
and the multiple polylogarithms
$$
\operatorname{Li}_{s_1, \ldots, s_r}(x_1, \ldots, x_r) = \sum_{n_1 > \cdots > n_r \ge 1} \frac{x_1^{n_1} \cdots x_r^{n_r}}{n_1^{s_1} \cdots n_r^{s_r}}
$$
are computed in [`acb_mzeta`](acb_mzeta) from their representations as
iterated integrals. For multiple zeta values FLINT uses Akhilesh's
double tails [@Akh2017], a recursion with $O(k^2)$ values per level and
about $P/2$ levels for weight $k$ and $P$ bits, carried out in
fixed-point arithmetic with exactly counted rounding errors; at very
high precision, Hölder convolution [@BBBL2001] reduces the value to
multiple polylogarithms at $\tfrac12$, which are evaluated by binary
splitting in quasi-linear time. General multiple polylogarithms are
written as Goncharov's iterated integrals [@Gon2001] and summed as power
series, directly or after splitting the path at $\tfrac12$.

## Words and iterated integrals {#words}

With the differential forms $\omega_0 = dt/t$ and $\omega_1 = dt/(1-t)$,
a word $w = e_1 \cdots e_k$ in the letters $x_0, x_1$ defines the
iterated integral
$$
\zeta(w) = \int_{1 > t_1 > \cdots > t_k > 0} \omega_{e_1}(t_1) \cdots \omega_{e_k}(t_k),
$$
which converges when $w$ is *admissible*: it starts with $x_0$ and ends
with $x_1$. Expanding $1/(1 - t)$ as a geometric series and integrating
term by term gives
$$
\zeta(s_1, \ldots, s_r) = \zeta(x_0^{s_1 - 1} x_1 x_0^{s_2 - 1} x_1 \cdots x_0^{s_r - 1} x_1),
$$
so the weight $k = s_1 + \cdots + s_r$ is the length of the word and the
depth $r$ is the number of letters $x_1$. The admissible words of weight
$k$ are $x_0 y x_1$ for all $2^{k-2}$ words $y$ of length $k - 2$;
`acb_mzeta_mzv_index` numbers them by the binary integer
$m = (1y)_2$, writing $x_0 = 0$ and $x_1 = 1$, so that the indices of
weight $k$ are $2^{k-2} \le m < 2^{k-1}$. The substitution
$t \mapsto 1 - t$ shows that $\zeta(w) = \zeta(w^*)$ for the *dual* word
$w^*$, obtained by reversing $w$ and exchanging the letters; it is
admissible when $w$ is.

## Double tails {#double-tails}

Akhilesh's double tails of a word $w$ (not necessarily admissible) are
$$
\zeta(w)_{m,n} = \int_{1 > t_1 > \cdots > t_k > 0} (1 - t_1)^m\, \omega_{e_1}(t_1) \cdots \omega_{e_k}(t_k)\, t_k^n,
\qquad m, n \ge 0,
$$
with $\zeta(w)_{0,0} = \zeta(w)$ for admissible $w$.

::: lemma Recurrences {#lem-tails}
For words $u$,
$$
\zeta(x_1 u)_{m,n} = \frac{\zeta(u)_{m,n}}{m}, \qquad
\zeta(u x_0)_{m,n} = \frac{\zeta(u)_{m,n}}{n}, \qquad
\zeta(\varnothing)_{m,n} = \frac{m!\, n!}{(m+n)!}
$$
(for $m \ge 1$, respectively $n \ge 1$), and
$$
\zeta(x_0 u x_1)_{m,n} = \zeta(x_0 u x_1)_{m+1,n+1}
+ \frac{\zeta(u x_1)_{m+1,n}}{m+1} + \frac{\zeta(x_0 u)_{m+1,n+1}}{n+1}.
$$
:::

::: proof
For the first identity, integrate out $t_1$:
$\int_{t_2}^1 (1 - t_1)^{m-1} dt_1 = (1 - t_2)^m / m$; the second is the
same with $t_k$ and $\int_0^{t_{k-1}} t_k^{n-1} dt_k = t_{k-1}^n / n$, and
the third is the beta integral. Writing $(1 - t_1)^m = (1 - t_1)^{m+1} + t_1 (1 - t_1)^m$,
$$
\zeta(x_0 v)_{m,n} = \zeta(x_0 v)_{m+1,n} + \frac{\zeta(v)_{m+1,n}}{m+1},
$$
since $t_1 \cdot dt_1/t_1 = dt_1$ and $\int_{t_2}^1 (1-t_1)^m dt_1 = (1 - t_2)^{m+1}/(m+1)$.
In the same way, from $t_k^n = t_k^{n+1} + (1 - t_k) t_k^n$,
$\zeta(v x_1)_{m,n} = \zeta(v x_1)_{m,n+1} + \zeta(v)_{m,n+1}/(n+1)$.
Applying the first relation to $v = u x_1$ and the second to
$\zeta(x_0 u x_1)_{m+1,n}$ gives the diagonal step.
:::

Every word reduces by these rules to shorter words or to the same word
with $(m, n)$ increased by $(1, 1)$, so all values are determined by
those at a large level $n = N$. Moreover a subword $w[i..j)$ of $w$
appears in the recursion only with $m - n = d_i$, the number of letters
$x_0$ in $w[0..i)$. The algorithm therefore stores, at level $n$, one
value $\zeta(u)_{n + d, n}$ per *state* $(u, d) = (w[i..j), d_i)$, at
most $O(k^2)$ values, and computes the levels $n = N, N-1, \ldots, 0$
with the values of level $N + 1$ replaced by zero; within a level the
states are processed by increasing length. Only two levels are kept in
memory. At level 0 only the states reachable from the requested words
are computed.

::: lemma Size of the tails {#lem-tail-size}
For a word $u$ with $a$ letters $x_1$ and $b$ letters $x_0$,
$\zeta(u)_{m,n} \le \max_{0 \le t \le 1} (1-t)^{m-a} t^{n-b}$
when $m \ge a$ and $n \ge b$. In particular, at level $n \ge k$ every
value is at most $4^{-(n - k)} (2n + 1)$.
:::

::: proof
On the simplex $t_1 \ge t_j \ge t_k$, so the density of a letter $x_1$
at position $j$ is $1/(1 - t_j) \le 1/(1 - t_1)$ and that of a letter
$x_0$ is $1/t_j \le 1/t_k$. The integrand is therefore at most
$(1-t_1)^{m-a} t_k^{n-b} \le (1-t_1)^{m-a} t_1^{n-b}$, and the simplex
has volume at most 1. At level $n$
the exponents $m - a$ and $n - b$ are at least $n - k$ and
$\max_t ((1-t) t)^{n-k} = 4^{-(n-k)}$; for the empty word,
$m!\, n!/(m+n)! \le 1/\binom{2n}{n} \le (2n+1)/4^n$.
:::

So $N \approx P/2 + O(\log P)$ levels suffice for $P$ bits, and in the
fixed-point format of about $2N$ fraction bits the values at level $n$
have only about $2(N - n)$ significant bits: the precision of each
level is tapered accordingly.

::: theorem Truncation {#thm-tail-trunc}
If every value dropped at level $N + 1$ is bounded by $G$ (for example
by [](#lem-tail-size)), the error of the values at level 0 caused by
replacing them by zero is at most $3 (N+1)(N+2)\, G$.
:::

::: proof
Let $\epsilon_n$ bound the errors of all values at level $n$, so
$\epsilon_{N+1} = G$. In the diagonal step at level $n \ge 1$, the
terms $\zeta(x_0 u x_1)_{m+1,n+1}$ and $\zeta(x_0 u)_{m+1,n+1}/(n+1)$
belong to level $n + 1$ and contribute at most
$(1 + \frac{1}{n+1}) \epsilon_{n+1}$, while $\zeta(u x_1)_{m+1,n}/(m+1)$
belongs to level $n$ (its offset is one larger) and contributes at most
$\epsilon_n / (n+1)$, since $m \ge n$. The other rules divide by
$m, n \ge 1$ and do not increase errors. Hence
$\epsilon_n \le (1 + \frac{1}{n+1}) \epsilon_{n+1} + \frac{\epsilon_n}{n+1}$,
that is $\epsilon_n \le \frac{n+2}{n} \epsilon_{n+1}$. At level 0 the
same count gives $E \le 2 \epsilon_1 + E/2$, so $E \le 4 \epsilon_1$, for
the states with offset $m \ge 1$, and then $2 \epsilon_1 + E \le 6 \epsilon_1$
for the words themselves ($m = 0$). So
$\epsilon_0 \le 6 \prod_{n=1}^{N} \frac{n+2}{n}\, G = 3 (N+1)(N+2)\, G$.
:::

**Arithmetic.** All quantities are nonnegative, and the recursion is
run in fixed point with $B$ fraction bits (one integer limb on top)
with truncating divisions by the small integers $m$, $n$, $m + 1$,
$n + 1$. Every computed value is then a lower bound for the exact one,
and its error is an integer number of ulps that is counted exactly
alongside the values (one ulp per truncating operation, propagated
through the same linear recurrences). The empty-word values are
carried as $4^n m!\, n!/(m+n)!$, whose step ratio is close to 1, and
shifted down when used. The cost is dominated by the divisions by one
limb: the two divisions of a diagonal step are merged into one by
$(m+1)(n+1)$, and several independent divisions (states of the same
length, whose children are shorter) are interleaved with precomputed
inverses [@GraMol2010] to hide the latency of the remainder chain. The
total cost is $O(k^2 P^2)$ bit operations.

## Several words at once {#batches}

Different words share subwords, and the recursion can carry the states
of several words together. This is used for:

- **Star values.** $\zeta^\star(s) = \sum_{n_1 \ge \cdots \ge n_r} n_1^{-s_1} \cdots n_r^{-s_r}$
  is the sum of $\zeta(p)$ over the $2^{r-1}$ indices $p$ obtained from
  $s$ by replacing some commas by plus signs.
- **Yamamoto's interpolation** [@Yam2013]
  $\zeta^t(s) = \sum_p t^{r - \operatorname{depth}(p)} \zeta(p)$ over the
  same indices: the coefficients $c_j = \sum_{\operatorname{depth}(p) = r - j} \zeta(p)$
  are computed together and the polynomial is evaluated at $t$ by
  Horner's rule, so for $t$ far from $[0, 1]$ the error is relative to
  the largest term.
- **All values of a weight.** `acb_mzeta_mzv_all` evaluates all
  $2^{k-2}$ admissible words of weight $k$ (or all weights up to $k$) in
  batches limited by memory, so that each state is computed once per
  batch, and forms star values by summing over the merges.

## Hölder convolution {#holder}

At very high precision the quadratic cost of the double tails is
replaced by a quasi-linear algorithm [@BBBL2001]. Splitting the
integration domain of $\zeta(w)$ at $\tfrac12$ (Chen's formula) and
substituting $t \mapsto 1 - t$ in the variables above $\tfrac12$ gives
$$
\zeta(w) = \sum_{j=0}^{k} I(w^*[k-j..k))\, I(w[j..k)),
\qquad
I(u) = \int_{\frac12 > t_1 > \cdots > t_L > 0} \omega_{u_1}(t_1) \cdots \omega_{u_L}(t_L).
$$
Every suffix of the admissible words $w$ and $w^*$ ends with $x_1$, so
each $I$ is a multiple polylogarithm at $\tfrac12$: for
$w = x_0^{c_1 - 1} x_1 \cdots x_0^{c_d - 1} x_1$ and a suffix starting
inside block $i$ with first entry $c'$,
$$
I = \operatorname{Li}_{c', c_{i+1}, \ldots, c_d}(\tfrac12)
= \sum_{n \ge 1} \frac{2^{-n}}{n^{c'}}\, H_{i+1}(n - 1),
\qquad
H_i(n) = \sum_{n \ge n_i > \cdots > n_d \ge 1} n_i^{-c_i} \cdots n_d^{-c_d}.
$$
With $G_i(n) = 2^{-n} H_i(n)$ and one accumulator $A_j$ per suffix,
$$
G_i(n) = \tfrac12 G_i(n-1) + \tfrac12 n^{-c_i} G_{i+1}(n-1), \qquad
A_j(n) = A_j(n-1) + \tfrac12 n^{-c'} G_{i+1}(n-1),
$$
a lower triangular linear recurrence $X(n) = M(n) X(n-1)$ of dimension
$d + k$. The matrices $2 n^{C} M(n)$, $C = \max_i c_i$, have integer
entries; their product over $n = 1, \ldots, N$ is computed by binary
splitting, exactly at the lower levels and in ball arithmetic at the
working precision above, the products for $w$ and $w^*$ in parallel.

::: proposition Truncation {#prop-holder-trunc}
Each $A_j$ omits at most
$2^{-N} (1 + \log(N+1))^d / (1 - q)$, $q = (1 + 1/N)^d / 2$, when the
products stop at $N$.
:::

::: proof
$H_i(n) \le (\sum_{m \le n} 1/m)^{d} \le (1 + \log n)^d$, so the tail is
at most $\sum_{n > N} 2^{-n} (1 + \log n)^d$, and the ratio of
consecutive terms is at most $q$.
:::

With $N \approx P$ terms the cost is $O(M(P) \log^2 P)$ for fixed weight.
`acb_mzeta_mzv` switches to this algorithm from 25000 bits (for weight at
most 500), the measured crossover.

## Multiple polylogarithms {#polylog}

Goncharov's iterated integrals [@Gon2001]
$$
G(a_1, \ldots, a_L; z) = \int_{z > t_1 > \cdots > t_L > 0} \frac{dt_1}{t_1 - a_1} \cdots \frac{dt_L}{t_L - a_L}
\qquad (a_L \ne 0),
$$
along the segment $[0, z]$, give
$$
\operatorname{Li}_{s_1, \ldots, s_r}(x_1, \ldots, x_r) = (-1)^r G(0^{s_1-1}, 1/y_1, \ldots, 0^{s_r-1}, 1/y_r; 1),
\qquad y_i = x_1 \cdots x_i,
$$
which follows by expanding each $dt/(t - a) = -a^{-1} \sum_k (t/a)^k dt$
and integrating term by term.

**Power series at 0.** Write $F_l(t) = G(a_l, \ldots, a_L; t) = \sum_n f_l(n) t^n$.
Then $f_l(n) = f_{l+1}(n)/n$ if $a_l = 0$, and $f_l(n) = -S_l(n-1)/(n a_l)$
with $S_l(n) = S_l(n-1)/a_l + f_{l+1}(n)$ otherwise. After scaling by
$z^n$, the values $T_l(n) = z^n S_l(n)$ for the nonzero letters and the
partial sums $A_l(n) = \sum_{m \le n} f_l(m) z^m \to F_l(z)$ for all
$l$ satisfy a lower triangular linear recurrence (with $g_{l+1}$
expanded down to the next nonzero letter). Multiplying the rows by
$n^C$, $C$ one more than the longest run of zeros, the product of the
matrices is computed by binary splitting in complex ball arithmetic
(exactly while the entries are exact and short). This gives
$G(a_l, \ldots, a_L; z)$ for all suffixes at once.

::: theorem Truncation {#thm-mpl-trunc}
Let $R$ be the smallest absolute value of a nonzero letter, $d$ the
number of nonzero letters, and $q = |z|/R < 1$. Truncating the series
after $n = N$ changes each $F_l(z)$ by at most
$$
\log(N+1)^d\, e\, (1 + 1/N)\, \frac{q^{N+1}}{1 - q (1 + 1/N)}
\qquad (q (1 + 1/N) < 1).
$$
:::

::: proof
The coefficients of $-1/(t - a) = a^{-1} \sum_k (t/a)^k$ are bounded in
absolute value by those of $1/(|a| - t)$, and integration preserves
coefficientwise majorization, so $|f_l(n)| \le \tilde f_l(n)$, the
coefficients of the series $\tilde F_l$ with the letters $|a_l|$. For
$0 < \rho < R$, by induction from $\tilde F_{L+1} = 1$:
a letter 0 does not increase $\tilde F(\rho)$ (for a series with
nonnegative coefficients vanishing at 0, $\sum_n c_n \rho^n / n \le \sum_n c_n \rho^n$),
and a nonzero letter $b$ multiplies it by at most
$\int_0^\rho dt/(|b| - t) \le \log(R/(R - \rho))$, since $\tilde F$ is
increasing. Hence $\tilde F_l(\rho) \le \log(R/(R-\rho))^d$ and, by
Cauchy's inequality, $|f_l(n) z^n| \le \log(R/(R-\rho))^d (|z|/\rho)^n$.
With $\rho = R N/(N+1)$, $\log(R/(R - \rho)) = \log(N+1)$ and
$|z|/\rho = q (1 + 1/N)$; sum the geometric tail and use
$(1 + 1/N)^{N+1} \le e (1 + 1/N)$.
:::

**Splitting at $\tfrac12$.** The series at 1 converges geometrically
with ratio $1/R$. If $R \ge 2$ it is used directly; otherwise, if every
nonzero letter $a$ and every nonzero $1 - a$ has absolute value at least
$\tfrac34$, the path is split at $\tfrac12$ (Chen's formula) and
$t \mapsto 1 - t$ is substituted above it:
$$
G(a_1, \ldots, a_L; 1) = \sum_{j=0}^{L} (-1)^j\, G(1 - a_j, \ldots, 1 - a_1; \tfrac12)\, G(a_{j+1}, \ldots, a_L; \tfrac12),
$$
where both families are suffixes of one word (the second of $a$, the
first of the reversed word $1 - a_L, \ldots, 1 - a_1$), so two binary
splitting runs with ratio at most $\tfrac23$ give all the factors.
Arguments outside these regions, such as letters close to $\tfrac12$ or
points $x$ outside the unit disk, are not supported.

**Colored multiple zeta values.** For $x_i = \zeta^{e_i}$ with
$\zeta = e^{2\pi i/N}$, the partial products $y_i$ are formed exactly
from the sums of the exponents. The letters $1/y_i$ are roots of unity
and the letters $1 - 1/y_i$ have absolute value $|1 - \zeta^j| \ge 2 \sin(\pi/N)$
for $j \not\equiv 0$, so the split applies for $N \le 8$ (with ratio at
most $\tfrac12$ for $N \le 6$).

## Use in FLINT {#use}

- `acb_mzeta_mzv`, `acb_mzeta_mzv_star`, `acb_mzeta_mzv_t`,
  `acb_mzeta_mzv_all`: double tails; Hölder convolution in
  `acb_mzeta_mzv` from 25000 bits.
- `acb_mzeta_polylog`, `acb_mzeta_mzv_colored`: Goncharov's iterated
  integrals by binary splitting.
