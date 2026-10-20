---
title: Discrete Fourier transforms over a generic ring
summary: The DFT algorithms of `gr_dft` (naive, radix-2 and split-radix Cooley–Tukey, four-step, mixed-radix with direct prime kernels, Good–Thomas and Bluestein), the conditions each one needs on the root of unity, and their costs.
complexity: $O(n \log n)$ ring operations for every $n$ over rings with roots of unity of all orders (Bluestein for prime factors above 19); otherwise $O(n \sum_i p_i)$ for $n = \prod_i p_i$, using direct kernels for the prime factors.
order: 60
---

The module `gr_dft` computes the discrete Fourier transform

$$
X_k = \sum_{j=0}^{n-1} x_j w^{jk}, \qquad 0 \le k < n,
$$

of a vector $x \in R^n$ over an arbitrary commutative ring $R$, given an
element $w \in R$ with $w^n = 1$. Unlike the transforms of `fft_small` and
`fft`, which are specialised to particular rings and lengths, it works with
any length and only uses the ring operations of `gr`, so the same code
serves complex balls, exact algebraic numbers, $\Z/m\Z$, finite fields and
fixed-point numbers. The classical references are the survey of Duhamel and Vetterli [@DV1990]
and the book of Nussbaumer [@Nus1982]; transforms over general rings are
treated by von zur Gathen and Gerhard [@GG2003, §8.2], Pollard [@Pol1971]
and Bernstein [@Ber2001].

## Roots of unity in a ring {#principal}

Over $\mathbb{C}$ the fast algorithms use identities of the roots of unity
such as $w^{n/2} = -1$ and $1 + v + \cdots + v^{p-1} = 0$ for a root $v$ of
order $p$. In a general ring these identities do not follow from $w$
having order $n$, because the ring can have zero divisors.

::: definition {#def-principal}
An element $w \in R$ is a *principal $n$-th root of unity* if $w^n = 1$
and $w^k - 1$ is not a zero divisor for $1 \le k < n$.
:::

::: lemma {#lem-principal}
$w$ is a principal $n$-th root of unity if and only if $w^n = 1$ and
$w^{n/p} - 1$ is not a zero divisor for every prime $p \mid n$. If $R$ is
an integral domain, this holds if and only if $w$ has order exactly $n$.
:::

::: proof
One direction is clear. Conversely, let $1 \le k < n$ and
$d = \gcd(k, n)$, a proper divisor of $n$, so $d \mid n/p$ for some prime
$p \mid n$. If $a = bc$ in $R$ and $a$ is not a zero divisor, then neither
is $b$, since $bx = 0$ implies $ax = 0$. Now $w^d - 1$ divides
$w^{n/p} - 1$, so it is not a zero divisor. Choose $u \ge 1$ with
$uk \equiv d \pmod n$; then $w^d - 1 = w^{uk} - 1$ is a multiple of
$w^k - 1$, which is therefore not a zero divisor either. In an integral
domain the nonzero elements are exactly the non-zero-divisors.
:::

Over $\Z/m\Z$ the lemma says that $w$ must have order exactly $n$ modulo
every prime factor of $m$. Over $\Z/p\Z$ with $p$ prime, such roots exist
exactly when $n \mid p - 1$.

::: lemma Orthogonality {#lem-orth}
If $w$ is a principal $n$-th root of unity and $k \not\equiv 0 \pmod n$,
then $\sum_{j=0}^{n-1} w^{jk} = 0$. For every prime $p \mid n$, the
element $v = w^{n/p}$ is a principal $p$-th root of unity, so
$1 + v + \cdots + v^{p-1} = 0$; in particular $w^{n/2} = -1$ if $n$ is
even.
:::

::: proof
Let $S = \sum_j w^{jk}$ with $0 < k < n$. Then
$(w^k - 1) S = w^{nk} - 1 = 0$, and $w^k - 1$ is not a zero divisor, so
$S = 0$. For $1 \le k < p$, $v^k - 1 = w^{kn/p} - 1$ with $kn/p < n$, so
$v$ is principal, and the first part applies to $v$ with $k = 1$.
:::

::: theorem Inversion {#thm-inverse}
If $w$ is a principal $n$-th root of unity, then so is $w^{-1} = w^{n-1}$,
and applying the transform with $w^{-1}$ to the transform of $x$ with $w$
gives $n x$.
:::

::: proof
The first claim follows from the definition, since
$w^{-k} - 1 = -w^{-k}(w^k - 1)$. The composition maps $x$ to
$\sum_k w^{-jk} \sum_l x_l w^{lk} = \sum_l x_l \sum_k w^{(l - j)k}$, and
by [](#lem-orth) the inner sum is $n$ for $l = j$ and 0 otherwise.
:::

The inverse transform of `gr_dft` is therefore the transform with
$w^{-1}$ followed by a division by $n$, and it fails with `GR_DOMAIN` when
the ring cannot divide by $n$. The same computation shows that the
transform maps cyclic convolution to pointwise multiplication.

The algorithms below need less than principality to compute the forward
transform as defined. Their derivations use only the following
identities:

| Algorithm | Identities used |
|---|---|
| naive, prime factor | $w^n = 1$ |
| radix-2, split-radix, four-step | $w^{n/2} = -1$ |
| mixed-radix, direct kernel for the prime $p$ | $w^n = 1$ and $1 + v + \cdots + v^{p-1} = 0$ for $v = w^{n/p}$ |
| Bluestein | $w^n = 1$ and $n$ odd; division by $N$ and a principal root of order $N = 2^k \ge 2n - 1$ |

The prime factor algorithm passes the roots $w^{n_2}$ and $w^{n_1}$ to its
sub-plans, which need their own identities. For $n = 2^d$, the condition
$w^{n/2} = -1$ already implies orthogonality: for $k = 2^s k'$ with $k'$
odd and $s < d$, the terms $j$ and $j + n/2^{s+1}$ of
$\sum_j w^{jk}$ cancel, because $w^{k n / 2^{s+1}} = (w^{n/2})^{k'} = -1$.
So for powers of two $w^{n/2} = -1$ is all that is needed.

`gr_dft_precomp_init_root` checks $w^{n/2} = -1$ for even $n$ and
$w^{n/p} \ne 1$ for odd primes $p \mid n$, as far as the ring can decide
them with `gr_is_neg_one` and `gr_is_one`. The second check is necessary
but not sufficient: over a ring with zero divisors, $w^{n/p} \ne 1$ does
not imply $1 + v + \cdots + v^{p-1} = 0$, and the mixed-radix kernel then
returns a wrong result without an error.
[](gr_dft/plans#example-zero-divisor) gives an example modulo 91.

## The naive algorithm {#naive}

`GR_DFT_ALG_NAIVE` evaluates the definition, accumulating
$x_j w^{jk \bmod n}$ with the exponent updated incrementally. It performs
$n(n-1)$ calls to the root multiplication (those with exponent 0 are
copies), is correct for any $w$ with $w^n = 1$, and serves as the
reference in the tests. For powers of two it also supports the
bit-reversed order of `GR_DFT_SCRAMBLED`.

## Radix-2 Cooley–Tukey {#radix2}

For $n = 2^d$, splitting the input into halves (decimation in frequency,
[@GS1966]) gives

$$
X_{2k} = \sum_{j < n/2} (x_j + x_{j + n/2}) (w^2)^{jk},
\qquad
X_{2k+1} = \sum_{j < n/2} (x_j - x_{j + n/2}) w^j (w^2)^{jk},
$$

because $w^{(j + n/2) \cdot 2k} = w^{2jk}$ and
$w^{(j + n/2)(2k + 1)} = -w^{j(2k+1)}$ by $w^n = 1$ and $w^{n/2} = -1$.
Since $(w^2)^{n/4} = w^{n/2} = -1$, the two halves are transforms of the
same kind with root $w^2$, and induction gives the algorithm of Cooley and
Tukey [@CT1965].

::: algorithm Radix-2 decimation in frequency {#alg-dif}
**Input:** $x_0, \ldots, x_{n-1}$ with $n = 2^d$, and $w$ with
$w^{n/2} = -1$.

**Output:** $X_{\mathrm{rev}(k)}$ in position $k$, where $\mathrm{rev}$
reverses the $d$ bits of $k$.

1. For $s = d, d - 1, \ldots, 1$, with $m = 2^s$ and $r = n/m$: for each
   block $[b, b + m)$ with $b \equiv 0 \pmod m$ and each $0 \le j < m/2$,
   replace $(x_{b+j}, x_{b+j+m/2})$ by
   $(x_{b+j} + x_{b+j+m/2},\ w^{jr}(x_{b+j} - x_{b+j+m/2}))$.
2. Return $x$.
:::

The forward transform of `gr_dft` uses [](#alg-dif) followed by a
bit-reversal permutation, which `GR_DFT_SCRAMBLED` skips. The inverse
transform runs the transposed algorithm (decimation in time,
$(a, b) \mapsto (a + w^{-jr} b, a - w^{-jr} b)$ for $s = 1, \ldots, d$)
on bit-reversed input, so a scrambled forward transform followed by a
scrambled inverse transform gives $n$ times the identity without any
permutation. There are $\frac{n}{2} d$ butterflies, each with one
multiplication by a root; those with $j = 0$ need none, and in complex
mode those with $jr \in \{n/8, n/4, 3n/8\}$ are cheap
([](#complex)). The integer FFTs of `fft`
([](algorithms/fft/schonhage-strassen#ring)) and `fft_small`
([](algorithms/fft/small-prime-fft)) are built on the same recursion,
specialised to their rings, and truncate it so that only the needed
values are computed ([](algorithms/fft/truncated-fft)).

## Split-radix {#split}

The split-radix algorithm of Yavne [@Yav1968] and Duhamel and Hollmann
[@DH1984] splits a transform of length $m$ with root $v$ into one
transform $U$ of length $m/2$ (the even entries) and two transforms $Z$,
$Z'$ of length $m/4$ (the entries $\equiv 1$ and $\equiv 3 \pmod 4$).
With $t = v^k Z_k$, $t' = v^{3k} Z'_k$ and $\iota = v^{m/4}$,

$$
X_k = U_k + (t + t'), \quad
X_{k + m/2} = U_k - (t + t'), \quad
X_{k + m/4} = U_{k + m/4} + \iota (t - t'), \quad
X_{k + 3m/4} = U_{k + m/4} - \iota (t - t'),
$$

for $0 \le k < m/4$, which again uses only $v^{m/2} = -1$. Each step
costs two multiplications by general roots per $k$ plus one by $\iota$,
which plays the role of $-i$ and is free in complex mode. The total is
about $\frac{n}{3} \log_2 n$ general multiplications instead of
$\frac{n}{2} \log_2 n$. In complex mode, with three real multiplications
per general root and two per eighth root, it uses $n \log_2 n - 3n + 4$
real multiplications, which is minimal for $n \le 16$ [@Win1978] and
close to the lower bound of Heideman and Burrus [@HB1986] for larger
$n$. Later variants reduce the number of additions further [@JF2007];
they are not used. `gr_dft` computes split-radix recursively, out of
place, in natural order.

## The four-step algorithm {#four-step}

For $n = n_1 n_2$ with $n_1 = 2^{\lfloor d/2 \rfloor}$, Bailey's four-step
algorithm [@Bai1990] writes $j = j_1 + n_1 j_2$, $k = k_2 + n_2 k_1$ and

$$
X_{k_2 + n_2 k_1} = \sum_{j_1} (w^{n_2})^{j_1 k_1}\, w^{j_1 k_2}
\sum_{j_2} x_{j_1 + n_1 j_2} (w^{n_1})^{j_2 k_2},
$$

which gives $n_1$ transforms of length $n_2$, $n$ twiddle
multiplications by $w^{j_1 k_2}$, $n_2$ transforms of length $n_1$ (all
radix-2 sub-plans) and a transposition, which `GR_DFT_SCRAMBLED` skips.
The method and its variants are described in
[](algorithms/fft/matrix-fourier). In `gr_dft` it is available as
`GR_DFT_ALG_BAILEY` but not chosen automatically.

## Mixed-radix Cooley–Tukey {#mixed}

For a transform of length $\ell = pm$ with $p$ prime and root $u$, split
the input into its $p$ residue classes modulo $p$ and let $Y_b$ be the
transform of length $m$, with root $u^p$, of the class $b$. Then

$$
X_{q + m r} = \sum_{b=0}^{p-1} v^{br} \bigl(u^{bq} Y_b(q)\bigr),
\qquad v = u^m,\ 0 \le q < m,\ 0 \le r < p,
$$

so after twiddling by $u^{bq}$, each $q$ requires a transform of prime
length $p$ with root $v$. `GR_DFT_ALG_MIXED` applies this recursively,
one prime factor of $n$ per level, with all exponents taken in the root
table of the top-level plan. The transforms of prime length use the
following kernel.

::: lemma Direct prime kernel {#lem-kernel}
Let $p$ be prime, $v^p = 1$ and $1 + v + \cdots + v^{p-1} = 0$. For
$0 < r < p$ let $b^*$ be the index with $b^* r \equiv p - 1 \pmod p$ and
$e_b = br \bmod p$. Then

$$
\sum_{b=0}^{p-1} v^{br} z_b = (z_0 - z_{b^*}) +
\sum_{b \ne 0,\, b^*} v^{e_b} (z_b - z_{b^*}).
$$
:::

::: proof
As $b$ runs over $0, \ldots, p-1$, so does $e_b$, hence
$\sum_b v^{e_b} = 0$, and subtracting $z_{b^*} \sum_b v^{e_b}$ from the
left side changes nothing. The term $b = 0$ has $e_0 = 0$, and the term
$b = b^*$ vanishes.
:::

The output for $r = 0$ is $\sum_b z_b$, so the kernel takes
$(p-1)(p-2)$ multiplications instead of $(p-1)^2$. For $p = 3$, for
example, $X_1 = (x_0 - x_2) + v(x_1 - x_2)$ and
$X_2 = (x_0 - x_1) + v(x_2 - x_1)$. With the $(p - 1)$ twiddles per
kernel, a transform of length $p^e$ costs at most $e\, p^{e-1} (p-1)^2$
multiplications, which is $O(n \log n)$ only for bounded $p$. For
$n = p^e$ with $p > 19$, the transforms of length $p$ use Bluestein's
algorithm instead when it is available.

## The prime factor algorithm {#pfa}

If $n = n_1 n_2$ with $\gcd(n_1, n_2) = 1$, the Chinese remainder theorem
removes the twiddle factors altogether. This is the algorithm of Good
[@Goo1958] and Thomas [@Tho1963].

::: theorem Good–Thomas {#thm-pfa}
Let $n = n_1 n_2$ with $\gcd(n_1, n_2) = 1$ and $w^n = 1$. Let $a, b$ be
the integers in $[0, n)$ with $a \equiv 1 \pmod{n_1}$,
$a \equiv 0 \pmod{n_2}$, $b \equiv 0 \pmod{n_1}$ and
$b \equiv 1 \pmod{n_2}$. Then for all $k_1 < n_1$ and $k_2 < n_2$,

$$
X_{(a k_1 + b k_2) \bmod n} = \sum_{j_1 < n_1} \sum_{j_2 < n_2}
x_{(n_2 j_1 + n_1 j_2) \bmod n}\, (w^{n_2})^{j_1 k_1} (w^{n_1})^{j_2 k_2},
$$

and both index maps are bijections between $[0, n_1) \times [0, n_2)$ and
$[0, n)$.
:::

::: proof
The maps are bijections by the Chinese remainder theorem, since
$n_2 j_1 + n_1 j_2 \equiv n_2 j_1 \pmod{n_1}$ and $n_2$ is invertible
modulo $n_1$, and symmetrically. Let $k = a k_1 + b k_2$, so that
$k \equiv k_1 \pmod{n_1}$ and $k \equiv k_2 \pmod{n_2}$. Then
$w^{(n_2 j_1 + n_1 j_2) k} = (w^{n_2})^{j_1 k} (w^{n_1})^{j_2 k}$, and
since $(w^{n_2})^{n_1} = w^n = 1$ the exponent $j_1 k$ can be reduced
modulo $n_1$ to $j_1 k_1$, and likewise $j_2 k$ to $j_2 k_2$.
:::

The right side is a two-dimensional transform: $n_1$ transforms of
length $n_2$ with root $w^{n_1}$ followed by $n_2$ transforms of length
$n_1$ with root $w^{n_2}$, with a gather before and a scatter after but no
multiplications in between. If $w$ is principal, so are $w^{n_2}$ and
$w^{n_1}$, since $(w^{n_2})^{n_1/p} - 1 = w^{n/p} - 1$. `GR_DFT_ALG_PFA`
takes $n_1$ to be the full power of the smallest prime factor of $n$ and
plans $n_2$ automatically, so a general length is reduced to its prime
power factors. With the default root, its sub-plans build their own
tables and the plan itself stores none.

## Bluestein's algorithm {#bluestein}

Bluestein's algorithm [@Blu1970], also known as the chirp z-transform
[@RSR1969], reduces a transform of any odd length $n$ to a cyclic
convolution, which can then be computed with transforms of a convenient
length.

::: theorem {#thm-bluestein}
Let $n$ be odd, $w^n = 1$ and $t = (n + 1)/2$. With $c_j = w^{t j^2}$ and
$h_j = w^{-t j^2}$,

$$
X_k = c_k \sum_{j=0}^{n-1} (c_j x_j)\, h_{(k - j) \bmod n},
$$

and $h$ is even and $n$-periodic: $h_{-j} = h_j = h_{j+n}$.
:::

::: proof
Since $2t = n + 1 \equiv 1 \pmod n$,
$jk \equiv t \cdot 2jk = t(j^2 + k^2 - (k - j)^2) \pmod n$, and $w^n = 1$
gives $w^{jk} = c_j c_k h_{k-j}$. Periodicity holds because
$(j + n)^2 \equiv j^2 \pmod n$.
:::

The cyclic convolution of length $n$ is embedded in one of length
$N = 2^{\lceil \log_2(2n - 1) \rceil}$: the vector $(c_j x_j)$ is padded
with zeros, and the kernel is $h'_j = h_j$ for $0 \le j < n$,
$h'_{N - j} = h_j$ for $1 \le j < n$ and 0 otherwise. Because
$N \ge 2n - 1$ the two parts do not overlap, and for $0 \le k < n$ the
cyclic convolution of length $N$ agrees with that of length $n$ by the
evenness of $h$.

::: algorithm Bluestein {#alg-bluestein}
**Precomputation:** a radix-2 plan of length $N$ with a principal root of
order $N$, in scrambled order, and $\hat h = \frac{1}{N} \mathrm{DFT}_N(h')$.

1. $z_j \gets c_j x_j$ for $j < n$, and $z_j \gets 0$ for $n \le j < N$.
2. $Z \gets \mathrm{DFT}_N(z)$, in scrambled order.
3. $Z_k \gets Z_k \hat h_k$ for $k < N$.
4. $y \gets$ the unscaled inverse transform of $Z$, from scrambled order.
5. Return $X_k = c_k y_k$ for $k < n$.
:::

The chirps $c_j$ and $h_j$ are entries of the root table of $w$, and
$\hat h$ is in the same scrambled order as $Z$, so no permutation is
needed. The inverse transform of length $n$ is computed as the forward
transform of $x_{(n - j) \bmod n}$, since
$\sum_j x_j w^{-jk} = \sum_j x_{(n-j) \bmod n} w^{jk}$. The cost is $2n$
chirp multiplications, $N$ pointwise products and two transforms of
length $N < 4n$; the precomputation costs one more transform. The
convolution needs a root of order $N$, unrelated to $w$, which `gr_dft`
obtains from `gr_dft_default_root`. Over rings without default roots,
such as $\Z/m\Z$ (even when $m$ has roots of order $N$), the algorithm is
therefore not available, and mixed-radix plans fall back to the direct
kernel. Over the fixed-point contexts the kernel is scaled by $1/(2N)$
instead of $1/N$ to keep the intermediate values below 1, and the outputs
are doubled at the end.

::: remark
Rader's algorithm [@Rad1968] reduces a transform of prime length $p$ to a
cyclic convolution of length $p - 1$ by means of a primitive root modulo
$p$, and Winograd's algorithm [@Win1978] combines it with short
convolution algorithms. `gr_dft` uses neither: Bluestein's convolution
length is always a power of two, which suits the radix-2 transforms.
:::

## Products of cyclic groups and Dirichlet characters {#prod}

The DFT on $\Z/c_0\Z \times \cdots \times \Z/c_{r-1}\Z$,
$X_k = \sum_j x_j \prod_a w_a^{j_a k_a}$ with $w_a$ a principal $c_a$-th
root of unity, factors over the axes. It is computed by transforming all
lines along each axis in turn with a cyclic plan of length $c_a$, with no
twiddle factors, in $\sum_a (n / c_a) C(c_a)$ operations for
$n = \prod c_a$, where $C(m)$ is the cost of a transform of length $m$.
Plans are shared between axes of equal length, and the inverse divides
by $n$ once. By [](#thm-pfa), a cyclic transform of coprime length
$n_1 n_2$ is a product transform up to permutations of the input and the
output.

The group of Dirichlet characters modulo $q$ is isomorphic to
$(\Z/q\Z)^\times$, which the Conrey decomposition writes as a product of
cyclic groups with generators of orders $\varphi_k$. If $x$ has
logarithms $b_k$ and the character $\chi_a$ is indexed by the logarithms
$a_k$ of $a$, then $\chi_a(x) = \prod_k e^{2\pi i a_k b_k / \varphi_k}$,
and the transform $\hat v_a = \sum_x \overline{\chi_a(x)} v_x$ is the
product transform over the $\varphi_k$ in the lexicographic order of the
logarithms. `gr_dft_dirichlet` adds the permutation between this order
and the order of the residues $x$.

## Complex mode {#complex}

When the ring is $C = R[i]$ with elements stored as pairs of elements of
$R$ and $w = e^{-2\pi i/n}$, multiplications by roots can be specialised
from the exponent alone, without inspecting values (which over balls
could not prove that a computed value is exactly $-i$):

- $w^{jn/4} \in \{1, -i, -1, i\}$: swaps and negations of the parts;
- $w^{jn/8}$ for odd $j$ equals $c(1 \pm i)$ with $c = \pm\sqrt{2}/2$, so
  that, for example, $c(1 + i)(x + iy) = c(x - y) + i\, c(x + y)$ takes
  two real multiplications;
- for a general root $c + di$, the product with $x + iy$ is
  $k_1 - y(d + c) + i\,(k_1 + x(d - c))$ with $k_1 = c(x + y)$: three real
  multiplications and three additions with the precomputed $c$, $d - c$
  and $d + c$, instead of four multiplications and two additions. This
  form is attributed to Gauss [@Knu1997, §4.6.4].

The last trick pays when multiplications in $R$ are expensive compared
to additions, so `gr_dft` uses it over exact rings and when the real
precision is at least 1024 bits; below that, a native complex
multiplication is faster. In fixed point the entries $d \pm c$ can reach
$\sqrt 2$ and do not fit the representation, so complex products use the
schoolbook formula below 22 limbs and a carry-extended three-product
variant from 22 limbs.

## Error handling {#errors}

The algorithms only call `gr` operations and combine their status codes
with `|`. Over exact rings the forward transforms cannot fail, but the
division by $n$ of the inverse can, with `GR_DOMAIN`, as can the division
by $N$ in the Bluestein precomputation. Over inexact rings the root
checks may be inconclusive (`T_UNKNOWN`), in which case the root is
accepted. Over the fixed-point contexts, an overflow during a transform
returns `GR_UNABLE`; [](algorithms/fft/fixed-point-dft) shows how to
scale the input so that this cannot happen and bounds the rounding
errors.

## Summary {#summary}

| Algorithm | Lengths | Multiplications by roots | Conditions |
|---|---|---|---|
| naive | any | $(n - 1)^2$ nontrivial | $w^n = 1$ |
| radix-2 | $2^d$ | $\frac{n}{2} d$ | $w^{n/2} = -1$ |
| split-radix | $2^d$, $d \ge 2$ | about $\frac{n}{3} d$ general ones | $w^{n/2} = -1$ |
| four-step | $2^d$, $d \ge 2$ | radix-2 counts plus $n$ | $w^{n/2} = -1$ |
| mixed-radix | any $n \ge 2$ | $\sum_{\text{levels}} n (p_i - 1)^2 / p_i$ | $\sum_{i<p} v^i = 0$ for $v = w^{n/p}$ and each prime $p \mid n$ |
| prime factor | $\ge 2$ distinct primes | those of the sub-plans | those of the sub-plans |
| Bluestein | odd $n \ge 3$ | $O(N \log N)$ with $2n - 1 \le N < 4n$ | default root of order $N$; division by $N$ |
