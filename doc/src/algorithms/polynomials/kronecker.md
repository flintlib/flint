---
title: Kronecker substitution
summary: Multiplies polynomials over $\Z/n\Z$ by packing their coefficients into large integers, multiplying the integers and unpacking; Harvey's KS2 and KS4 variants evaluate at ±2^b and ±2^(-b) to split the integer product into two or four smaller ones.
complexity: One integer product of about $(\ell_1 + \ell_2)(2\beta + \lg \ell_2)$ bits (KS), or two (KS2) or four (KS4) products of a half or a quarter of that size, plus $O(\ell_1 + \ell_2)$ word operations.
order: 10
---

Kronecker substitution reduces polynomial multiplication to integer
multiplication. The coefficients of each factor are written into the
bits of one large integer, with enough space between them that the
coefficients of the product do not overlap, and one integer product
yields all coefficients of the polynomial product. The idea goes back to
Kronecker [@Kro1882]; Schönhage used it to multiply polynomials with
integer and complex coefficients [@Sch1982], and it is discussed in
[@GG2003, §8.4] and [@Ber2008]. Its strength is that it inherits
all the work that has gone into integer multiplication: FLINT's and
GMP's basecase and Toom–Cook code, and the FFTs.

Harvey [@Har2009a] observed that evaluating at several points, $\pm 2^b$
and $\pm 2^{-b}$, replaces the single integer product by two or four
products of a half or a quarter of the size. Since integer
multiplication costs more than linear time, this is faster. FLINT's
implementations of these variants, KS2 and KS4, come from Harvey's
`zn_poly` library.

Throughout, $f$ and $g$ have lengths $\ell_1 \ge \ell_2$ and coefficients
in $[0, n)$, $h = fg$ has length $\ell_1 + \ell_2 - 1$, and $\beta$ is
the bit length of $n$ (`NMOD_BITS`). We write $h_k$ for the coefficients
of the product of $f$ and $g$ as polynomials with integer coefficients.
Each is a sum of at most $\ell_2$ products, so

$$
0 \le h_k \le \ell_2 (n - 1)^2 < 2^{2\beta + \lceil \log_2 \ell_2 \rceil}.
$$

All three methods compute the integers $h_k$ exactly and reduce them
modulo $n$ at the end, with one-, two- or three-word reductions
(`NMOD_RED`, `NMOD2_RED2`, `NMOD_RED3`).

## Plain Kronecker substitution {#ks}

::: lemma {#lem-ks}
Let $b \ge 1$ and let $h(x) = \sum h_k x^k$ have integer coefficients
$0 \le h_k < 2^b$. Then $h_k$ is the integer formed by bits
$bk, \ldots, bk + b - 1$ of $h(2^b)$.
:::

::: proof
The binary expansions of the terms $h_k 2^{bk}$ occupy disjoint bit
ranges, so their sum is formed without carries.
:::

::: algorithm KS {#alg-ks}
**Input:** $f$, $g$ as above.

**Output:** $fg \bmod n$.

1. Choose $b$ with $\ell_2 (n - 1)^2 < 2^b$.
2. Form the integers $F = f(2^b)$ and $G = g(2^b)$ by writing the
   coefficients into $b$-bit fields (bit packing).
3. Compute the integer product $H = FG$.
4. Read the $h_k$ from the $b$-bit fields of $H$ (bit unpacking) and
   reduce each modulo $n$.
:::

Correctness is [](#lem-ks) applied to the integer polynomial $fg$, since
$FG = (fg)(2^b)$. The integers $F$ and $G$ have about $\ell_1 b$ and
$\ell_2 b$ bits, so the cost is one integer product of
$(\ell_1 + \ell_2) b$ bits plus linear work. The field size
$b \approx 2\beta + \log_2 \ell_2$ grows only logarithmically with the
length, but for small moduli a large part of each field holds the bits
of the sums rather than of the inputs, which is the waste that KS2 and
KS4 reduce.

`_nmod_poly_mulmid_KS`, which implements the full, low and middle
products by substitution, takes $b = 2\beta + \lfloor \log_2 \ell_2
\rfloor + 1$, and if this fits in a word replaces it by the exact bit
length of $\ell_2 (n-1)^2$. It packs with `_nmod_poly_bit_pack`,
multiplies with `flint_mpn_mul` (or `flint_mpn_sqr` for squares) and
extracts only the requested range of coefficients with
`_nmod_poly_bit_unpack`, which reduces each field of one to three words.
Before packing it drops the low coefficients of the operands that cannot
contribute to the requested range, but otherwise it forms the whole
integer product.

## Evaluation at $2^b$ and $-2^b$ (KS2) {#ks2}

Write $h(x) = h_e(x^2) + x h_o(x^2)$ for the even and odd parts, and
similarly for $f$ and $g$. Let $B = 2^b$.

::: lemma {#lem-ks2}
$h(B) + h(-B) = 2 h_e(B^2)$ and $h(B) - h(-B) = 2 B h_o(B^2)$. If the
coefficients of $h$ are less than $B^2$, the coefficients of $h_e$ and
$h_o$ are the $2b$-bit fields of $h_e(B^2)$ and $h_o(B^2)$.
:::

::: proof
The identities follow from $(-B)^k = (-1)^k B^k$. The second statement
is [](#lem-ks) with $2b$ in place of $b$.
:::

With $b$ about half the field size of KS, KS2 computes two integer
products $f(B) g(B)$ and $f(-B) g(-B)$, each of about
$(\ell_1 + \ell_2) b$ bits, half the size of the single product of KS.
The evaluations are cheap: $f_e(B^2)$ and $B f_o(B^2)$ are packings of
the even and odd coefficients into $2b$-bit fields, the second shifted
by $b$ bits, and $f(\pm B)$ is their sum or difference. Since $f(-B)$
may be negative, FLINT stores absolute values and tracks the sign of the
product.

::: algorithm KS2 {#alg-ks2}
**Input:** $f$, $g$ as above, $\ell_2 \ge 2$.

**Output:** $fg \bmod n$.

1. $c \gets 2\beta + \lceil \log_2 \ell_2 \rceil$, $b \gets \lceil c / 2 \rceil$, $B \gets 2^b$.
2. Pack $f_e(B^2)$, $B f_o(B^2)$, $g_e(B^2)$ and $B g_o(B^2)$
   (`_nmod_poly_KS2_pack`); form $f(\pm B)$ and $g(\pm B)$ with an
   addition and a signed subtraction each.
3. $P \gets f(B) g(B)$, $M \gets f(-B) g(-B)$.
4. Unpack the $2b$-bit fields of $(P + M)/2$ to obtain $h_e$, and of
   $(P - M)/(2B)$ to obtain $h_o$ (`_nmod_poly_KS2_unpack`); reduce modulo
   $n$ into the even and odd positions of the output
   (`_nmod_poly_KS2_reduce`).
:::

The bound $h_k < 2^c \le B^2$ makes [](#lem-ks2) applicable. For a square
only two integer squarings are needed.

## Adding the reciprocal points $\pm 2^{-b}$ (KS4) {#ks4}

The reversal $\operatorname{rev}(f)(x) = x^{\ell_1 - 1} f(1/x)$
satisfies $\operatorname{rev}(fg) = \operatorname{rev}(f)
\operatorname{rev}(g)$, so evaluating reversed polynomials at $B$
amounts to evaluating at $1/B$ in integer arithmetic. KS4 runs the
KS2 evaluation twice, on $f, g$ and on their reversals, with $b$ about a
quarter of the field size of KS. Now a coefficient of $h_e$ may have up
to $4b$ bits while the fields of $h_e(B^2)$ are only $2b$ bits wide, so
neighbouring coefficients overlap. The low halves of the coefficients
can be read from the normal evaluation, starting at the low end, and
the high halves from the reversed one, starting at the high end.

::: theorem Recovery from two overlapping evaluations {#thm-recover}
Let $D \ge 4$ and let $c_0, \ldots, c_{m-1}$ be integers with
$c_i = \alpha_i + D \beta_i$, $0 \le \alpha_i < D$ and
$0 \le \beta_i \le D - 2$. Let

$$
X = \sum_{i=0}^{m-1} c_i D^i, \qquad Y = \sum_{i=0}^{m-1} c_i D^{m-1-i},
$$

with base-$D$ digits $x_p$ and $y_p$ for $0 \le p \le m$. Set
$\alpha_{-1} = \beta_{-1} = \kappa_0 = 0$. Then for $i = 0, 1, \ldots, m-1$,

$$
\alpha_i = (x_i - \beta_{i-1} - \kappa_i) \bmod D, \qquad
\kappa_{i+1} = [\, x_i < \beta_{i-1} + \kappa_i \,],
$$

$$
\beta_i = \bigl((y_{m-i} - \alpha_{i-1}) \bmod D\bigr) - [\, y_{m-1-i} < \alpha_i \,],
$$

where $[\cdot]$ is 1 if the condition holds and 0 otherwise. In
particular the $c_i$ are determined by $X$ and $Y$ and are computed in
increasing order of $i$ with $O(m)$ digit operations.
:::

::: proof
*Digits of $X$.* Coefficient $c_i$ contributes $\alpha_i$ to position
$i$ and $\beta_i$ to position $i + 1$. Let $\kappa_i$ be the carry into
position $i$. Then $\alpha_i + \beta_{i-1} + \kappa_i = x_i + D
\kappa_{i+1}$; the left side is at most $2D - 2$, so carries are 0 or 1.
Reducing modulo $D$ gives $\alpha_i$, and $\kappa_{i+1} = 1$ exactly when
the subtraction $x_i - \beta_{i-1} - \kappa_i$ borrows.

*Digits of $Y$.* Here $c_i$ contributes $\alpha_i$ to position
$p = m - 1 - i$ and $\beta_i$ to position $p + 1$. Let $\lambda_p$ be
the carry into position $p$. Position $p$ receives
$\alpha_i + \beta_{i+1} + \lambda_p$ (with $\beta_m = 0$). Since
$\beta_{i+1} + \lambda_p \le D - 1$, this sum reaches $D$, which makes
$\lambda_{p+1} = 1$, exactly when its residue $y_p$ is less than
$\alpha_i$. Position $p + 1 = m - i$ receives
$\alpha_{i-1} + \beta_i + \lambda_{p+1}$, hence
$(y_{m-i} - \alpha_{i-1}) \bmod D = \beta_i + \lambda_{p+1}$, with no
wrap-around because $\beta_i + \lambda_{p+1} \le D - 1$. Subtracting
$\lambda_{p+1} = [y_{m-1-i} < \alpha_i]$ gives $\beta_i$. Each step
needs only $\beta_{i-1}$ and $\alpha_i$, which are known by then.
:::

In KS4, $D = 2^{2b}$, the $c_i$ are the coefficients of $h_e$ (and then
of $h_o$), $X$ is $h_e(B^2)$ from the normal evaluation and $Y$ is the
corresponding value from the reversed one. The hypothesis on $\beta_i$
holds with FLINT's choice of $b$:

::: proposition {#prop-ks4-bound}
Let $c = 2\beta + k$ with $k = \lceil \log_2 \ell_2 \rceil$,
$b = \lceil c/4 \rceil$ and $D = 2^{2b}$, as in `_nmod_poly_mul_KS4`.
Then every $h_j$ satisfies $\lfloor h_j / D \rfloor \le D - 2$.
:::

::: proof
Since $n - 1 \le 2^\beta - 1$ and $\ell_2 \le 2^k$,

$$
h_j \le 2^k (2^{2\beta} - 2^{\beta + 1} + 1) = 2^c - 2^k (2^{\beta+1} - 1).
$$

If $4b = c$, then $D^2 = 2^c$ and $D = 2^{\beta + k/2} < 2^k
(2^{\beta+1} - 1)$, so $h_j < D^2 - D$, which is the claim. If $4b > c$,
then $h_j < 2^c \le D^2/2$, so $\lfloor h_j / D \rfloor < D/2 \le D - 2$.
:::

::: algorithm KS4 {#alg-ks4}
**Input:** $f$, $g$ as above, $\ell_2 \ge 2$.

**Output:** $fg \bmod n$.

1. $c \gets 2\beta + \lceil \log_2 \ell_2 \rceil$, $b \gets \lceil c / 4 \rceil$, $B \gets 2^b$.
2. As in [](#alg-ks2), compute $f(\pm B) g(\pm B)$ and from them
   $2 h_e(B^2)$ and $2 B h_o(B^2)$.
3. Do the same with $\operatorname{rev}(f)$ and $\operatorname{rev}(g)$.
   This gives the even and odd parts of $\operatorname{rev}(h)$ at $B^2$,
   which are the reversals of $h_e$ and $h_o$, exchanged if the length of
   $h$ is even (FLINT compensates with shifts by $b$ bits).
4. Split the four integers into $2b$-bit digits and apply [](#thm-recover)
   to obtain the even and the odd coefficients, reducing each modulo $n$
   (`_nmod_poly_KS2_recover_reduce`).
:::

## Cost {#cost}

Let $I(N)$ be the cost of multiplying two $N$-bit integers and let
$N = (\ell_1 + \ell_2) c$ be the size of the product formed by KS. Up to
linear terms,

$$
\text{KS: } I(N), \qquad \text{KS2: } 2\, I(N/2), \qquad \text{KS4: } 4\, I(N/4).
$$

Since $I$ grows faster than linearly, KS4 is the fastest for large $N$.
The gain is large where integer multiplication is quadratic or
Toom–Cook, and shrinks to a modest constant factor in the FFT range
[@Har2009a]. The linear overhead (more packing, signed additions and
the recovery) is larger for KS2 and KS4, so plain KS is preferred for
short products; the crossovers used by `_nmod_poly_mul` are given in
[](algorithms/polynomials/multiplication#cutoffs). For truncated and
middle products, substitution offers no saving beyond shortening the
operands, since the whole integer product is formed.

## Use in FLINT {#use}

- `_nmod_poly_mul_KS`, `_nmod_poly_mullow_KS` and `_nmod_poly_mulmid_KS`
  implement [](#alg-ks); `_nmod_poly_mul_KS2` and `_nmod_poly_mul_KS4`
  implement [](#alg-ks2) and [](#alg-ks4). `_nmod_poly_mul`,
  `_nmod_poly_mullow` and `_nmod_poly_mulmid` choose among them by length
  and modulus size, and above tuned lengths use the small-prime FFT
  instead ([](algorithms/polynomials/multiplication)). `_nmod_poly_mulhigh`
  goes through `_nmod_poly_mulmid` for all but tiny inputs.
- The integer products go through `flint_mpn_mul`, which uses FLINT's
  basecase and Toom-2 code, GMP's `mpn_mul` and, above a tuned size, the
  small-prime FFT of `fft_small` or, when that is not available, the
  Schönhage–Strassen FFT ([](algorithms/fft/schonhage-strassen)).
- The low-level helpers `_nmod_poly_bit_pack`, `_nmod_poly_bit_unpack`,
  `_nmod_poly_KS2_pack`, `_nmod_poly_KS2_unpack`, `_nmod_poly_KS2_reduce`
  and `_nmod_poly_KS2_recover_reduce` are public.
- In other modules, `_fmpz_poly_mul_KS` and its truncated and middle forms
  substitute integer polynomials with signed coefficients;
  `fq_poly_mul_KS` packs elements of $\mathbb{F}_q$ into integers and
  multiplies the resulting integer polynomials ([](fq_poly/mul));
  `fq_nmod_poly_mul_univariate` substitutes $x = y^{2d-1}$ to turn a
  product in $\mathbb{F}_{p^d}[x]$ into one product in $\mathbb{F}_p[y]$;
  and `nmod_poly_mat_mul_KS` packs polynomial matrices into integer
  matrices.

::: functions _nmod_poly_mul_KS _nmod_poly_mul_KS2 _nmod_poly_mul_KS4 _nmod_poly_mulmid_KS :::
