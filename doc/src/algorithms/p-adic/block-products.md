---
title: Error terms of dot products, matrix products and polynomial products
summary: How padic_radix computes the precision of a sum of products term by term, how it computes the values of whole matrix and polynomial products with integer, modular or Kronecker-substitution kernels, and how it computes the per-entry precisions by min-plus products and min-plus convolutions, exactly and in quasi-linear time in the regular cases.
complexity: Precision bookkeeping $O(mnk)$ for an $m \times k$ by $k \times n$ product (often $O(mn)$), and from $O(n)$ to $O(n^2)$ for a product of polynomials of length $n$ (quasi-linear for regular precision patterns); values by $\Z/B^S\Z$ or integer matrix and polynomial multiplication.
order: 12
---

Sums of products dominate $p$-adic linear algebra and polynomial
arithmetic: an entry of a matrix product or a coefficient of a polynomial
product is a dot product $\sum_i a_i b_i$. In the ball model of
`padic_radix` ([](algorithms/p-adic/precision#balls)) such a sum has a
value and an error term. This page describes how both are computed for a
single dot product and for whole matrix and polynomial products, where
the error terms of all outputs form a min-plus product or a min-plus
convolution. Throughout, $p$ is a word-size prime, $B = p^e$ is the limb
radix of `padic_radix` (the largest power of $p$ that fits in a word),
and an element is $p^v u + O(p^N)$.

## The error term of a sum of products {#dot}

::: proposition Precision of a dot product {#prop-dot}
Let $a_i = x_i + O(p^{A_i})$ and $b_i = y_i + O(p^{B_i})$ for
$i = 1, \ldots, n$. Define $v_i' = \operatorname{ord}_p(x_i)$ if $x_i$
is nonzero, $v_i' = A_i$ if $a_i$ is an inexact zero $0 + O(p^{A_i})$,
and $w_i'$ likewise for $b_i$. Then
$$
\sum_i a_i b_i \in \sum_i x_i y_i + O(p^{N_t}), \qquad
N_t = \min_i \min(A_i + w_i', \; B_i + v_i'),
$$
where terms with an exact factor contribute only the error term of the
other factor, and terms with an exact zero factor nothing.
:::

::: proof
By [](algorithms/p-adic/precision#prop-ring) a product of nonzero
balls has error term $O(p^{\min(A_i + w_i, B_i + v_i)})$. If $a_i$ is an
inexact zero and $y_i \ne 0$, the product is $\delta (y_i + \varepsilon)$
with $\operatorname{ord}_p(\delta) \ge A_i$, of valuation at least
$A_i + w_i = A_i + w_i'$, and $B_i + v_i' = B_i + A_i$ is larger since
$B_i > w_i$. If both factors are inexact zeros, the product has
valuation at least $A_i + B_i$, which is both terms of the minimum. The
sum then follows from the ultrametric inequality.
:::

The convention $v' = A$ for an inexact zero matters: using the stored
valuation $0$ of a zero instead would claim too much precision for
negative $A$ and too little for positive $A$.

`padic_radix_dot_strided_delayed` computes $N_t$ and the minimal
valuation $v_{\min}$ of the nonzero terms in a first pass, and caps the
precision by the context: $N = \min(N_t, a, v_{\min} + r)$ for the
absolute and relative context precisions $a$ and $r$. The cap uses
$v_{\min}$ rather than the valuation of the sum, which is only known
after the sum has been computed; when leading digits cancel, the result
has less relative precision than $r$, as it would with a sequence of
additions. The second pass adds the products modulo $p^N$ into three-word
accumulators, one per residue of the digit offset modulo $e$ and per
sign, so that each limb product costs one multiplication and one
three-word addition and all carries and digit shifts are done once at the
end (as in the dot products modulo $n$ of
[](algorithms/linear-algebra/dot-products)).

If every term is exact, the sum is exact, and it is kept exact (as by
`padic_radix_mul` and `padic_radix_add`) when it fits in the context
precision. To decide this before computing the sum, the dot product uses
the a priori extent $H = \max_i (v_i + e|x_i| + w_i + e|y_i|)$ of the
terms, where $|x|$ is the number of limbs of the unit: if $H$ is at most
one limb beyond the capped precision, the sum is computed in full and
finalised as an exact element (which rounds it if it does not fit after
all); otherwise it is rounded to $N$.

## The limb-aligned form {#limb-aligned}

Within an element, the unit $u$ is stored in radix $B$, with the
valuation $v$ in digits. For block products it is better to align the
values with limbs: write $v = ew + s$ with $0 \le s < e$ and store the
residue $U = u p^{s} \bmod B^L$ together with the *limb frame* $w$. A
product $p^{v_1 + v_2} u_1 u_2$ then equals $U_1 U_2 B^{w_1 + w_2}$, so
all terms of a sum land at whole-limb offsets of a single accumulator and
no digit shifts remain.

::: lemma Window of an output {#lem-window}
Let a sum of products have minimal term valuation $v_{\min}$, capped
precision $N \le v_{\min} + r$ and minimal frame $W = \min_i (f_i + g_i)$
over its terms, where $f_i = \lfloor v_i / e \rfloor$ and
$g_i = \lfloor w_i / e \rfloor$ are the frames of the factors. Then the
digits of the result relative to $p^{eW}$ that are needed lie in
$[0, N - eW)$, and $N - eW \le r + 2e - 2$.
:::

::: proof
For each term $v_i + w_i \ge v_{\min}$, and $e(f_i + g_i) \ge v_i + w_i
- 2(e - 1)$. Hence $eW \ge v_{\min} - 2(e - 1)$ and
$N - eW \le (N - v_{\min}) + 2e - 2 \le r + 2e - 2$.
:::

An output therefore needs at most $\lceil (r + 2e - 2)/e \rceil$ limbs
relative to its frame, and its digits relative to the frame are those of
the sum of the residue products modulo $B^s$ for its window of $s$ limbs.
Any algorithm that computes the sums modulo $B^s$ can therefore compute
the values: the packed classical kernels accumulate the limb products
directly (in registers when every window has at most two limbs), and the
other kernels below reduce modulo a power of $B$. Outputs whose terms are
all exact are computed in a window covering their extent plus a carry
limb and a sign limb, so that a residue $R$ modulo $B^s$ with a nonzero
top limb is decoded as the negative sum $R - B^s$.

## Matrix products {#matrices}

For $C = AB$ the error term of $C_{ij}$ is the minimum over $t$ of the
error terms of $A_{it} B_{tj}$, so the precisions of all entries form a
min-plus matrix product of the exponent arrays (valuations, precisions
and frames, stored as separate arrays). The packed kernel
`padic_radix_mat_mul_packed` computes it with $O(mnk)$ minimum
operations next to the arithmetic, and uses the exact frame $W_{ij}$ of
every output.

For larger matrices the values are computed by integer or modular matrix
multiplication. Row $i$ of $A$ is scaled to the frame $r_i = \min_t
w(A_{it})$ and column $j$ of $B$ to $c_j = \min_t w(B_{tj})$, so that
$C_{ij} = (A' B')_{ij} B^{r_i + c_j}$ for the integer matrices $A'_{it}
= U_{it} B^{w(A_{it}) - r_i}$ and $B'$ likewise, truncated modulo
$B^{S}$ for the widest window $S$. The product is computed by
`nmod_mat_mul` modulo $B$ when every window has one limb, by
`mpn_mod_mat_mul` modulo $B^S$ when $B^S$ has at most 16 binary limbs
(classical multiplication with delayed reduction, Waksman's algorithm or
multimodular multiplication, see
[](algorithms/linear-algebra/matrix-multiplication)), and by
`fmpz_mat_mul` otherwise, followed by a reduction of each entry modulo
$B^{s_{ij}}$. Converting between radix $B$ and binary costs $O(mk + kn +
mn)$ conversions, negligible against the product for the dimensions
where these kernels are used.

The precision bookkeeping is reduced by the following observation.

::: lemma Blocked minima {#lem-blocked}
Split the inner index into blocks. If the exponents $x_{it}$ of row $i$
of $A$ are constant on each block $\beta$, with value $x_{i\beta}$, then
$\min_t (x_{it} + y_{tj}) = \min_\beta (x_{i\beta} + \min_{t \in \beta}
y_{tj})$, and symmetrically for the columns of $B$.
:::

The proof is immediate. In general the right-hand side with $\min_{t
\in \beta} x_{it}$ is a lower bound. The large-matrix routes compute the
blockwise minima of the rows and columns ($O((m + n) k)$ operations),
use them for all outputs whose row or column is block-constant (in
$O(mn)$ operations per block; matrices with uniform valuations and
precisions qualify), and compute the exact $O(k)$ minimum only for the
other outputs. The results are identical to those of the dot products.

A cheap upper bound for the windows decides between the routes before
any of this work: the valuation of an output is at most $\max_{t \in
\beta} v(A_{it}) + \max_{t \in \beta} v(B_{tj})$ for every block
$\beta$, which bounds the windows by separable quantities. The modular
routes win over the packed kernel from about dimension 8 for one-limb
windows (modulo $B$), 32 at two limbs, 16 at 3 to 8 limbs and about 30
to 50 beyond; spread-out frames and irregular exponents raise these
crossovers.

## Polynomial products {#polynomials}

For $c_k = \sum_i a_i b_{k-i}$ the minimal valuation and the error term
of each coefficient are min-plus convolutions
$$
(\alpha \star \beta)_k = \min_i (\alpha_i + \beta_{k-i})
$$
of exponent sequences, with $+\infty$ for absent terms (exact zeros, and
exact coefficients in the error terms). The values are computed by
Kronecker substitution [@Kro1882] in radix $B$, with the residues
relative to the frames of the two polynomials in slots of $2S + O(1)$
limbs and one call to `radix_mulmid`, which uses the small-prime FFT
directly in radix $B$ and needs no conversion to binary; or, without the
padding of the slots, by multiplying the residues modulo $B^S$ with
`_nmod_poly_mulmid` (for $S = 1$) or `_mpn_mod_poly_mulmid` (Karatsuba or
the small-prime FFT, as tuned for `mpn_mod`). The modular route is the
fastest for windows of up to about 10 limbs, from length 4 for one-limb
windows; the radix Kronecker substitution is the fastest beyond.

The min-plus convolution is the difficult part: no algorithm is known
that computes it in time $O(n^{2 - \varepsilon})$ for general integer
sequences, and several problems are equivalent to it [@CMWW2019].
`padic_radix` computes it exactly by combining methods that are fast
in the cases that occur in practice:

::: lemma Convex sequences {#lem-convex}
If $\alpha$ and $\beta$ are convex on intervals of integers (finite
everywhere on them), then $\alpha \star \beta$ is convex, and its
sequence of slopes is the sorted merge of the slopes of $\alpha$ and
$\beta$.
:::

::: proof
This is the discrete form of the Minkowski sum of the epigraphs: walking
from $\alpha_{i_0} + \beta_{j_0}$, the minimum at the next index is
obtained by advancing in the sequence with the smaller next slope, and
by convexity this greedy choice is optimal at every step.
:::

Constant and linear sequences (uniform valuations and precisions, or
precisions growing linearly with the degree) are covered in $O(n)$.

::: lemma Kronecker substitution for minima {#lem-ks-minplus}
Let the finite values of $\alpha$ and $\beta$ lie in $[0, R]$, and let
$F = \sum_i x^i y^{\alpha_i}$ and $G = \sum_j x^j y^{\beta_j}$. Then
$(\alpha \star \beta)_k$ is the lowest power of $y$ that occurs in the
coefficient of $x^k$ of $FG$.
:::

::: proof
The coefficient of $x^k y^d$ in $FG$ counts the pairs $(i, k - i)$ with
$\alpha_i + \beta_{k-i} = d$, a nonnegative integer that is nonzero
exactly for the attained sums.
:::

Packing the counts (at most $n$) into fields of $\lceil \log_2 n \rceil
+ 1$ bits and $2R + 1$ fields per power of $x$ turns this into one
integer multiplication of $O(nR\log n)$ bits. After subtracting a common
linear trend $\lambda i$ from both sequences (which subtracts $\lambda k$
from the result), the range $R$ is small whenever the exponents are a
linear function of the degree plus small fluctuations. In practice this
beats the quadratic algorithm only for $R$ well below $n / \log n$.

For irregular sequences, `padic_radix` splits the values into 32 bands
of equal width $W$ and processes the pairs of bands $(s, t)$ by
increasing $s + t$. All sums from the pair $(s, t)$ are at least
$\alpha_{\min} + \beta_{\min} + (s + t) W$, so an output whose current
value is at most this bound is final, and a pair of bands whose output
range contains only final outputs (or outputs without any finite term,
found by a Kronecker substitution of the supports) is skipped. Each
remaining pair is done directly over its finite entries or by
Kronecker substitution over the actual range of values in the two bands,
whichever is cheaper. Clustered values (block-dependent precisions,
alternating or periodic patterns, sparse sequences with exact zeros, a
few outliers) give few nonempty bands with small internal ranges, and
the early termination skips most pairs; values spread at random give at
worst about the work of the quadratic algorithm. For lengths 4000 to
16000 the measured cost was below 0.1 of the quadratic algorithm for
periodic patterns, alternating patterns with exact zeros and multivariate
patterns, about 0.1 for random values over a large range, and 0.5 to 0.8
for piecewise linear sequences and sequences with a few outliers.

With the context flag `PADIC_RADIX_FAST_BOUNDS`, when the exact
computation would cost more than a Kronecker substitution of range 64,
the result may instead be a lower bound: the maximum of the convex hull
bound of [](#lem-convex) (applied to the lower convex hulls, which is a
lower bound for arbitrary sequences) and the Kronecker substitution of
the values rounded down to multiples of $\delta = \lceil R/64 \rceil$,
which is within $2\delta$ of the minimum. Outputs without finite terms
are still recognised exactly, so exact coefficients stay exact.

Exact precisions matter for structured inputs. Multiplying bivariate
polynomials in $\mathbf{Q}_p[x][y]$ by Kronecker substitution
(`_gr_poly_mulmid_bivariate_KS`) produces sequences with exact zeros as
padding, exact leading coefficients and precisions that depend on the
inner polynomial; with exact min-plus convolutions and exact sums of
exact terms, every coefficient of the product has the same precision as
with the nested classical product, and the exact coefficients remain
exact. For the exact sums the extents of [](#dot) form a max-plus
convolution, computed by the same algorithm on negated sequences.

## Spread exponents {#spread}

The packed classical kernels work in the frame of every output, so their
temporaries are bounded by the output windows, that is, by the precision.
The matrix and polynomial routes through integer or modular arithmetic
instead use base frames shared by many outputs (the row and column
frames $r_i + c_j$, or the frames of the two polynomials), and an output
whose terms all have high valuations then needs a window reaching from
the base frame up to its precision. When the valuations vary by more than
the precision along the outputs (alternating valuations $0, 10^6, 0,
10^6, \ldots$ at 20 digits of precision, valuations growing linearly
with the degree, or valuations at random in a wide range), these windows
are as wide as the spread of the valuations, and the products are
arbitrarily more expensive than the classical algorithm, or impossible.
High valuations are harmless when they lie beyond the precision of every
output they contribute to (for example a few coefficients of huge
valuation among coefficients of valuation zero): they are truncated
away.

The routes therefore compare the window $S_k$ of every output in the
base frames with the largest window $S_{\text{own}}$ of an output in its
own frame, both available after the precision computation. The matrix
routes give up when $\max_k S_k > 4 S_{\text{own}} + 8$, and the
dispatcher then uses the packed kernel. The polynomial routes compute
the outputs with $S_k > 4 S_{\text{own}} + 8$ separately by dot products
when there are few of them (this covers coefficients next to a few
coefficients of huge valuation, such as a polynomial whose first or last
coefficients have huge valuations), and give up otherwise. Since the
min-plus convolutions are not free, the polynomial dispatcher first
estimates the windows at 16 sampled outputs (from the minimal valuations
of their terms and the context precision) when the frames are spread
out, which detects the common spread patterns in $O(n)$ operations. A
product of length $10^4$ with alternating, trending or random valuations
spread over $10^6$ digits at 20 digits of precision takes the time of the
packed kernel (0.2 seconds), and one in which a fifth of the coefficients
have valuation $10^6$ is done by the transform in 0.03 seconds.

With these rules, the temporaries of all routes are bounded in terms of
the precision of the outputs in contexts with a finite precision, as in
ball arithmetic at a given working precision. In a context without
finite precision, an exact result can itself be arbitrarily large (an
exact sum with a huge valuation gap, an exact power with a huge
exponent); such results are bounded only by the `size_limit` of the
context. The delayed dot product then falls back to the chain of
additions when its window exceeds the limit: that chain only builds the
partial sums, and so still handles exact cancellations of terms with
huge valuations.

## Use in FLINT {#use}

- `padic_radix_dot_strided_delayed` implements [](#dot); the generic
  dot products of `padic_radix` use it above a few terms.
- `padic_radix_mat_mul_packed`, `padic_radix_mat_mul_fmpz` and
  `padic_radix_mat_mul_mpn_mod` implement the matrix products, and
  `padic_radix_mat_mul` (the `gr_mat_mul` method) chooses between them.
- `_padic_radix_poly_mulmid_packed`, `_padic_radix_poly_mulmid_KS` and
  `_padic_radix_poly_mulmid_mod` implement the polynomial products, and
  `_padic_radix_poly_mulmid` (the `gr_poly` method) chooses between them.
