---
title: Multiplication with several primes and CRT
summary: How `fft_small` multiplies integers and polynomials whose product coefficients exceed one prime: chunk sizes and the number of primes, the two-prime and multi-prime pipelines, Chinese remaindering with the normalisation folded in, windows computed with wrap-around, sums of products in transform space, precomputed and shared operands, and transposed products for middle products.
complexity: For two $n$-limb integers, transforms of about $128 n / \beta$ points modulo $k$ primes with $\beta \approx 25k$, about $5n$ points in all, so $O(n \log n)$, plus $O(k^2)$ word operations per coefficient for the reconstruction.
order: 50
---

A number-theoretic transform modulo a 50-bit prime $p$ multiplies
polynomials over $\Z/p\Z$. To multiply integers, or polynomials over $\Z$
or $\Z/n\Z$, the product is first computed over $\Z$: its coefficients
are bounded, so they can be computed modulo several primes $p_0, \ldots,
p_{k-1}$ whose product $P$ exceeds the bound, and recovered with the
Chinese remainder theorem. This multimodular approach goes back to
Pollard [@Pol1971] (see also [@Knu1997, §4.3.3] and [@GG2003]); it is the
method of FLINT's `fft_small` module, whose transforms are described in
[](algorithms/fft/small-prime-fft).

The primes are the eight primes of an `mpn_ctx_t`, all of 50 bits, with
$2^{41} \mid p_i - 1$; FLINT's default primes are listed with
`mpn_ctx_struct`. $P_k = p_0 \cdots p_{k-1}$ has 50, 100, 150, 199, 249,
299, 348 and 398 bits for $k = 1, \ldots, 8$.

## Exactness {#exactness}

::: theorem {#thm-exact}
Let $a = \sum_{i < \ell_a} a_i x^i$ and $b = \sum_{j < \ell_b} b_j x^j$
with integer coefficients $0 \le a_i, b_j < 2^\beta$, and let $K$ such
products be added. Every coefficient $c$ of the sum satisfies
$0 \le c < K \min(\ell_a, \ell_b) \, 2^{2\beta}$, so it is determined by
its residues modulo $p_0, \ldots, p_{k-1}$ if
$K \min(\ell_a, \ell_b) \, 2^{2 \beta} \le P_k$. If the products are
added and subtracted and $|c| < P_k / 2$, $c$ is determined as the
representative of its residue class in $(-P_k/2, P_k/2]$.
:::

::: proof
A coefficient of $ab$ is a sum of at most $\min(\ell_a, \ell_b)$ products
$a_i b_j < 2^{2\beta}$. The Chinese remainder theorem determines $c$
modulo $P_k$, and there is only one representative in an interval of
length $P_k$.
:::

The drivers use slightly different bounds, all of this form:

| Operands | Bound $\ell \, 2^{2\beta}$ used |
|---|---|
| integers cut into $\beta$-bit chunks | $\lceil 64 b_n / \beta \rceil 2^{2\beta}$ ($b_n$ limbs of the second operand) |
| polynomials over $\Z/n\Z$, $n$ of $k$ bits | $b_n 2^{2k}$ |
| polynomials over $\Z$ with $\beta_a$- and $\beta_b$-bit coefficients | $b_n 2^{\beta_a + \beta_b + 1}$ (symmetric range) |
| sums of products in transform space | `len_bound` $\cdot 2^{\mathtt{prod\_bits}}$ |

## Reconstruction {#crt}

::: lemma Chinese remaindering with folded normalisation {#lem-crt}
Let $Q_i = P_k / p_i$ and $e_i = Q_i \bmod p_i$. If
$y_i \equiv c \, e_i^{-1} \pmod{p_i}$ with $0 \le y_i < p_i$, then
$s = \sum_{i < k} y_i Q_i$ satisfies $s \equiv c \pmod{P_k}$ and
$0 \le s < k P_k$.
:::

::: proof
Modulo $p_i$ every term but the $i$-th vanishes, and the $i$-th is
$y_i Q_i \equiv c e_i^{-1} e_i = c$. Each term is less than
$p_i Q_i = P_k$.
:::

So $c = s - jP_k$ for the unique $j < k$ with $s - jP_k < P_k$, which the
code finds by repeated subtraction. The factor $e_i^{-1}$ costs nothing:
it is multiplied into the constant $m_i = (e_i 2^L)^{-1} \bmod p_i$ that
the pointwise product applies anyway to undo the factor $2^L$ of the
inverse transform. A context stores, for each $k$, the cofactors $Q_i$
and $P_k$ as multi-limb integers of a common length large enough to hold
$k P_k$ (1, 2, 3, 4, 4, 5, 6 and 7 limbs), and the reduced cofactors
$e_i$. Reconstructing one coefficient costs $k$ products of a limb by a
cofactor of up to $k - 1$ limbs, and a few subtractions: $O(k^2)$ word
operations, specialised at compile time for each $k$.

For two primes the reconstruction is Garner's formula [@Gar1959]: with
$r_i \equiv c \pmod{p_i}$ in $[0, p_i)$,

$$
c = r_0 + p_0 \left( (r_1 - r_0) \, p_0^{-1} \bmod p_1 \right),
$$

which lies in $[0, p_0 p_1)$. The two-prime pipeline computes it in
floating point on vectors of eight, with the normalisations $2^{-L}$
applied in the same step, and then needs only one $50 \times 50$-bit
product and an addition per coefficient.

The result is then converted to the output: for integers, the
coefficients are added into the product with carries; for polynomials
over $\Z/n\Z$, the integer $c$ of up to seven limbs is reduced modulo $n$
(with two- and three-limb reductions for $k \le 3$); for polynomials over
$\Z$, $c$ is lifted to $(-P_k/2, P_k/2]$ and stored as an `fmpz`. With one
prime $p$ and $n \ne p$, $c < p$ is simply reduced modulo $n$.

## Integers {#integers}

To multiply integers of $a_n$ and $b_n$ limbs, both are cut into chunks
of $\beta$ bits, $a = \sum_i a_i 2^{\beta i}$, giving polynomials of
lengths $\ell_a = \lceil 64 a_n / \beta \rceil$ and
$\ell_b = \lceil 64 b_n / \beta \rceil$. Their product
$\sum_k c_k x^k$ evaluated at $x = 2^\beta$ is $ab$. Larger chunks give
shorter transforms but need more primes, since the chunk products have
$2\beta$ bits; roughly $2 \beta + \log_2 \ell_b \le \log_2 P_k$, so
$\beta \approx 25 k$ and the number of transform points per limb of
the product, $64 k / \beta$, is nearly independent of $k$ (it slowly
decreases with $k$). The real trade-off is between the cost of the
transforms, which grows with $k \log \ell$, and that of the conversions,
which is smaller for fewer, smaller chunks.

### The two-prime pipeline {#two-prime}

Up to 40000 limbs (`FLINT_MPN_MUL_NP2_MAX_AN`) FLINT uses $k = 2$ and
$\beta \in [40, 45]$, the depth $d$ being the least with
$\ell_a + \ell_b - 1 \le 2^d$ and $\beta$ then the largest with
$\ell_b 2^{2\beta} < P_2$ (`_fft_small_np2_choose`; $P_2$ has 99 bits).
The chunks are smaller than both primes, so the input conversion is plain
bit slicing (one unaligned 64-bit load per chunk on x86-64); the
reconstruction is
the vectorised Garner step; and each coefficient, below $2^{100}$, is
added at bit offset $k \beta$ into a sliding window of four limbs. In this
range conversion and reconstruction dominate the cost of the wider
pipeline, which this avoids.

### Profiles {#profiles}

Above 40000 limbs the chunk size and the number of primes come from a
table of fifteen *profiles*, each with a vectorised conversion that
reduces the $\beta$-bit chunks modulo all $k$ primes at once, using
tables of $2^j \bmod p_i$ for $j < 256$ and four primes per vector. A
profile applies to second operands of at most `bn_bound` limbs, the
largest $b_n$ with $\lceil 64 b_n / \beta \rceil 2^{2\beta} \le P_k$:

| $k$ | $\beta$ | `bn_bound` | | $k$ | $\beta$ | `bn_bound` |
|---|---|---|---|---|---|---|
| 4 | 84 | 2536637511 | | 6 | 144 | 2337 |
| 4 | 88 | 10380583 | | 7 | 160 | 613493872 |
| 4 | 92 | 42390 | | 7 | 164 | 2456369 |
| 5 | 112 | 32822700 | | 7 | 168 | 9826 |
| 5 | 116 | 132790 | | 8 | 184 | 2177184315 |
| 5 | 120 | 534 | | 8 | 188 | 8689506 |
| 6 | 136 | 144789875 | | 8 | 192 | 34662 |
| 6 | 140 | 582218 | | | | |

For each $k$ the profile with the largest admissible $\beta$ is a
candidate, and the candidate with the smallest estimated cost

$$
\left(1 - \frac{r}{4}\right) k \, d \, t, \qquad
t = \text{the product length rounded up to a multiple of 256}, \quad
2^d \ge t, \quad r = t / 2^d,
$$

is chosen; the factor $1 - r/4$ rewards truncated transforms that fill
little of their length. When the work is spread over several threads,
only profiles whose $k$ is a multiple of the number of threads qualify,
so that each thread transforms the same number of primes. If $b_n$
exceeds even the largest bound, a scalar conversion with a computed chunk
size is used instead. For balanced products with one thread this gives:

| $a_n = b_n$ (limbs) | 500 | 5000 | 40000 | 50000 | $10^5$ | $10^6$ | $10^7$ | $10^8$ |
|---|---|---|---|---|---|---|---|---|
| $k$ | 2 | 2 | 2 | 5 | 5 | 6 | 7 | 6 |
| $\beta$ | 44 | 43 | 41 | 116 | 116 | 136 | 160 | 136 |
| $d$ | 11 | 14 | 17 | 16 | 17 | 20 | 23 | 27 |

### Reassembly {#reassembly}

Coefficient $c_k$ starts at bit $k \beta$ and occupies the limb count of
the Chinese remainder data, so consecutive coefficients overlap and are
added with carries. The reconstruction is split into ranges of whole
blocks of 256 coefficients, one per thread: each writes the limbs that
only its coefficients reach and returns the few limbs that spill into
the next range, which are added afterwards in one serial pass.

## Windows and wrap-around {#windows}

Often only part of a product is needed: the high limbs of an integer
product, the low or middle coefficients of a polynomial product
[](algorithms/polynomials/middle-product). Only the coefficients in the
window, or for integers the chunk coefficients $c_k$ with
$c_{\mathrm{lo}} \le k < c_{\mathrm{hi}}$ whose bits reach the limb
window, need to be reconstructed. Moreover, the transform can be shorter
than the product.

::: theorem Wrap-around {#thm-wrap}
Let $c = ab$ have length $z_n$, and let $2^d \ge \ell_a, \ell_b$. The
cyclic convolution of length $2^d$ of $a$ and $b$ has coefficients
$c_k + c_{k + 2^d}$ for $k < 2^d$. If $z_n \le z_l + 2^d$, these equal
$c_k$ for all $z_l \le k < \min(2^d, z_n)$; so if also $z_h \le 2^d$,
the window $[z_l, z_h)$ of $c$ is computed by a transform of length $2^d$.
:::

::: proof
Modulo $x^{2^d} - 1$, the term $c_k x^k$ for $k \ge 2^d$ becomes
$c_k x^{k - 2^d}$; since $z_n = \ell_a + \ell_b - 1 < 2^{d+1}$, no
coefficient wraps twice.
For $k \ge z_l$ the wrapped coefficient $c_{k + 2^d}$ has index at least
$z_l + 2^d \ge z_n$, so it is zero.
:::

The plans use it when $\max(\ell_a, \ell_b, z_h) \le 2^d \le z_n \le
z_l + 2^d$ for $2^d$ the largest power of two not exceeding $z_n$. For a
middle product of lengths $2m$ and $m$, the window $[m-1, 2m)$ of a
product of length $3m - 1$, this applies when a power of two lies between
$2m$ and $3m - 1$, and the transform then has that length instead of a
truncated length of about $3m$; for a high product of integers it can
halve the transform (`_mpn_ctx_mpn_mul_range` uses it on the
multi-prime pipeline only, above 40000 limbs). Transposed products
([](#transposed)) give the shorter length for middle products without
this condition. For integers the coefficients below
$c_{\mathrm{lo}}$ are dropped. Their contribution to the window is a
carry of at most a few units into its lowest limb, so the window is a
lower approximation, as for every high or middle product of integers in
FLINT.

## Polynomials over $\Z/n\Z$ and $\Z$ {#polynomials}

For polynomials over $\Z/n\Z$ the coefficients are converted to residues
modulo each prime directly (with one reduction if $n > p_i$). The number
of primes is the least $k$ with $b_n 2^{2 \lceil \log_2 n \rceil} \le
P_k$: one up to about 16 bits, two up to about 44 bits and three above,
depending slowly on the length ([](fft_small/poly_mul#primes)); the
lengths from which `nmod_poly` uses these products, with the number of
primes at each cutoff, are tabulated in
[](algorithms/polynomials/multiplication#cutoffs), and
[](fft_small/usage#primes) summarises the choices for integers and
polynomials. There are two shortcuts:

- If $n$ is itself one of the context primes, a single transform modulo
  $n$ gives the result directly.
- If $n$ is a prime of at least 20 bits that `fft_small` can use, that is
  $2^d \mid n - 1$ for the transform length $2^d$ and $n$ passes the
  arithmetic check, and the second operand has at least 100
  coefficients, the transform is computed modulo $n$ itself, with
  twiddle tables computed for the call. This replaces two or three
  primes by one.

For very small moduli `nmod_poly` reduces the work further by Kronecker
substitution [](algorithms/polynomials/kronecker): for $n \le 23$ it
packs two coefficients into one, $a_{2i} + a_{2i+1} M$ with $M = 2^{16}$
or $2^{17}$, provided norm bounds show that the packed product does not
overflow modulo $p_0$, and multiplies the half-length polynomials modulo
$p_0$ with one prime.

For polynomials over $\Z$, signed coefficients are reduced to residues,
and the result is lifted to the symmetric range, which needs one more
bit. With eight primes, products whose coefficients exceed about 397 bits
are beyond reach; FLINT's `fmpz_poly` multiplication then uses Kronecker
substitution or Schönhage–Strassen [](algorithms/fft/schonhage-strassen),
whose pointwise products modulo $2^N + 1$ are themselves computed with
this module ([](fft_small/mulmod_2expp1)).

## Sums of products and fixed operands {#reuse}

The transform is linear, so a sum of products $\sum_j a_j b_j$ can be
formed pointwise and transformed back once: each operand is transformed
once, each product costs one pointwise multiplication, and only one
inverse transform and one reconstruction are needed. The bound of
[](#thm-exact) must then allow $K$ products; for differences of products
the result must be read in the symmetric range, which costs one more
bit. Products carry the normalisation constant $m_i$, so a product must
not be multiplied again; sums of primal operands, on the other hand,
represent sums of the inputs, with correspondingly larger coefficients.
This is the interface of plans and operands
([](fft_small/transform-space)), used by FLINT for complex integer
products ($a_r b_r - a_i b_i$ and $a_r b_i + a_i b_r$ from four
transforms instead of eight), matrix products and the transformed rings
of `gr`.

When one operand $b$ is used many times, its transforms can be kept. A
product with $b$ then costs one forward transform, one pointwise product
and one inverse transform per prime instead of two forward transforms.
FLINT stores such an operand for cyclic convolutions of a fixed length
$2^d$ (`_mul_precomp_init`) and computes each window of a later product
with [](#thm-wrap). Division with remainder by a fixed divisor $b$ of
length $b_n$ uses two stored operands, $b$ and the reversed power series
inverse $B$ of the reversed divisor to precision $B_n$
[@GG2003]:

::: proposition Division with stored operands {#prop-divrem}
Let $a$ have length $a_n$ with $q_n = a_n - b_n + 1 \le B_n$. The
coefficients $[B_n - 1, B_n - 1 + q_n)$ of the product of the top $q_n$
coefficients of $a$ with $B$ are the quotient $q$ of $a$ by $b$, and the
remainder is $r = (a - bq) \bmod (x^N - 1)$ for any $N \ge b_n - 1$.
:::

::: proof
Write $\operatorname{rev}_m f = x^{m-1} f(1/x)$ for polynomials of
length $m$. From $a = bq + r$ with $\deg r < b_n - 1$ we get
$\operatorname{rev} a \equiv \operatorname{rev} b \cdot \operatorname{rev} q
\pmod{x^{q_n}}$, so $\operatorname{rev} q = \operatorname{rev} a \cdot
(\operatorname{rev} b)^{-1} \bmod x^{q_n}$. Since
$\operatorname{rev}_{B_n} B$ is that inverse to precision $B_n \ge q_n$,
reversing back expresses $q$ as the stated window of a product of the
top $q_n$ coefficients of $a$ and $B$. For the remainder, $r$ has length
less than $b_n - 1 \le N$, so it equals its own reduction modulo
$x^N - 1$, which is the reduction of $a - bq$.
:::

The remainder therefore needs only a cyclic product of length $N$, the
least power of two at least $\max(b_n - 1, 256)$, instead of a product of
length $a_n$. The condition $q_n \le B_n$ is the caller's responsibility
(see [](fft_small/precomp)).

### Integers with a common factor {#shared}

For integers the same saving applies to a factor $x$ shared by several
products. `_fft_small_mpn_mul_shared` forms $ax$ and $bx$ from three
forward and two inverse transforms per prime instead of four and two, and
`_fft_small_mpn_mul_precomp` keeps the transform of $x$ for any number of
products, saving one conversion and forward transform per product (about
30% of a two-prime product). The transform of $x$ is computed for the
longest product, and each partner is transformed only to the truncation
of its own product, so a short partner costs less; only the transform
points that the longest product reads need to be stored. In the two-prime
range each product is computed in one pass over its data (conversion,
forward transform, pointwise product and inverse transform modulo one
prime, then the other, then the Garner step). Making these stages
separate passes over all operands, as the plan and operand interface
does, was measured to cost as much as the shared transform saves. The
binary splitting of `mp_real` uses both functions; see
[](fft_small/mpn_shared).

### Transposed products {#transposed}

A middle product is the transpose of a product
([](algorithms/polynomials/middle-product#transpose)). Transposing the
polynomial pipeline ([](#polynomials)) step by step gives a pipeline for it: the
conversion of the operand that slides, the transposed inverse transform,
the same pointwise product with the same normalisation, the transposed
forward transform truncated to the output length, and the same
reconstruction ([](algorithms/polynomials/middle-product#thm-transposed-fft)).
The coefficient bound is that of the product it transposes: an output
coefficient $\sum_i a_i v_{i+j}$ is a sum of at most $a_n$ products, so
the plan of the product of lengths $a_n$ and $z_n$ applies unchanged,
including the transform modulo $n$ itself. `_nmod_poly_mul_tr_mpn_ctx`
implements this with the convolution engine in a transposed mode.

A general window $[z_l, z_h)$ of a product is such a transposed product
of the reversed shorter operand (length $s$) and the longer one, shifted
by $z_l - (s - 1)$ coefficients, at the cost of a product of lengths $s$
and $z_h - z_l$; the drivers use it when the window is shorter than the
longer operand and the transform is strictly shorter than the one of
[](#windows). For integers (`_mpn_ctx_mpn_mul_range`, on both the
two-prime and the multi-prime pipeline) the shift is rounded down so that
it is a whole number of limbs, and the shorter operand is padded with
zero chunks, so that the reconstruction sees the same chunk coefficients
at the same bit positions as for the full product and returns the same
window; the chunk size and primes are those chosen for the full product.
`radix_mulmid_fft_small` does the same with the radix-$B$
reconstruction. With a fixed long operand $v$, the images of $v$ through
the transposed inverse transform can be stored (`_mul_precomp_init_tr`);
each product then costs the transform of the short operand and the
transposed forward transform, of a length about $v_n$
(`_nmod_poly_mul_tr_precomp`, the convolution engine in its second
transposed mode).

## Parallelism {#threads}

The work for different primes is independent until the reconstruction:
the multi-prime drivers transform the operands for different primes on
different threads (one thread per prime when enough are available, and
for large operands a second thread per prime that transforms the other
operand), and then reconstruct ranges of coefficients in parallel. The
two-prime integer pipeline runs on one thread. The thresholds differ by
coefficient type; see [](fft_small/usage#resources).
