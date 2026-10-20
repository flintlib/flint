---
title: Multiplying polynomials with ball coefficients
summary: How arb_poly and acb_poly multiply polynomials with ball coefficients by scaling, splitting into blocks of similar magnitude and exact integer polynomial products, with separate products of nonnegative polynomials for the radii; with proofs that the results are enclosures.
complexity: For well-scaled operands of length $n$ at precision $p$, one product of integer polynomials with coefficients of $O(p)$ bits, $O(M(n(p + \log n)))$ bit operations; $O(n^2 M(p))$ for the classical algorithm and in the worst case.
order: 40
---

Let $A = \sum_i a_i x^i$ and $B = \sum_j b_j x^j$ be polynomials whose
coefficients are balls $a_i = [m_i \pm r_i]$, $b_j = [m'_j \pm r'_j]$.
A ball product $C$ must contain $\sum_{i+j=k} \alpha_i \beta_j$ for all
$\alpha_i \in a_i$, $\beta_j \in b_j$ and every $k$. The *classical*
algorithm evaluates each such sum as one ball dot product
([](algorithms/ball-arithmetic/dot-products)), which costs $O(n^2)$
multiplications for length $n$ but gives nearly optimal enclosures.
Fast algorithms (Karatsuba, FFT) cannot simply be run in ball arithmetic:
they rely on cancellations that ball arithmetic does not see, so the
radii of the output grow like the sum of the absolute values of all
intermediate quantities, which for polynomials with coefficients of very
different magnitudes destroys the accuracy entirely
([@vdH2008], [@Hig2002]). FLINT uses the following algorithm
instead: the midpoints are multiplied *exactly* with fast integer
polynomial arithmetic, after a scaling and a splitting into blocks that
keep the integers small, and the radii are bounded by a separate product
of nonnegative polynomials of low precision, following van der Hoeven
[@vdH2008] and Johansson [@Joh2017].

## Error bound for the product {#bound}

Write $A = A_m + E$ with $A_m = \sum m_i x^i$ and $E = \sum \varepsilon_i
x^i$, $|\varepsilon_i| \le r_i$, and likewise $B = B_m + E'$. For a
polynomial $P$, let $|P|$ denote the polynomial whose coefficients are
the absolute values of those of $P$, and write $P \preceq Q$ if every
coefficient of $P$ is at most the corresponding coefficient of $Q$.

::: lemma {#lem-bound}
With the notation above, $A_r = \sum r_i x^i$ and $B_r = \sum r'_j x^j$,
$$|AB - A_m B_m| \preceq |A_m| B_r + A_r (|B_m| + B_r), \qquad
|A^2 - A_m^2| \preceq A_r (2|A_m| + A_r).$$
:::

::: proof
$AB - A_m B_m = A_m E' + E B_m + E E' = A_m E' + E(B_m + E')$, and the
coefficients of a product of polynomials satisfy $|PQ| \preceq |P|\,|Q|$
by the triangle inequality, with $|E| \preceq A_r$, $|E'| \preceq B_r$
and $|B_m + E'| \preceq |B_m| + B_r$. For the square, $A^2 - A_m^2 =
2 A_m E + E^2 = E(2A_m + E)$.
:::

So it suffices to compute the product $A_m B_m$ of the midpoints with
rigorously bounded rounding error, and upper bounds for one or two
products of polynomials with nonnegative coefficients. These products
need only a few bits of relative accuracy, because they bound errors.

## Scaling and blocks {#blocks}

The coefficients of a power series often decay or grow geometrically:
the coefficients of $e^x$ are $1/k!$, those of $1/(1 - 100x)$ are
$100^k$. Representing all of them as integers with a common exponent
would need $O(n \log n)$ or $O(n)$ bits per coefficient. The
substitution $x \to 2^{-c} x$, which multiplies the coefficient of $x^k$
by $2^{-ck}$, is exact in binary floating-point arithmetic and commutes
with multiplication: if $\tilde A(x) = A(2^{-c} x)$, then $\tilde
A(x) \tilde B(x) = (AB)(2^{-c} x)$, and the coefficients of $AB$ are
those of $\tilde A \tilde B$ multiplied by $2^{ck}$.

FLINT chooses $c$ as the nearest integer to the average slope of the
exponents between the first and the last nonzero midpoint of the two
operands: if the first and last nonzero midpoints of $A$ have indices
$a_1 < a_2$ and exponents $e_1, e_2$, and those of $B$ indices $b_1 < b_2$
and exponents $f_1, f_2$, then
$$c = \left\lfloor \frac{(e_2 - e_1) + (f_2 - f_1)}{(a_2 - a_1) + (b_2 - b_1)} + \frac12 \right\rfloor.$$
This is cheap and gives slowly varying scaled coefficients for power
series with similar decay rates; it is not optimal in general.

Next, the scaled midpoints of each operand are split into maximal runs
of consecutive indices, called *blocks*, such that within a block all
nonzero midpoints are integer multiples of a common power of two
$2^{e}$ and have absolute value below $2^{e + H}$, where the *height* $H$
is at most $3p + 512$ for the working precision $p$. A block $t$
starting at index $s_t$ thus defines an integer polynomial $P_t$ with
$\tilde A_m = \sum_t 2^{e_t} x^{s_t} P_t(x)$, exactly, and similarly
$\tilde B_m = \sum_u 2^{e'_u} x^{s'_u} Q_u(x)$. Then
$$\tilde A_m \tilde B_m = \sum_{t, u} 2^{e_t + e'_u} x^{s_t + s'_u} P_t Q_u,$$
where each $P_t Q_u$ is a product of integer polynomials, computed
exactly by `_fmpz_poly_mulmid` (Kronecker substitution, Schönhage–Strassen
or the small-prime FFT, see [](algorithms/polynomials/multiplication-integers)).
Each output coefficient is the sum of the contributions of all block
pairs; each contribution is rounded to $p$ bits when it is added, with
the rounding error added to the radius. Well-scaled operands form a
single block each, and then the midpoint product is one exact integer
polynomial product followed by one rounding per coefficient.

The products of nonnegative polynomials in [](#lem-bound) are
computed in the same way. Their coefficients are `mag_t` numbers with
30-bit mantissas ([](algorithms/ball-arithmetic/representation)); after
the same scaling they are split into blocks of height at most $3 \cdot
30 + 512 = 602$ bits. For a pair of blocks of which the shorter has
fewer than 1000 coefficients, the product is computed in
double-precision floating-point arithmetic, each coefficient being
converted exactly to a double with a fixed exponent offset (the block
height guarantees that neither the operands nor the products overflow or
underflow); otherwise it is computed exactly with integers. The results
are converted to upper bounds with `mag_set_d_2exp_fmpz` and
`mag_set_fmpz_2exp_fmpz`, which round upwards.

::: lemma Floating-point sums of nonnegative products {#lem-double}
Let $s = \sum_{i=1}^{L} x_i y_i$ with nonnegative doubles $x_i, y_i$, $L
\le 1000$, and let $\tilde s$ be the value computed by multiplying and
adding in any order in IEEE double precision with rounding to nearest,
without overflow or underflow. Then $s \le (1 + 10^{-9}) \tilde s$.
:::

::: proof
Each product and each partial sum is rounded once, with a relative error
factor $1 + \delta$, $|\delta| \le u = 2^{-53}$. Since all quantities are
nonnegative, each term $x_i y_i$ enters $\tilde s$ multiplied by a
product of at most $L$ such factors, each at least $1 - u$, so $\tilde s
\ge (1 - u)^L s$. Hence $s \le (1 - u)^{-L} \tilde s \le (1 + 1.2 \cdot
10^{-13}) \tilde s$ for $L \le 1000$. (A fused multiply-add only removes
roundings.)
:::

The implementation multiplies the double results by $1 + 10^{-9}$,
a comfortable margin over [](#lem-double).

::: algorithm Block multiplication {#alg-block}
Input: $A$, $B$ with ball coefficients, a window $[\ell, h)$ of output
indices, precision $p$. Output: balls $c_\ell, \ldots, c_{h-1}$.

1. Drop the coefficients of each operand that only contribute below
   $x^\ell$ or above $x^{h-1}$, and the trailing zero coefficients.
2. If some coefficient is not finite, use the classical algorithm.
3. Compute the scaling exponent $c$.
4. Radii: compute upper bounds for the coefficients $\ell, \ldots, h-1$
   of $|\tilde A_m| \tilde B_r + \tilde A_r(|\tilde B_m| + \tilde B_r)$
   (or of $\tilde A_r(2|\tilde A_m| + \tilde A_r)$ for a square), omitting
   the terms with $A_r = 0$ or $B_r = 0$, by block products as described
   above, and store them as the radii of $c_k$.
5. Midpoints: split $\tilde A_m$, $\tilde B_m$ into blocks of height at
   most $3p + 512$; for each pair of blocks, compute the exact middle
   product of the integer polynomials and add its coefficients, times
   the power of two, to the balls $c_k$ with rounding to $p$ bits.
6. Multiply $c_k$ by $2^{ck}$.
:::

::: theorem {#thm-block}
For every choice of $\alpha_i \in a_i$ and $\beta_j \in b_j$, the
coefficient $\sum_{i+j=k} \alpha_i \beta_j$ of the product lies in the
ball $c_k$ computed by [](#alg-block), for $\ell \le k < h$.
:::

::: proof
Step 1 only removes terms that do not occur in the requested
coefficients. The scaling in steps 3 and 6 is exact. In step 5, the
block representation of the midpoints is exact and the block products
are exact, so the only errors are the roundings when adding the
contributions, which ball addition adds to the radii. The midpoint of
the result after step 5 therefore differs from the coefficient of
$\tilde A_m \tilde B_m$ by at most the radius accumulated in step 5. By
[](#lem-bound), the coefficient of the exact product differs from that
of $\tilde A_m \tilde B_m$ by at most the coefficient of the bound
polynomial, for which step 4 computed an upper bound: the conversions of
midpoints to `mag_t` in step 4 round upwards, the integer products are
exact, the double products are bounded by [](#lem-double), and the final
conversions round upwards.
:::

The block algorithm is about as accurate as the classical one when the
scaled coefficients of each operand vary by much less than $2^{3p}$, as
is typical for power series of analytic functions: each coefficient of
the product is then rounded once, to $p$ bits relative to the largest
term in its sum. For badly scaled operands it degenerates gracefully to
many small blocks, and in the worst case (every coefficient a block of
its own) to the classical algorithm with more overhead.

## Selection and cost {#selection}

`_arb_poly_mullow` uses the classical algorithm if the output or one of
the operands has length at most 7, or if the output length or twice the
length of the shorter operand is at most a cutoff that is 110 for $p
\le 128$ and $10^4/(\ln p)^3$ clamped to $[8, 60]$ above (increased by
25% for squares at $p \ge 256$ and again at $p \ge 4096$). At high
precision, the integer products are long, and the block algorithm wins
earlier.

For operands of length $n$ that form single blocks, the midpoint product
is one integer polynomial product with $O(p + \log n)$-bit coefficients,
that is, $O(M(n(p + \log n)))$ bit operations, and the radius products
cost $O(n^2)$ double operations for $n < 1000$ (fast in practice, since
they are simple dot products) and $O(M(n \log n))$ beyond. The
classical algorithm costs $O(n^2 M(p))$.

## Complex polynomials {#complex}

For complex polynomials, `_acb_poly_mullow` and `_acb_poly_mulmid`
first try `_acb_poly_mulmid_fmpzi`, when `_acb_poly_mulmid_use_fmpzi`
predicts that it is faster (roughly, when the number of coefficient
products per input coefficient times $p^{0.57}$ is at least 240, or 330
for squares). It scales both operands as above, using the larger of
the real and imaginary midpoints for the slope, writes all real and
imaginary midpoints of each operand as Gaussian integers times a single
power of two, and fails if this needs more than $3p + 512$ bits. The
midpoints of the product are then one product of Gaussian integer
polynomials (`_fmpzi_poly_mulmid`), and the radii of the real and
imaginary parts are bounded separately: for $(a_r + i a_i)(b_r + i b_i)$
with radii $\rho$,
$$\operatorname{rad}(\operatorname{Re}) \le |a_r| \rho_{b_r} + \rho_{a_r}(|b_r| + \rho_{b_r}) + |a_i| \rho_{b_i} + \rho_{a_i}(|b_i| + \rho_{b_i})$$
and symmetrically for the imaginary part, which is [](#lem-bound)
applied to $a_r b_r - a_i b_i$ and $a_r b_i + a_i b_r$; the eight (four
for a square) products of nonnegative polynomials are computed by
`_arb_poly_addmulmid_mag`. If the operands are not well-scaled, the
*transposed* algorithm computes the real and imaginary parts of the
product with four real products `_arb_poly_mulmid`, each of which
uses its own blocks. The *Gauss* variant with three real products,
$(a + bi)(c + di) = (ac - bd) + ((a + b)(c + d) - ac - bd) i$, is
faster but loses accuracy when $|a|$ and $|b|$ (or $|c|$ and $|d|$)
differ greatly, since the error of the large products then swamps the
small part; it is available but not used by default.

## Use in FLINT {#use}

- `_arb_poly_mullow_block` and `_arb_poly_mulmid_block` implement
  [](#alg-block); `_arb_poly_mullow_classical` and
  `_arb_poly_mulmid_classical` the classical algorithm with `arb_dot`;
  `_arb_poly_mullow`, `_arb_poly_mulmid` and `_arb_poly_mul` choose
  between them. `_arb_poly_addmulmid_mag` exposes the radius products.
- `_acb_poly_mulmid_fmpzi`, `_acb_poly_mullow_transpose`,
  `_acb_poly_mullow_transpose_gauss` and `_acb_poly_mullow_classical` are
  the complex algorithms, selected by `_acb_poly_mullow`.
- Over the generic ball contexts of [](gr_poly), `_gr_poly_mullow` and
  `_gr_poly_mulmid` call these functions, so that every generic
  polynomial algorithm over `arb` and `acb` coefficients uses them.
- Matrix multiplication uses the same idea of exact integer products of
  scaled blocks ([](algorithms/ball-arithmetic/matrix-multiplication)).
