---
title: Rigorous numerical DFTs in fixed point
summary: How `gr_dft` computes complex DFTs in multi-limb fixed-point arithmetic with an a priori bound on every intermediate value and on the output error, and how this yields ball enclosures for transforms of acb vectors.
complexity: A transform in $k$-limb fixed point costs the multiplications of the chosen algorithm, each $O(M(k))$ word operations; the bounds cost $O(\log n)$ double operations per plan. For $n = 2^d$ the output error is $O(n^2)$ ulp at an input scale of about $1/n$.
order: 65
---

Ball arithmetic [@Hoe2009; @Joh2017] gives rigorous results for any
computation, but for a DFT it is wasteful: every operation propagates a
radius, the radii grow with the depth of the transform, and each
multiplication costs several floating-point operations. Since all roots
of unity have modulus 1, the magnitudes in a DFT are easy to bound in
advance, which makes fixed-point arithmetic a better fit. `gr_dft` runs
transforms over two `gr` contexts of fixed-point numbers, computes a
bound for every intermediate value and for the output error from the
structure of the plan alone, and uses these bounds to produce enclosures
for transforms of `acb` vectors (`gr_dft_acb`, which also implements
`acb_dft`). For floating-point FFTs, error bounds of a
similar kind are classical [@Hig2002, ch. 24; @Per2003].

## Arithmetic {#arithmetic}

With $k$ limbs of $\beta = \bits$ bits, a real fixed-point number is
$t = \pm f\, 2^{-k\beta}$ with an integer $0 \le f < 2^{k\beta}$, so
$|t| < 1$, and we write $\mathrm{ulp} = 2^{-k\beta}$. A complex number is
a pair of real ones. All error bounds below are in ulp and, for complex
values, bound the modulus of the error unless stated otherwise.

::: lemma Elementary operations {#lem-ops}
Addition, subtraction and negation are exact when the exact result has
absolute value less than 1; otherwise the operation returns `GR_UNABLE`.
A real product is truncated toward zero with error less than 1 ulp. A
complex product has an error of at most $\kappa$ ulp per part, with
$\kappa = 2$ for $k \le 2$, $\kappa = 4$ for $3 \le k < 22$ (four real
products) and $\kappa = 6$ for $k \ge 22$ (three real products), hence of
modulus at most $K = \sqrt 2\, \kappa$.
:::

The error of a real product is less than 1 ulp, but the analysis allows
2, which covers `flint_mpn_mulhigh_n`; for $k \le 2$ the products are
exactly truncated, and a complex product combines two of them exactly.
The three-product formula from 22 limbs computes
$(x_r + x_i)(y_r + y_i) - x_r y_r - x_i y_i$ with the sums held in an
extra limb, so its imaginary part combines three rounded products.
Overflow is detected from the carries, so a computation that returns
`GR_SUCCESS` has kept every value in $(-1, 1)$. The constants $\pm 1$ are
not representable and are replaced by $\pm(1 - \mathrm{ulp})$; the
transforms never multiply by them, since in complex mode the powers
$\pm 1$, $\pm i$ of the root are applied by exchanging and negating
parts.

## The table of roots {#roots}

The plan's table holds approximations $\hat w_j$ of
$w^j = e^{-2\pi i j/n}$. The primitive root is computed with `arb` at
$k\beta + 64$ bits and truncated, with error at most 2 ulp per part. The
entries up to $L \approx n/8$ (or $n/4$, or $n/2$, if $4 \nmid n$ or
$2 \nmid n$) are computed as
$\hat w_j = \hat w_{\lfloor j/2 \rfloor} \hat w_{\lceil j/2 \rceil}$, and
the others are exact copies with exchanged and negated parts, by the
symmetries $w^{n/4 - j} = -i\, \overline{w^j}$,
$w^{n/2 - j} = -\overline{w^j}$ and $w^{n - j} = \overline{w^j}$.

::: proposition {#prop-roots}
Let $E(j)$ bound the error of $\hat w_j$, with $E(1) = 2$ and
$E(j) = E(\lfloor j/2 \rfloor) + E(\lceil j/2 \rceil) + K$. Then
$E(j) \le (2 + K) j - K$ for all $j \ge 1$.
:::

::: proof
For $j = 1$ both sides are 2. For $j \ge 2$, with $a = \lfloor j/2 \rfloor$
and $b = \lceil j/2 \rceil$, induction gives
$E(j) \le (2 + K)(a + b) - 2K + K = (2 + K) j - K$.
:::

The recurrence holds up to second-order terms, because
$|\hat w_a \hat w_b - w^a w^b| \le |\hat w_a - w^a| |\hat w_b| +
|\hat w_b - w^b|$ and $|\hat w_b| \le 1 + E(b)\,\mathrm{ulp}$; the
implementation inflates the bound by a factor $1 + 10^{-7}$ to cover
them. The resulting bound $D_w$ for the whole table, about $(2 + K) L$,
is stored in the plan and is $O(n)$ ulp. Sub-plans whose roots are powers of the parent
root take strided copies of the parent table and its bound.

## Propagating bounds through a plan {#propagation}

The bound computation `gr_dft_precomp_nfixed_bound` follows the plan and
tracks two numbers through the transform: a bound $B$ for the modulus of
the current values and a bound $D$ for the modulus of their errors, in
ulp. A running maximum, the *peak*, records the largest modulus reached
by any input, intermediate or output value.

::: lemma Rules {#lem-rules}
For values $a$, $b$ with bounds $(B_a, D_a)$, $(B_b, D_b)$ and a table
root $\hat w$ with error at most $D_w$:

1. $a \pm b$ has bounds $(B_a + B_b, D_a + D_b)$;
2. the computed $\hat w a$ has bounds $(B_a, D_a + D_w + K)$, if
   $B_a < 1$, and its computation involves values of modulus at most
   $\sqrt 2\, B_a$.
:::

::: proof
The first rule holds since additions are exact. For the second,
$\hat w \hat a - w a = \hat w (\hat a - a) + (\hat w - w) a$, where
$|\hat w| \le 1$ up to second-order terms and $|a| \le B_a < 1$, and the
product adds a rounding error of modulus at most $K$. The value bound
$B_a$ is that of $wa$ (the computed value is within $D$ ulp of it). In
the two-multiplication classes $c(1 \pm i)$, the intermediate values
$x_r \pm x_i$ have modulus at most $\sqrt 2\, |x|$.
:::

Every algorithm is a sequence of these steps, so its bounds are affine
functions of the input bounds $(B_0, D_0)$, and for recursive algorithms
they are computed with memoised coefficients in $O(\log n)$ operations
on doubles, each inflated by a factor $1 + 10^{-10}$ to cover the
rounding of the doubles themselves.

::: theorem Radix-2 {#thm-radix2}
For a radix-2 transform of length $n = 2^d$ with inputs bounded by
$(B_0, D_0)$, the outputs are bounded by $(n B_0, D)$ with
$D \le n D_0 + (n - 1)(D_w + K)$, and the peak is at most
$\sqrt 2\, n B_0$.
:::

::: proof
Each level of [](algorithms/fft/dft-algorithms#alg-dif) doubles $B$ and
$D$ (rule 1) and then multiplies by roots (rule 2, applied to all values
for simplicity), so $D_s = 2 D_{s-1} + D_w + K$, which gives the claim by
induction. The peak is reached inside the last rotations.
:::

The other algorithms are handled in the same way: the split-radix
recursion by memoised affine coefficients per depth, the direct prime
kernel of [](algorithms/fft/dft-algorithms#lem-kernel) by bounding its
$p - 1$ differences of modulus $2B$ and their rotations, the prime factor
and four-step algorithms by composing the bounds of their sub-plans,
and the naive algorithm as $n$ rotated terms. Since $D_w = O(n)$, the
output error of a radix-2 transform is $O(n^2)$ ulp for inputs of
modulus about $1/n$, compared with the $O(\log n)$ relative error of a
floating-point FFT with accurate roots [@Hig2002, ch. 24]; the cost is a
few guard bits per doubling of $n$, about 3 in practice.

::: remark
Computing each root to within 2 ulp, as is done over `acb` by a running
product at a higher precision, would make $D_w$ constant and save about
$\log_2 n$ guard bits at the cost of a more expensive plan.
:::

### Bluestein's algorithm {#bluestein}

A generic bound on Bluestein's algorithm would be too weak: the inverse
convolution transform doubles the bound at each of its $\log_2 N$ levels,
starting from values that are already $N$ times the input. Here the
kernel $\hat h$ is scaled by $1/(2N)$, so $|\hat h_k| \le (2n - 1)/(2N)
< 1/2$, and the values inside the inverse transform are bounded directly:
each intermediate value of the unscaled inverse transform at level $s$ is
$2^s$ times a sum of $N/2^s$ aliased values of the exact convolution
$u = z \ast h$, and $|u_k| \le \|z\|_1 \max|h| \le n B_z/(2N)$. All
intermediate values of the inverse transform are therefore at most
$n B_z / 2$, and the final doubling and chirp multiplication give outputs
bounded by $n B_z$, as for any DFT. The errors follow the generic rules.

## From bounds to enclosures {#enclosures}

The bounds are linear in $B_0$, so for unit inputs the plan has a peak
$P_1 \ge n$, and inputs of modulus less than $1/P_1$ cannot overflow.
`gr_dft_acb` uses this as follows, for working precision $p$.

1. With $P_1 < 2^{e_1}$, the input scale is $m = 2^{-e_1 - 1}$, and $E$ is
   the output error bound for inputs bounded by $m$ with errors of 1 ulp
   (the truncation of the input).
2. The number of limbs $k$ is the least with $k\beta \ge p$ and
   $E\, 2^{-k\beta} < 2^{-p-2} m$, so that the roundoff is below
   $2^{-p}$ relative to the largest input.
3. For an input whose midpoint parts are less than $2^{E_M}$ in absolute
   value, the midpoints are multiplied by $2^{-e}$, $e = E_M + e_1 + 1$,
   truncated to $k$ limbs and transformed. The outputs are multiplied by
   $2^e$, rounded to $p$ bits, and get the radius
   $E\, 2^{e - k\beta} + \rho$ plus their rounding error, where $\rho$ is
   the sum of all input radii.

::: theorem Input radii {#thm-radii}
If $|\Delta x_j| \le r_j$ for all $j$, then the transforms of $x$ and
$x + \Delta x$ differ by at most $\sum_j r_j$ in every entry, for the
forward transform and the inverse transform without division by $n$.
:::

::: proof
$|\sum_j \Delta x_j w^{\pm jk}| \le \sum_j |\Delta x_j|$, since
$|w| = 1$.
:::

The output therefore contains the exact transform of every point of the
input balls: the transform of the midpoints is within the roundoff
bound, and moving from the midpoints to any point of the balls changes
the result by at most $\rho$. The bound $\rho$ is computed once with
upward rounding, which is both cheaper and, for inputs with radii of
similar size, sharper than propagating radii through $\log n$ levels of
ball arithmetic. For the inverse transform, the radius is divided by
$n$, exactly by an exponent shift for powers of two and otherwise
together with a single rounded division of each output midpoint.

::: remark
Step 3 bounds each *part* of the scaled input by $m$, so the modulus of a
scaled input can be up to $\sqrt 2\, m$, while the bound $E$ is evaluated
for modulus $m$. This leaves the overflow analysis valid (the peak is at
most $\sqrt 2/2 < 1$) and does not affect the error bounds of most
algorithms, which do not depend on $B_0$. For Bluestein's algorithm the
error bound has a term proportional to $B_0$ (the error of the kernel
times the values it multiplies), and evaluating it at $\sqrt 2\, m$
increases $E$ by 4 to 10 percent for $n$ between 101 and 4099. The
bounds are pessimistic by much larger factors in practice, but the
proof requires the modulus.
:::
