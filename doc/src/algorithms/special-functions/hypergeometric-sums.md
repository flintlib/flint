---
title: Sums of hypergeometric series with rational parameters
summary: Partial sums of hypergeometric series whose parameters are small rational numbers and whose argument is a real ball — the forward recurrence with a common denominator, rectangular splitting with exact integer block coefficients and decreasing precision, binary splitting for short arguments, and the variant that splits real and imaginary parts.
complexity: N terms at p bits — O(N) multiplications by small integers plus O(N^{1/2}) full multiplications by rectangular splitting; O(M(p) log² p) by binary splitting when the argument has O(1) bits.
order: 42
---

Many special functions and constants reduce to partial sums
$$
S_N(w) = \sum_{n=0}^{N-1} \frac{(a_1)_n \cdots (a_r)_n}{(b_1)_n \cdots
(b_s)_n}\, w^n
$$
with small rational parameters $a_i = p_i/q_i$, $b_j = u_j/v_j$ and a
real argument $w$ of full precision: the error function ($r = 0$, $b =
(3/2)$), the sine and cosine integrals (${}_1F_2$ and ${}_2F_3$ series
and asymptotic series $\sum_n n!\,(i/z)^n$), Bessel functions of
rational order, and the incomplete gamma function. The functions
`arb_hypgeom_sum_fmpq_arb` and `arb_hypgeom_sum_fmpq_imag_arb` of
[`arb_hypgeom`](arb_hypgeom) compute such sums exactly as written (no
truncation error), with $w = z$ or $w = 1/z$. Their usage is described
in [](arb_hypgeom/sums). In each term ratio
$$
\frac{T(n+1)}{T(n)} = w \prod_i \frac{p_i + n q_i}{q_i}
\prod_j \frac{v_j}{u_j + n v_j},
$$
the numerator and denominator are products of small integers: the cost
of the sum is dominated by the operations that involve the
full-precision $w$.

## Forward recurrence {#forward}

`arb_hypgeom_sum_fmpq_arb_forward` keeps the partial sum as a fraction
$R/D$ with an integer-valued ball $D$: with $c_n$ and $d_n$ the integer
numerator and denominator of the term ratio,
$$
T \gets T\, w\, c_n, \qquad R \gets R\,d_n + T, \qquad D \gets D\, d_n,
$$
and a single division $R/D$ at the end. Each step costs one full
multiplication by $w$ and multiplications by word-size integers when
the products of the parameters fit in a word (otherwise by `fmpz`).
This is the fastest method for a handful of terms.

## Rectangular splitting {#rs}

The sum is evaluated backwards as a polynomial in $w^m$ with $m \approx
\sqrt N$, after computing $w, w^2, \ldots, w^m$:
$$
S_N(w) = B_0(w) + T(m)\big(B_1(w) + \tfrac{T(2m)}{T(m)}\big(B_2(w) +
\cdots\big)\big),
$$
where each block $B_k$ is a polynomial of degree less than $m$ in $w$
with rational coefficients. Within a block, up to 8 consecutive term
ratios are combined over a common denominator, so that the block is
computed as one dot product of the powers of $w$ with exact integer
coefficients (`arb_dot_fmpz`, with single-word products when the
parameters are small) followed by one division by an integer; the
multiplication by $w^m$ between blocks is the only full multiplication
per block.

::: algorithm Rectangular splitting for a hypergeometric sum {#alg-rs}
**Input:** the parameters, $w$, $N$, the precision $p$.

**Output:** $S_N(w)$.

1. $m \gets \max(2, \lfloor\sqrt N\rfloor)$; compute $w^2, \ldots,
   w^m$ (with the denominators of the parameters absorbed into $w$).
2. $s \gets 0$; for the indices $k = N-1, N-2, \ldots, 0$ in groups
   that do not cross multiples of $m$ and have at most 8 elements:
   compute the integer coefficients of the group over a common
   denominator, add their dot product with the corresponding powers of
   $w$ to $s$, and divide by the integer denominator (the product of
   the lower-parameter factors of the group); when $k$ reaches a
   multiple of $m$, multiply $s$ by $w^m$.
3. Return $s$.
:::

The cost is about $N/m + m$ full multiplications and $O(N)$
multiplications and divisions of balls by small integers, which are
linear in $p$. When the terms decrease (for example $r < s$ and $|w|$
moderate), the blocks with small terms are computed at a precision
reduced by the logarithm of their magnitude relative to the largest
term, estimated in double precision from $\log_2|w|$ and the growth of
$(a)_n/(b)_n$. This method is due to Smith [@Smi1989] and is analysed
with these refinements in [@Joh2014b; @Joh2014c].

## Binary splitting {#bs}

When $w = P/Q$ itself has few bits, the matrix product form of
[](algorithms/special-functions/hypergeometric#summation) is computed
as a product tree of exact (integer-valued) balls; `arb_hypgeom_sum_fmpq_arb_bs`
puts $z$ into the numerators, or into the denominators if `reciprocal`
is set. The cost is quasi-linear in the size of the result
[@Bre1976; @CC1988]. The default choice uses binary splitting from
8192 bits when $z$ has at most $p/1000$ bits.

## Real and imaginary parts {#imag}

`arb_hypgeom_sum_fmpq_imag_arb` computes $S_N(iw)$ with real arithmetic:
the even and odd terms form the real and imaginary parts, with signs
$i^{2k} = (-1)^k$ and $i^{2k+1} = (-1)^k i$. The forward and
rectangular splitting versions track the two parts separately (in
rectangular splitting, the dot products take every second power of $w$
and coefficient). This saves the complex multiplications that `acb`
arithmetic would perform. Binary splitting, where the leaves are cheap,
uses complex balls with the factor $i$ in the leaf matrices.

## Use in FLINT {#use}

- `arb_hypgeom_sum_fmpq_arb` and `arb_hypgeom_sum_fmpq_imag_arb` with
  their `_forward`, `_rs` and `_bs` variants.
- `arb_hypgeom_si` and `arb_hypgeom_ci` sum their ${}_1F_2$ and
  ${}_2F_3$ series with `arb_hypgeom_sum_fmpq_arb` and their asymptotic
  series with `arb_hypgeom_sum_fmpq_imag_arb`.
- The special cases $\sum z^k/(a)_k$ and $\sum (a)_k z^{-k}$ with $a =
  p/q$ are implemented separately by `_arb_hypgeom_gamma_lower_sum_rs_1`
  and `_arb_hypgeom_gamma_upper_sum_rs_1` for the error function.
- Multiple-precision floating-point versions of hypergeometric
  binary splitting for constants are in [`mp_real`](mp_real) (see
  [](algorithms/floating-point/overview)).
