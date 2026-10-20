---
title: Dot products of balls
summary: How arb_dot sums many products of balls in one fixed-point accumulator with a single final rounding, bounds all radius contributions with 30-bit word arithmetic, and lowers the working precision to what the radii justify.
complexity: About len multiplications of the midpoint mantissas truncated to the working precision, plus O(len) word operations for the radius; one rounding at the end.
order: 30
---

Dot products $s + \sum_{i} x_i y_i$ are the inner loop of matrix
multiplication, linear solving, polynomial multiplication and series
evaluation, so their cost and accuracy matter more than those of any
other ball operation. A loop of `arb_addmul` has three weaknesses: each
step rounds the midpoint and allocates or resizes an `arf_t`, each step
computes and rounds a radius with `mag_t` operations including special
case checks, and the intermediate sums are computed to the full working
precision even when the radii make most of those bits meaningless.
`arb_dot` [@Joh2017] avoids all three.

## The algorithm {#algorithm}

Let the terms be $t_i = x_i y_i$ with midpoints $m_{x_i} m_{y_i}$ and
let $s$ be the initial value. Write $e_i$ for the exponent of
$m_{x_i} m_{y_i}$ (the sum of the exponents of the factors) and
$\rho_i = |m_{x_i}| r_{y_i} + |m_{y_i}| r_{x_i} + r_{x_i} r_{y_i}$ for the
propagated error of the term ([](algorithms/ball-arithmetic/error-propagation)).

::: algorithm Ball dot product {#alg-dot}
**Input:** balls $s$, $x_i$, $y_i$ ($0 \le i < n$) with exponents in the
inline `fmpz` range, and a precision $p$.

**Output:** a ball containing $s + \sum x_i y_i$ for all values in the
input balls.

1. Scan the input once and compute $E = \max_i e_i$ (the largest
   possible exponent of a midpoint product), $E_r$, the largest exponent
   of any radius contribution, and the smallest exponent $E_{\min}$ of the
   lowest bit of any midpoint product. If an exponent is outside the
   inline range or a value is not finite, fall back to the simple loop.
2. Reduce the precision: $p \gets \min(p, E - E_r + 30, E - E_{\min} +
   30)$. Bits of the result below the radius are not computed, and if all
   products are short, no more bits than they have are used.
3. Allocate a two's complement fixed-point accumulator of $w = p + g +
   \lceil \log_2(n + 1) \rceil + 1$ bits, rounded up to whole limbs, with
   top exponent $E + \lceil \log_2(n + 1) \rceil + 1$; the guard bits $g = 4 +
   \lceil \log_2 n \rceil$ make the truncation errors of the terms small
   compared with the final rounding.
4. For each term, compute the product of the mantissas truncated to the
   accumulator width (one or two limbs inline, otherwise with
   `mpn_mul` or the MPFR short product `mpfr_mulhigh_n`), shift it to its
   position and add or subtract it. Count one unit of the last accumulator
   bit for each term that is truncated (the arithmetic error counter
   $n_a$).
5. In the same pass, add upper bounds for $\rho_i$ to a 64-bit integer
   accumulator $R$ scaled by $2^{E_r - 30}$, using the top 30 bits of
   each midpoint (plus one) and the 30-bit radius mantissas, and adding one
   unit for each truncation of a shifted product (`RAD_ADDMUL`).
6. Round the accumulator once to $p$ bits (`ARF_RND_DOWN`), and set the
   radius to an upper bound for the sum of three errors: the final rounding
   ($< 1$ unit in the last place), $n_a$ units of the accumulator, and $R
   \cdot 2^{E_r - 30}$.
:::

::: proposition {#prop-dot}
[](#alg-dot) returns a valid enclosure. The error of the midpoint is less
than one unit in the last place of the result plus $n_a \le 3n$ units of
$2^{E + \lceil\log_2 (n+1)\rceil + 1 - w}$, which is a few units of
$2^{E - p - 4}$ at most.
:::

::: proof
The accumulator cannot overflow: there are at most $n + 1$ nonzero
terms, each of absolute value below $2^E$, so the exact sum is below $(n
+ 1) 2^E$, and the extra $\lceil \log_2(n + 1) \rceil + 1$ bits at the top
hold it in two's complement. Each term added to the accumulator differs
from the exact midpoint product only by truncation below the last bit of
the accumulator (when the shifted product does not fit, or when a short
product is used, whose error is at most the length of the product in
units of its last limb, which the extra limb of the short product brings
below one unit of the accumulator); each such event is counted in $n_a$.
The radius accumulator adds, for every term, upper bounds of the three
products in $\rho_i$ rounded up after shifting, so $R \cdot 2^{E_r - 30}
\ge \sum \rho_i$. The rounding of the final sum is a truncation with
error below one unit in the last place. By the triangle inequality the
exact value of $s + \sum x_i y_i$ for any choice of points in the input
balls lies within the sum of the three bounds from the computed
midpoint. The bound on $n_a$ holds because each term (and its radius
contribution) is counted at most three times.
:::

Since all products are added exactly up to the accumulator width, the
result is usually much more accurate than with repeated `arb_addmul`,
where each step rounds to $p$ bits and the errors add up to $n$ units of
the largest partial sum. The radius computation uses only word
operations: a 30-bit by 30-bit product, a shift and an addition per
contribution, without any normalisation of intermediate `mag_t` values.

Lowering the precision in step 2 is what makes `arb_dot` fast for inexact
input: if the entries of a matrix are only known to 50 bits, the products
are computed to about 80 bits, whatever `prec` is.

## Variants {#variants}

- `arb_dot_simple` is the loop of `arb_addmul`, used as a fallback for
  huge exponents and special values, and for testing.
- `arb_dot_precise` forms all products exactly with `arf_dot` and rounds
  once; it is slow and serves as a reference.
- `arb_approx_dot` and `arf_approx_dot` run the same midpoint
  accumulation without the radius part, for floating-point algorithms.
- `arb_dot_ui`, `arb_dot_si`, `arb_dot_uiui`, `arb_dot_siui` and
  `arb_dot_fmpz` convert the integer vector into shallow `arb_struct`
  entries (exact balls whose mantissas point to the integers) and call
  `arb_dot`.
- `acb_dot` runs two accumulators, for the real and imaginary parts, over
  the four real products of each complex product, with the same
  precision reduction and error accounting (`_arb_dot_addmul_generic` and
  `_arb_dot_output`), and the integer variants of `acb_dot` are built
  like those of `arb_dot`.

## Use in FLINT {#use}

`arb_dot` and `acb_dot` are used by the classical and block matrix
multiplication ([](algorithms/ball-arithmetic/matrix-multiplication))
and the triangular solvers of `arb_mat` and `acb_mat` (see [](algorithms/linear-algebra/dot-products) for the analogous
techniques over exact rings), by the basecase polynomial multiplication
and series functions of `arb_poly` and `acb_poly`, by the generic ring
contexts for real and complex balls (the `gr` vector dot product of
[gr_domains](gr_domains) calls `arb_dot` directly), and by many special
function implementations, for example theta functions and the Platt
multi-evaluation of the zeta function.

::: functions arb_dot acb_dot arb_approx_dot arf_approx_dot :::
