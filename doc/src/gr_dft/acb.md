---
title: Transforms of complex ball vectors
functions: gr_dft_acb_pre_t gr_dft_acb gr_dft_acb_precomp_init gr_dft_acb_prod gr_dft_acb_prod_precomp_init gr_dft_acb_dirichlet_index
---

These functions compute DFTs of vectors of complex balls (`acb`) with the
default root $e^{-2\pi i/n}$, in the convention of `acb_dft`. Their output
balls contain the exact transform of every vector of complex numbers in
the input balls. Instead of ball arithmetic, they run the transform on
the midpoints in fixed-point arithmetic with a precision chosen from an
a priori error bound, and add that bound and the effect of the input radii
to the output radii. This is faster than ball arithmetic and, for inputs
of similar magnitude, gives slightly tighter enclosures. `acb_dft`,
`acb_dft_inverse`, the `acb_dft_precomp` functions, `acb_dirichlet_dft`
and `acb_dirichlet_dft_index` are thin wrappers around them.

## Method {#method}

A transform with working precision $\mathit{prec}$ proceeds as follows;
[](algorithms/fft/fixed-point-dft) gives the proofs.

1. *Plan.* A layout of the plan is made in complex mode
   (`_gr_dft_precomp_init_layout`), and `gr_dft_precomp_nfixed_bound`
   gives the bound $\mathit{peak}(1)$ of all intermediate values for
   inputs of modulus 1. With $2^{e_1 - 1} \le \mathit{peak}(1) < 2^{e_1}$,
   inputs scaled to modulus $2^{-e_1 - 1}$ cannot overflow. The bound
   $E$ of the output error in ulp at this scale determines the number of
   limbs $k$: the smallest with $k \bits \ge \mathit{prec}$ and
   $E \, 2^{-k \bits} < 2^{-\mathit{prec} - 2} \cdot 2^{-e_1 - 1}$.
   The layout is then realised over the fixed-point contexts with $k$
   limbs.
2. *Scan.* One pass over the input checks that all midpoints are finite,
   finds the largest absolute value $2^{E_M}$ (upper bound) of a
   midpoint part, and sums all radii into $\rho$.
3. *Transform.* The midpoints are multiplied by $2^{-e}$ with
   $e = E_M + e_1 + 1$, truncated to $k$ limbs, transformed, and
   converted back with multiplication by $2^e$ and rounding to
   $\mathit{prec}$ bits.
4. *Radii.* Each part of each output gets the radius
   $E \, 2^{e - k \bits} + \rho$ plus its rounding error. The term
   $\rho$ bounds the effect of the input radii, since every entry of the
   transform matrix has modulus 1: $|\Delta X_k| \le \sum_j |\Delta x_j|$.
   For the inverse transform, the radius and the midpoints are divided
   by $n$; for powers of two the division is an exact exponent shift,
   otherwise the midpoints are divided during the conversion with a
   single rounding.

If all midpoints are zero, the outputs are zero with radius $\rho$ (or
$\rho/n$). Outputs whose magnitude is close to that of the largest input
have about $\mathit{prec}$ accurate bits; much smaller outputs (from
cancellation) have correspondingly fewer relative bits, but their
enclosures remain correct. The number of guard bits grows linearly with
$\log n$, by about 3 bits per doubling for powers of two, so at 64 bits the
method uses 2 limbs for power-of-two lengths up to at least $2^{16}$.
Lengths with large prime factors, which use Bluestein's algorithm, need
more guard bits.

::: note
Ball arithmetic propagates the radius of every input through every
operation, so the output radii grow with the depth of the transform. The
bound $\rho$ is independent of the algorithm, and a fixed-point
transform has a single roundoff budget, which is why the fixed-point path
usually gives tighter enclosures.
:::

## Falling back to ball arithmetic {#fallback}

The functions use ball arithmetic (a complex-mode plan over `acb` and
`arb`) when the fixed-point path cannot be used:

- the error bound requires more than `GR_DFT_NFIXED_MAX_NLIMBS` limbs,
  that is, $\mathit{prec}$ beyond about $2048 \cdot 64$ bits, or the
  bound computation fails;
- a midpoint is infinite or NaN, or the largest midpoint exponent exceeds
  `WORD_MAX`/4 in absolute value.

With the precomputed functions, the first case is decided once at
initialisation, while inputs of the second kind make a single transform
fall back to a one-shot ball transform. `_gr_dft_acb` with `which` = 1
forces ball arithmetic, and with `which` = 2 it forces the fixed-point
path, returning `GR_UNABLE` and indeterminate outputs where it does not
apply.

## Precomputation and threads {#precomp}

`gr_dft_acb_precomp_init` performs step 1 once: it stores the realised
plan, the fixed-point contexts, the number of limbs and the error bound,
or a ball-arithmetic plan. The power-of-two scale of step 3 is chosen
for each transform from the actual input, so one object serves inputs of
any magnitude. The `prec` argument of the transforms only sets the
rounding of the output midpoints; the accuracy is determined at
initialisation. The product versions do the same for
`gr_dft_acb_prod`.

The transform itself is threaded as described in
[](gr_dft/performance). In addition, when $2n$ is at least 8192 the
scan and the conversions are split over threads; the scan is reduced in a
fixed order, so the results do not depend on the number of threads.

## Products and Dirichlet characters {#prod}

`gr_dft_acb_prod` and its relatives apply the same method to the product
transforms of [](gr_dft/prod), with the bound composed axis by axis by
`gr_dft_prod_precomp_nfixed_bound`. `gr_dft_acb_dirichlet_index` and
`gr_dft_acb_dirichlet` are the Dirichlet transforms of
`gr_dft_dirichlet_index` and `gr_dft_dirichlet`, computed through
`gr_dft_acb_prod`, and implement the transforms of
[](acb_dirichlet#dft).

## Example {#example}

The transform of $x_j = 1/(j + 1)$, $0 \le j < 1000$, at 128 bits, with
the fixed-point method and with ball arithmetic. Here $1000 = 8 \cdot 125$
uses the prime factor algorithm. The fixed-point output is more accurate,
and both enclosures are consistent.

```c example name=gr_dft_acb_accuracy
#include <flint/flint.h>
#include <flint/acb.h>
#include <flint/gr_dft.h>

static slong
min_accuracy(acb_srcptr v, slong n)
{
    slong j, a = WORD_MAX;
    for (j = 0; j < n; j++)
        a = FLINT_MIN(a, acb_rel_accuracy_bits(v + j));
    return a;
}

int main(void)
{
    slong j, n = 1000, prec = 128;
    acb_ptr x = _acb_vec_init(n);
    acb_ptr y = _acb_vec_init(n);
    acb_ptr z = _acb_vec_init(n);

    /* x_j = 1/(j + 1), with rounding radii from the division */
    for (j = 0; j < n; j++)
    {
        acb_one(x + j);
        acb_div_ui(x + j, x + j, j + 1, prec);
    }

    gr_dft_acb(y, x, n, prec);
    GR_MUST_SUCCEED(_gr_dft_acb(z, x, n, 0, 1, prec));   /* balls */
    flint_printf("fixed point: %wd bits, ball arithmetic: %wd bits, %s\n",
        min_accuracy(y, n), min_accuracy(z, n),
        _acb_vec_overlaps(y, z, n) ? "consistent" : "inconsistent");

    /* the round trip contains the input */
    gr_dft_acb_inverse(z, y, n, prec);
    flint_printf("round trip contains input: %d, accuracy %wd bits\n",
        _acb_vec_contains(z, x, n), min_accuracy(z, n));

    acb_printn(y + 1, 30, 0);
    flint_printf("\n");

    _acb_vec_clear(x, n);
    _acb_vec_clear(y, n);
    _acb_vec_clear(z, n);
    return 0;
}
```
```output
fixed point: 123 bits, ball arithmetic: 117 bits, consistent
round trip contains input: 1, accuracy 113 bits
[5.05661018752866038020114767365 +/- 3.06e-30] + [-1.38326640178665733205654167633 +/- 1.27e-30]*I
```
