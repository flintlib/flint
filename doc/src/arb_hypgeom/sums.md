---
title: Hypergeometric sums with rational parameters
functions: arb_hypgeom_sum_fmpq_arb arb_hypgeom_sum_fmpq_imag_arb
---

`arb_hypgeom_sum_fmpq_arb` computes the partial sum
$$
S_N(w) = \sum_{n=0}^{N-1} \frac{(a_1)_n \cdots (a_r)_n}{(b_1)_n \cdots
(b_s)_n} w^n
$$
of a hypergeometric series with rational parameters $a_i, b_j$ given as
`fmpq` arrays and an argument $w$ given as a real ball (or its
reciprocal, if `reciprocal` is set). `arb_hypgeom_sum_fmpq_imag_arb`
computes the real and imaginary parts of $S_N(iw)$ with real
arithmetic only. These are the workhorses behind the real sine and
cosine integrals, the error function and the incomplete gamma function,
and a convenient way to evaluate any function or constant with a series
of this shape: Bessel and Struve functions of rational order,
${}_pF_q$ of rational parameters, $\cos$ and $\sin$ as in the example
below.

The functions sum exactly $N$ terms and add no truncation error. The
caller chooses $N$ and bounds the tail, typically by a geometric
series: if the ratio of consecutive terms is at most $\rho < 1$ in
absolute value for $n \ge N$, the tail is at most $|T_N|/(1 - \rho)$
(see [](algorithms/special-functions/hypergeometric#tail-bound)).

::: summary algorithms/special-functions/hypergeometric-sums :::

## Parameters {#parameters}

- There is no implicit factorial: the usual ${}_pF_q(a; b; w)$ is
  $S_N$ with the lower parameters $b_1, \ldots, b_q, 1$.
- `a` and `b` may be `NULL` if `alen` or `blen` is zero.
- The lower parameters must not be nonpositive integers $-m$ with $m <
  N - 1$, since the terms would divide by zero; an upper parameter that
  is a nonpositive integer $-m$ makes the sum terminate after $m+1$
  terms.
- If `reciprocal` is set, the sum is computed in powers of $w = 1/z$.
  Binary splitting then puts $z$ into the denominators of the term
  ratios instead of computing $1/z$ as a ball. This is how asymptotic
  expansions at a large $z$ are summed (for example $\sum_n (1)_n
  (iz)^{-n}$ in the asymptotic expansion of the sine integral).
- The numerators and denominators of the term ratios are products of
  `fmpz`; they should be small (a few words) for the methods to be
  efficient.

## Algorithms {#algorithms}

`arb_hypgeom_sum_fmpq_arb` chooses between three methods, described in
[](algorithms/special-functions/hypergeometric-sums):

- **forward recurrence** (`_forward`), for few terms: $N \le 2$, $N \le
  8$ up to 1024 bits or $N \le 4$ up to 4096 bits;
- **binary splitting** (`_bs`), from 8192 bits when $z$ has at most
  $\mathit{prec}/1000$ bits, which makes the cost quasi-linear in the
  precision for short arguments;
- **rectangular splitting** (`_rs`) otherwise: the powers $w, \ldots,
  w^m$ with $m \approx \sqrt N$ are computed once and the sum is
  evaluated backwards as a polynomial in $w^m$. Within each block, up
  to 8 consecutive term ratios are combined over a common denominator
  into exact integer coefficients (single-word products when the
  parameters are small), which are applied to the powers of $w$ with
  one `arb_dot_fmpz` and one division by the denominator. Only about
  $N/m + m$ multiplications are between full-precision balls, and when
  the terms decrease, later blocks are computed at a lower precision.

`arb_hypgeom_sum_fmpq_imag_arb` uses the forward recurrence for $N \le
10$ (or $N \le 16$ up to 1024 bits) and otherwise the same choice. Its
rectangular splitting variant forms the same blocks and applies the
even- and odd-indexed coefficients separately (dot products with
stride 2), with the signs of the powers of $i$.

## Example {#example}

The series $e^{ix} = \sum_n (ix)^n/(1)_n$ gives $\cos x$ and $\sin x$
in one call to `arb_hypgeom_sum_fmpq_imag_arb`, and $\cos x =
{}_0F_1(\tfrac12; -x^2/4)$ gives $\cos x$ again with half as many
terms. The tail after $N$ terms of $\sum_n |x|^n/n!$ is at most
$\frac{|x|^N}{N!} (1 - |x|/N)^{-1}$, which bounds both truncations:

```c example name=arb_hypgeom_sum_fmpq_cos_sin
#include <flint/fmpq.h>
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>

int main(void)
{
    arb_t x, c, s, t;
    fmpq b[2];
    mag_t err;
    slong N = 40, prec = 128;

    arb_init(x);
    arb_init(c);
    arb_init(s);
    arb_init(t);
    mag_init(err);
    fmpq_init(b + 0);
    fmpq_init(b + 1);

    arb_set_ui(x, 1);
    arb_div_ui(x, x, 3, prec);

    /* exp(ix) = sum (ix)^n / (1)_n: real part cos(x), imaginary part sin(x) */
    fmpq_one(b + 0);
    arb_hypgeom_sum_fmpq_imag_arb(c, s, NULL, 0, b, 1, x, 0, N, prec);

    /* tail bound |x|^N / N! * 1/(1 - |x|/N) */
    arb_get_mag(err, x);
    mag_pow_ui(err, err, N);
    {
        mag_t u;
        mag_init(u);
        mag_rfac_ui(u, N);
        mag_mul(err, err, u);
        arb_get_mag(u, x);
        mag_div_ui(u, u, N);
        mag_geom_series(u, u, 0);
        mag_mul(err, err, u);
        mag_clear(u);
    }
    arb_add_error_mag(c, err);
    arb_add_error_mag(s, err);
    flint_printf("cos(1/3) = "); arb_printn(c, 35, 0); flint_printf("\n");
    flint_printf("sin(1/3) = "); arb_printn(s, 35, 0); flint_printf("\n");

    /* cos(x) = 0F1(1/2, -x^2/4) = sum (-x^2/4)^n / ((1/2)_n (1)_n) */
    fmpq_set_si(b + 0, 1, 2);
    fmpq_one(b + 1);
    arb_sqr(t, x, prec);
    arb_mul_2exp_si(t, t, -2);
    arb_neg(t, t);
    arb_hypgeom_sum_fmpq_arb(c, NULL, 0, b, 2, t, 0, N / 2, prec);
    arb_add_error_mag(c, err);   /* the same tail bound is valid */
    flint_printf("cos(1/3) = "); arb_printn(c, 35, 0); flint_printf("\n");

    arb_cos(c, x, prec);
    flint_printf("arb_cos  = "); arb_printn(c, 35, 0); flint_printf("\n");

    arb_clear(x);
    arb_clear(c);
    arb_clear(s);
    arb_clear(t);
    mag_clear(err);
    fmpq_clear(b + 0);
    fmpq_clear(b + 1);
    flint_cleanup();
    return 0;
}
```
```output
cos(1/3) = [0.94495694631473766438828400767588061 +/- 2.16e-36]
sin(1/3) = [0.32719469679615224417334408526762061 +/- 3.95e-36]
cos(1/3) = [0.94495694631473766438828400767588061 +/- 2.16e-36]
arb_cos  = [0.94495694631473766438828400767588061 +/- 2.16e-36]
```

The ${}_0F_1$ sum with 20 terms consists of the even terms of the
exponential series up to index 40, so its tail is bounded by the same
quantity.
