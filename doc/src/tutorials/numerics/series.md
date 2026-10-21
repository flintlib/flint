---
title: Summing a series with a rigorous error bound
summary: A complete rigorous computation from scratch: truncating an infinite series, bounding the tail, choosing the number of terms and the working precision, and verifying the result; with the dilogarithm at 1/2 and Apéry's series for ζ(3).
order: 3
---

Most constants and special functions are defined by infinite series, and
summing a series is the most common way to compute a number that has no
library function. Ball arithmetic handles the rounding errors of the
summation automatically, but the error of stopping after finitely many
terms is invisible to it and has to be bounded by hand. This tutorial
goes through the complete process for two series:

1. find a bound for the tail $\sum_{k > N} t_k$ that holds for every $N$;
2. choose $N$ so that the tail is smaller than the target accuracy;
3. choose the working precision so that rounding errors do not dominate;
4. sum the terms in ball arithmetic and add the tail bound to the radius;
5. check the result independently.

## Tail bounds {#tails}

Two simple bounds cover a large proportion of the series met in practice.

::: lemma Geometric tail bound {#lem-geometric}
Let $t_k \ge 0$ and suppose that $t_{k+1} \le q\, t_k$ for all $k > N$,
with $q < 1$. Then $\sum_{k > N} t_k \le t_{N+1} / (1 - q)$.
:::

::: proof
By induction $t_{N+j} \le q^{j-1} t_{N+1}$ for $j \ge 1$; sum the
geometric series.
:::

::: lemma Alternating series bound {#lem-alternating}
If $t_k = (-1)^k u_k$ where $u_k \ge 0$ is decreasing with limit 0, then
$\left|\sum_{k > N} t_k\right| \le u_{N+1}$.
:::

::: proof
The partial sums of $\sum_{k > N} t_k$ alternate around the limit and the
distance between consecutive ones is $u_k$, so the limit lies between
$0$ and $t_{N+1}$.
:::

A bound need not be tight: overestimating the tail by a factor 2 costs
one more term at most for a series that converges geometrically. It must
be valid, though, for the actual $N$ used, which is why it should be
proved for all $N$ (or for all $N$ above an explicit threshold that the
program checks).

## The dilogarithm at 1/2 {#li2}

The series
$$
S = \sum_{k=1}^{\infty} \frac{1}{k^2 2^k} = \operatorname{Li}_2(\tfrac12)
  = \frac{\pi^2}{12} - \frac{(\log 2)^2}{2}
$$
has positive terms with ratio $t_{k+1}/t_k = k^2 / (2 (k+1)^2) < 1/2$, so
by [](#lem-geometric) the tail after $N$ terms is at most
$2 t_{N+1} \le 1 / ((N+1)^2 2^N)$. For a target of $p$ bits, $N = p$
terms suffice. The program sums the series at 128 bits with 20 terms,
with and without the tail bound, and with $N = 128$ terms, and compares
each result with the closed form.

```c example name=numerics_series_li2
#include <flint/arb.h>

/* S = sum_{k>=1} 1 / (k^2 2^k), with N terms and a bound for the tail */
void sum_li2_half(arb_t s, slong N, int with_tail, slong prec)
{
    arb_t t;
    mag_t tail;
    slong k;

    arb_init(t);
    mag_init(tail);

    arb_zero(s);
    for (k = 1; k <= N; k++)
    {
        arb_one(t);
        arb_mul_2exp_si(t, t, -k);
        arb_div_ui(t, t, k, prec);
        arb_div_ui(t, t, k, prec);
        arb_add(s, s, t, prec);
    }

    /* sum_{k>N} 1/(k^2 2^k) <= 1/((N+1)^2 2^N) */
    if (with_tail)
    {
        mag_one(tail);
        mag_mul_2exp_si(tail, tail, -N);
        mag_div_ui(tail, tail, N + 1);
        mag_div_ui(tail, tail, N + 1);
        arb_add_error_mag(s, tail);
    }

    arb_clear(t);
    mag_clear(tail);
}

int main(void)
{
    arb_t s, ref, t;
    slong prec = 128;
    slong N = prec;     /* terms < 2^-prec from here on */

    arb_init(s); arb_init(ref); arb_init(t);

    /* the closed form pi^2/12 - (log 2)^2 / 2 */
    arb_const_pi(ref, prec);
    arb_sqr(ref, ref, prec);
    arb_div_ui(ref, ref, 12, prec);
    arb_const_log2(t, prec);
    arb_sqr(t, t, prec);
    arb_mul_2exp_si(t, t, -1);
    arb_sub(ref, ref, t, prec);
    flint_printf("closed form:       "); arb_printn(ref, 30, 0); flint_printf("\n");

    sum_li2_half(s, 20, 0, prec);
    flint_printf("20 terms:          "); arb_printn(s, 30, 0);
    flint_printf("  overlap: %d\n", arb_overlaps(s, ref));

    sum_li2_half(s, 20, 1, prec);
    flint_printf("20 terms + tail:   "); arb_printn(s, 30, 0);
    flint_printf("  overlap: %d\n", arb_overlaps(s, ref));

    sum_li2_half(s, N, 1, prec);
    flint_printf("%wd terms + tail:  ", N); arb_printn(s, 30, 0);
    flint_printf("  overlap: %d\n", arb_overlaps(s, ref));

    arb_clear(s); arb_clear(ref); arb_clear(t);
    flint_cleanup();
    return 0;
}
```

```output
closed form:       [0.582240526465012505902656320160 +/- 3.20e-31]
20 terms:          [0.582240524473457414707293756711 +/- 2.79e-31]  overlap: 0
20 terms + tail:   [0.58224052 +/- 6.64e-9]  overlap: 1
128 terms + tail:  [0.582240526465012505902656320160 +/- 3.20e-31]  overlap: 1
```

The second line shows what goes wrong without the tail bound: the ball
has a radius of $3 \cdot 10^{-31}$ and 30 digits are printed, but only
eight of them are correct. Ball arithmetic proves that this ball contains
the 20th partial sum, which is a different number from $S$. With the
tail bound the result is correct and honest about its accuracy, and with
enough terms it is as accurate as the closed form.

Tail bounds are computed with [](mag) (`mag_t`), the type of radii: a
`mag_t` is an unsigned floating-point number with a 30-bit mantissa, and
all its operations round upwards, so the computed value is an upper
bound for the exact one. `arb_add_error_mag` adds it to the radius. The
`mag` operations are much cheaper than ball operations, and since an
upper bound is exactly what is needed, they are the natural type for
error bounds; for lower bounds, use the `mag_*_lower` functions.

## Apéry's series for ζ(3) {#zeta3}

The series
$$
\zeta(3) = \frac52 \sum_{k=1}^{\infty} \frac{(-1)^{k+1}}{k^3 \binom{2k}{k}}
$$
which Apéry used to prove the irrationality of $\zeta(3)$ [@Apery1979],
converges like $4^{-k}$. Its terms alternate in sign and decrease in
magnitude, so by [](#lem-alternating) the tail is bounded by the
magnitude of the first omitted term, which the program has computed
anyway. Writing $u_k = 1/(k^3 \binom{2k}{k})$, the ratio
$$
\frac{u_{k+1}}{u_k} = \frac{k^3}{2 (k+1)^2 (2k + 1)}
$$
gives each term from the previous one with three multiplications and two
divisions by small integers, which are cheap (`arb_mul_ui` and
`arb_div_ui`). Since $\binom{2k}{k} \ge 4^k / (2k)$, we have
$u_k \le 2 \cdot 4^{-k}$, so $N = p/2 + 1$ terms bring the tail below
$2^{-p}$.

The working precision needs a few guard bits: each of the $N$ additions
and $5N$ multiplications or divisions rounds, and the relative errors
of the term recurrence accumulate. Adding $\log_2 N + 10$ bits is
generous; the final `arb_set_round` rounds the result to the target
precision.

```c example name=numerics_series_zeta3
#include <flint/arb.h>

/* zeta(3) = (5/2) sum_{k>=1} (-1)^(k+1) / (k^3 binomial(2k,k)) */
void zeta3_apery(arb_t res, slong prec)
{
    arb_t s, t;
    slong k, N, wp;

    /* |t_k| <= 2 * 4^-k, so N = prec/2 + 1 terms leave a tail < 2^-prec */
    N = prec / 2 + 1;
    wp = prec + FLINT_BIT_COUNT(N) + 10;   /* guard bits for N roundings */

    arb_init(s);
    arb_init(t);

    arb_set_d(t, 0.5);                      /* t_1 = 1/(1^3 binomial(2,1)) */
    arb_zero(s);
    for (k = 1; k <= N; k++)
    {
        if (k % 2)
            arb_add(s, s, t, wp);
        else
            arb_sub(s, s, t, wp);

        /* t_{k+1} = t_k k^3 / (2 (k+1)^2 (2k+1)) */
        arb_mul_ui(t, t, k, wp);
        arb_mul_ui(t, t, k, wp);
        arb_mul_ui(t, t, k, wp);
        arb_div_ui(t, t, 2 * (k + 1) * (k + 1), wp);
        arb_div_ui(t, t, 2 * k + 1, wp);
    }

    /* alternating series with decreasing terms: |tail| <= |t_{N+1}| */
    arb_add_error(s, t);

    arb_mul_ui(s, s, 5, wp);
    arb_mul_2exp_si(res, s, -1);
    arb_set_round(res, res, prec);

    arb_clear(s);
    arb_clear(t);
}

int main(void)
{
    arb_t x, y;
    slong prec = 3330;          /* about 1000 digits */

    arb_init(x);
    arb_init(y);

    zeta3_apery(x, prec);
    arb_zeta_ui(y, 3, prec);

    arb_printn(x, 1000, ARB_STR_CONDENSE * 10);
    flint_printf("\naccuracy: %wd bits\n", arb_rel_accuracy_bits(x));
    flint_printf("overlaps arb_zeta_ui: %d\n", arb_overlaps(x, y));

    arb_clear(x);
    arb_clear(y);
    flint_cleanup();
    return 0;
}
```

```output
[1.2020569031{...979 digits...}7766181092 +/- 2.29e-1000]
accuracy: 3328 bits
overlaps arb_zeta_ui: 1
```

`arb_add_error(s, t)` adds $|t|$, more precisely an upper bound for
$|x|$ over all $x$ in the ball $t$, to the radius of $s$. Since the ball
$t$ contains $u_{N+1}$, this is a valid bound for the tail even though
$u_{N+1}$ itself is only known approximately.

## Beyond term-by-term summation {#faster}

Both programs perform $O(N) = O(p)$ operations on $p$-bit numbers. Each
of them adds two numbers or multiplies or divides by a small integer, at
a cost of $O(p)$ word operations, so the total cost is $O(p^2)$. That is
adequate up to a few thousand digits. For more digits,
the terms of a hypergeometric series like these (consecutive terms have a
rational ratio) are better summed by *binary splitting*, which computes
the partial sum exactly as a fraction by recursively combining halves of
the sum, so that the work is dominated by a few multiplications of
large integers; this gives quasi-linear complexity [@HP1998]. See
[](algorithms/integers/binary-splitting). The library constants
(`arb_const_pi`, `arb_const_apery` and the others) and the evaluation of
hypergeometric functions in [](arb_hypgeom) use binary splitting and
related techniques, with the same structure of rigorous truncation
bounds as here. For series whose terms are not hypergeometric, faster
convergence can often be obtained by convergence acceleration or by
Euler–Maclaurin summation, as in the zeta function code (see
[](algorithms/special-functions/zeta)).
