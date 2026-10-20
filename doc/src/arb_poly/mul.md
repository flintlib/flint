---
title: Multiplying polynomials with ball coefficients
functions: _arb_poly_mullow arb_poly_mullow _arb_poly_mul _arb_poly_mulmid
---

`arb_poly_mul`, `arb_poly_mullow` and `arb_poly_mulmid` compute the
full product, the low part (a truncated power series product) and an
arbitrary window of coefficients of the product of two polynomials with
real ball coefficients. The result contains the product of every pair of
polynomials in the input balls. Multiplication is the kernel of almost
everything else in [](arb_poly) and [](acb_poly): power series
inversion, division, the exponential and the other Newton iterations,
composition and the special functions of power series all reduce to it.

::: summary algorithms/ball-arithmetic/polynomial-multiplication :::

## Two algorithms {#algorithms}

The *classical* algorithm computes each output coefficient
$c_k = \sum_{i+j=k} a_i b_j$ as a single `arb_dot`, which rounds the
midpoint sum once and bounds the radius of the whole sum at once. It
costs $O(n^2)$ multiplications of `prec`-bit numbers for length $n$ and
gives the tightest enclosures of the two algorithms.

The *block* algorithm reduces the product to products of integer
polynomials, which FLINT multiplies in quasi-linear time by Kronecker
substitution and FFTs. It first replaces $x$ by $2^{-c} x$, with $c$ the
average slope of $\log_2 |a_i|$ and $\log_2 |b_j|$ between the first and
last nonzero midpoints, so that power series whose coefficients decay
(or grow) geometrically become polynomials with coefficients of similar
magnitude. It then splits each operand into maximal runs of consecutive
coefficients whose midpoints, written as integers times a common power
of two, need at most $3\,\mathit{prec} + 512$ bits, and multiplies every
pair of runs exactly with `_fmpz_poly_mulmid`. The radii are bounded by
the same procedure applied to polynomials with nonnegative coefficients:
for $a = a_m + a_r$ and $b = b_m + b_r$ (midpoint and error),
$$|ab - a_m b_m| \le |a_m| \, b_r + a_r (|b_m| + b_r)$$
coefficientwise, where $|a_m|$ denotes the polynomial of absolute
values, and for a square $|a^2 - a_m^2| \le a_r (2|a_m| + a_r)$. These
products of nonnegative polynomials have 30-bit coefficients and are
computed in double precision for short blocks (with a correction factor
$1 + 10^{-9}$ for rounding errors) and over $\Z$ otherwise. Finally the
scaling is undone exactly. The algorithm page proves that the output is
an enclosure.

The block algorithm is as accurate as the classical one for
well-scaled input, such as power series with geometrically decaying
coefficients, where a single block covers each operand. For badly
scaled input it falls back gradually to more blocks: in the worst case,
every coefficient is a block of its own and the algorithm does the same
work as the classical one, with more overhead.

## Algorithm selection {#selection}

`_arb_poly_mullow` uses the classical algorithm if $n = 1$, if $n \le
7$, or if either operand has length at most 7. Otherwise it computes a
cutoff
$$t = \begin{cases} 110 & \mathit{prec} \le 128, \\ \min(\max(10^4 / (\ln \mathit{prec})^3, 8), 60) & \text{otherwise}, \end{cases}$$
multiplied by $1.25$ for a square at $\mathit{prec} \ge 256$ and again
at $\mathit{prec} \ge 4096$, and uses the classical algorithm if $n \le
t$ or twice the length of the shorter operand is at most $t$, and the
block algorithm otherwise. For example, at 1024 bits the block algorithm
is used from about length 30. `_arb_poly_mulmid` applies the same rule
with $\mathit{nhi}$ in place of $n$.

Squaring is detected by pointer equality: the classical algorithm then
computes the symmetric sums with half the multiplications, and the
block algorithm converts the operand only once and uses the bound
$a_r(2|a_m| + a_r)$ for the radius. Passing two different copies of the
same polynomial is correct but slower and gives slightly larger radii.

The functions without underscore truncate the operands to length $n$
and handle zero lengths and aliasing.

## Middle products {#mulmid}

`arb_poly_mulmid(res, A, B, nlo, nhi, prec)` computes the coefficients
of $x^{\mathit{nlo}}, \ldots, x^{\mathit{nhi}-1}$ of $AB$. The block
version first drops the coefficients of each operand that only
contribute to lower terms and then computes the exact middle products of
the blocks, so that it costs about as much as a product of the length of
the window when the operands are much longer. The Newton iterations for
`arb_poly_inv_series`, `arb_poly_exp_series` and `arb_poly_sqrt_series`
compute the unknown half of a product whose low half is known in this
way ([](algorithms/polynomials/middle-product)).

## Examples {#examples}

The square of the truncated series of $e^x$ is the series of $e^{2x}$,
with coefficients $2^k / k!$ that range over more than 400 orders of
magnitude. Both algorithms give about 124 accurate bits for every
coefficient at a precision of 128 bits, and the middle product agrees
with the corresponding coefficients of the truncated product:

```c example name=arb_poly_mul_exp
#include <flint/arb_poly.h>

int main(void)
{
    arb_poly_t a, c, d, m;
    arb_t t;
    slong i, n = 300, prec = 128;

    arb_poly_init(a);
    arb_poly_init(c);
    arb_poly_init(d);
    arb_poly_init(m);
    arb_init(t);

    /* a = exp(x) + O(x^n), with coefficients 1/k! */
    for (i = 0; i < n; i++)
    {
        arb_fac_ui(t, i, prec);
        arb_inv(t, t, prec);
        arb_poly_set_coeff_arb(a, i, t);
    }

    /* a^2 = exp(2x) + O(x^n), with coefficients 2^k/k! */
    arb_poly_mullow_classical(c, a, a, n, prec);
    arb_poly_mullow_block(d, a, a, n, prec);

    for (i = 50; i < n; i += 100)
    {
        arb_fac_ui(t, i, prec);
        arb_inv(t, t, prec);
        arb_mul_2exp_si(t, t, i);     /* 2^k/k! */
        flint_printf("k = %wd: ", i);
        arb_printn(d->coeffs + i, 15, ARB_STR_NO_RADIUS);
        flint_printf("  accuracy %wd, %wd bits, correct: %d\n",
            arb_rel_accuracy_bits(c->coeffs + i),
            arb_rel_accuracy_bits(d->coeffs + i),
            arb_overlaps(d->coeffs + i, t));
    }

    /* coefficients 100, 101, 102 of a^2 */
    arb_poly_mulmid(m, a, a, 100, 103, prec);
    flint_printf("mulmid: %d %d %d\n",
        arb_overlaps(m->coeffs + 0, d->coeffs + 100),
        arb_overlaps(m->coeffs + 1, d->coeffs + 101),
        arb_overlaps(m->coeffs + 2, d->coeffs + 102));

    arb_poly_clear(a);
    arb_poly_clear(c);
    arb_poly_clear(d);
    arb_poly_clear(m);
    arb_clear(t);
    return 0;
}
```
```output
k = 50: 3.70190194189053e-50  accuracy 124, 125 bits, correct: 1
k = 150: 2.49807767793329e-218  accuracy 124, 124 bits, correct: 1
k = 250: 5.59644861483662e-418  accuracy 124, 124 bits, correct: 1
mulmid: 1 1 1
```

The midpoints $1/k!$ span about 2030 bits, which at this precision
would need three blocks per operand, since a block may span at most $3
\cdot 128 + 512 = 896$ bits. The average slope of $\log_2 (1/k!)$
between $k = 0$ and $k = 299$ is $-6.8$, so the algorithm substitutes $x
\to 2^{7} x$; the scaled coefficients $2^{7k}/k!$ lie between $1$ and
$2^{180}$, and each operand forms a single block.

Input radii propagate as in the bound above. Here the product of
$(1 + x)$ with a ball polynomial whose coefficients have radius
$10^{-10}$ gets radii $10^{-10}$, $2 \cdot 10^{-10}$ and $10^{-10}$
(`arb_poly_printd` prints each coefficient with `arb_printd`):

```c example name=arb_poly_mul_radii
#include <flint/arb_poly.h>

int main(void)
{
    arb_poly_t a, b, c;
    mag_t r;

    arb_poly_init(a);
    arb_poly_init(b);
    arb_poly_init(c);
    mag_init(r);

    arb_poly_set_coeff_si(a, 0, 1);
    arb_poly_set_coeff_si(a, 1, 1);            /* a = 1 + x, exact */

    arb_poly_set_coeff_si(b, 0, 2);
    arb_poly_set_coeff_si(b, 1, 3);            /* b = 2 + 3x +/- 1e-10 */
    mag_set_d(r, 1e-10);
    mag_set(arb_radref(b->coeffs + 0), r);
    mag_set(arb_radref(b->coeffs + 1), r);

    arb_poly_mul(c, a, b, 64);
    arb_poly_printd(c, 10);
    flint_printf("\n");

    arb_poly_clear(a);
    arb_poly_clear(b);
    arb_poly_clear(c);
    mag_clear(r);
    return 0;
}
```
```output
[(2.000000000 +/- 1.0000e-10), (5.000000000 +/- 2.0000e-10), (3.000000000 +/- 1.0000e-10)]
```
