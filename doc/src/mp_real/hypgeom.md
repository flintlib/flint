---
title: Hypergeometric series
functions: mp_real_hypgeom_int_struct mp_real_hypgeom_series_struct mp_real_hypgeom_series _mp_real_atan_frac_bsplit
---

`mp_real_hypgeom_series` sums a geometrically convergent hypergeometric
series with integer coefficients to a given precision by binary
splitting. It is a general engine: every constant of the module except
$\log 10$ is one call to it, and new constants can be added by giving
the polynomials of a formula.

## The format {#format}

The series is

$$
S = \left(\frac{c_Q + c_P \sum_{k \ge 1} \frac{P(k)}{Q(k)}
\prod_{j=1}^{k-1} \frac{R(j)}{Q(j)}}{c_D}\right)^{\pm 1}
$$

for integer polynomials $P$, $Q$, $R$ and integers $c_P$, $c_Q$, $c_D$.
It is the format of y-cruncher's `SeriesHypergeometric` formula files
[@Yee2025], whose fields `CoefficientP`, `CoefficientQ`,
`CoefficientD`, `PolynomialP`, `PolynomialQ`, `PolynomialR` and
`Power` translate directly. Prefactors that y-cruncher applies outside the series, such as
the $1/\sqrt{10005}$ of the Chudnovsky formula, are left to the caller.
The term and the ratio share the denominator $Q(k)$ and the product
stops at $k - 1$, so most formulas need only polynomials of small degree.

The requirements are: $Q(k) \ne 0$ for $k \ge 1$, $c_D \ne 0$, and
geometric convergence, that is $\deg R \le \deg Q$ with
$|\operatorname{lc}(R)| < |\operatorname{lc}(Q)|$ when the degrees are
equal. Coefficients may have any size (y-cruncher limits them to 64
bits): `mp_real_hypgeom_series` takes them as limb arrays with signs in
`mp_real_hypgeom_int_struct`, and `mp_real_hypgeom_series_int64` as
`int64_t`.

For example, $\pi$ by Chudnovsky's formula [@CC1988] is
$S/\sqrt{10005}$ with power $-1$, $c_P = 1$, $c_Q = 13591409$,
$c_D = 4270934400$,
$P = [-67957045, -2100495856, 23608573992, -57896553024, 39250089648]$,
$Q = [0, 0, 0, -10939058860032000]$ and $R = [-5, 46, -108, 72]$; and
$\operatorname{atan}(p/q)$ is $S$ with $c_P = c_Q = p$, $c_D = q$,
$P = R = [p^2, -2p^2]$ and $Q = [q^2, 2q^2]$ (for
$\operatorname{atanh}(p/q)$, $P = R = [-p^2, 2p^2]$).

## Accuracy and cost {#accuracy}

The result is a ball with about $\bits(n - 1)$ correct bits, so callers
pass one limb more than they need. The number of terms comes from a
rigorous bound for the tail, which is added to the radius; a series
whose tail cannot be bounded, or which does not converge geometrically,
throws. The algorithm is binary splitting with exact leaves of 24 to 32
terms, removal of common factors at the nodes, truncation of the values
to the working precision near the root, and a single division at the
end; the cost is $O(M(N)\log^2 N)$ for $N$ bits. The merges share the
`fft_small` transform of their common operand between two products
([](fft_small/mpn_shared)). The function uses the
threads of `flint_set_num_threads`, and its result does not depend on
their number. Below about 30 limbs the setup, 1 to 2 microseconds,
dominates.

::: summary algorithms/floating-point/hypergeometric-series :::

`_mp_real_atan_frac_bsplit` and `_mp_real_atanh_frac_bsplit` are a
specialised splitting of the Taylor series of
$\operatorname{atan}(p/q)$ and $\operatorname{atanh}(p/q)$, used by the
Machin-type formulas of the module.

## Example {#example}

```c example name=mp_real_hypgeom
#include <flint/flint.h>
#include <flint/mp_real.h>
#include <flint/arb.h>

int main(void)
{
    /* Chudnovsky: pi = S / sqrt(10005) */
    static const int64_t P[] = { -67957045, -2100495856, 23608573992,
                                 -57896553024, 39250089648 };
    static const int64_t Q[] = { 0, 0, 0, -10939058860032000 };
    static const int64_t R[] = { -5, 46, -108, 72 };
    /* atan(1/5) = S with p = 1, q = 5 */
    static const int64_t P2[] = { 1, -2 }, Q2[] = { 25, 50 };
    mp_real_t s, t;
    arb_t a;
    slong n = 5;

    mp_real_init(s);
    mp_real_init(t);
    arb_init(a);

    mp_real_hypgeom_series_int64(s, -1, 1, 13591409, 4270934400,
        P, 5, Q, 4, R, 4, n);
    mp_real_rsqrt_ui(t, 10005, n);
    mp_real_mul(s, s, t, n);
    mp_real_get_arb(a, s);
    flint_printf("pi         = ");
    arb_printn(a, 60, 0);
    flint_printf("\n");

    mp_real_hypgeom_series_int64(s, 1, 1, 1, 5, P2, 2, Q2, 2, P2, 2, n);
    mp_real_get_arb(a, s);
    flint_printf("atan(1/5)  = ");
    arb_printn(a, 60, 0);
    flint_printf("\n");

    mp_real_clear(s);
    mp_real_clear(t);
    arb_clear(a);
    flint_cleanup_master();
    return 0;
}
```

```output
pi         = [3.14159265358979323846264338327950288419716939937510582097494 +/- 4.60e-60]
atan(1/5)  = [0.197395559849880758370049765194790293447585103787852101517689 +/- 5.98e-62]
```
