---
title: The generalized hypergeometric function
functions: acb_hypgeom_pfq
---

`acb_hypgeom_pfq` computes
$$
{}_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;z) = \sum_{k=0}^\infty
\frac{(a_1)_k\cdots(a_p)_k}{(b_1)_k\cdots(b_q)_k} \frac{z^k}{k!}
$$
for complex ball parameters and argument, or the regularized function
${}_p\mathbf{F}_q = {}_pF_q/(\Gamma(b_1)\cdots\Gamma(b_q))$ if
`regularized` is nonzero. `arb_hypgeom_pfq` is the real version. The
function is a front end: it recognizes the orders for which this module
has an implementation with analytic continuation and asymptotic
expansions, and sums the series directly with a rigorous tail bound
otherwise.

::: summary algorithms/special-functions/hypergeometric :::

## Dispatch {#dispatch}

| $(p, q)$ | Function | Implemented by |
|---|---|---|
| $(0, 0)$ | $e^z$ | `acb_exp` |
| $(1, 0)$ | $(1-z)^{-a}$ | `acb_pow` |
| $(0, 1)$ | ${}_0F_1(b; z)$ | `acb_hypgeom_0f1` |
| $(1, 1)$ | ${}_1F_1(a; b; z)$ | `acb_hypgeom_m` |
| $(2, 1)$ | ${}_2F_1(a_1, a_2; b; z)$ | `acb_hypgeom_2f1` |
| other | the series | `acb_hypgeom_pfq_direct` |

For the first five orders, all $z$ are supported (with the branch cut
$[1, \infty)$ for $p = q + 1$). For the other orders the function sums
the defining series with `acb_hypgeom_pfq_direct`, choosing the number
of terms with `acb_hypgeom_pfq_choose_n` and bounding the remainder
with `acb_hypgeom_pfq_bound_factor`. This gives

- a finite result for $p \le q$ and all $z$ (the series is entire),
  although large $|z|$ needs many terms and a high precision to
  overcome the cancellation, since no asymptotic expansion is
  implemented for general orders;
- a finite result for $p = q + 1$ and $|z| < 1$; the closer $|z|$ is
  to 1, the more terms are needed, and no analytic continuation to
  $|z| \ge 1$ is implemented for general orders;
- a finite result for $p > q + 1$ only if the series terminates, that
  is, if some $a_i$ is a nonpositive integer; otherwise the series
  diverges and the result is indeterminate.

If an upper parameter is exactly 1, it is removed together with the
$k!$; other parameters that occur among both the $a_i$ and the $b_j$
are not cancelled, which costs time, and the caller can remove them
beforehand. With `regularized` set and a lower parameter that may be a
nonpositive integer, the series is summed with power series arithmetic
(`acb_hypgeom_pfq_series_direct`) so that the vanishing factors
$1/\Gamma(b_j)$ cancel the poles.

## Example {#example}

The dilogarithm is $\operatorname{Li}_2(z) = z\, {}_3F_2(1,1,1;2,2;z)$;
a ${}_3F_1$ with an upper parameter $-3$ is a polynomial even though
$p > q + 1$; and ${}_1F_2(1; -2, 3; z)$ is undefined while its
regularized version is not:

```c example name=acb_hypgeom_pfq_cases
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>

int main(void)
{
    acb_ptr a, b;
    acb_t z, r;

    a = _acb_vec_init(3);
    b = _acb_vec_init(2);
    acb_init(z);
    acb_init(r);

    /* Li_2(1/2) = (1/2) 3F2(1,1,1; 2,2; 1/2) */
    acb_set_ui(a + 0, 1); acb_set_ui(a + 1, 1); acb_set_ui(a + 2, 1);
    acb_set_ui(b + 0, 2); acb_set_ui(b + 1, 2);
    acb_set_d(z, 0.5);
    acb_hypgeom_pfq(r, a, 3, b, 2, z, 0, 64);
    acb_mul(r, r, z, 64);
    flint_printf("z 3F2(1,1,1;2,2;z) = "); acb_printn(r, 18, 0); flint_printf("\n");
    acb_hypgeom_dilog(r, z, 64);
    flint_printf("Li2(z)             = "); acb_printn(r, 18, 0); flint_printf("\n");

    /* a terminating series with p > q + 1 */
    acb_set_si(a + 0, -3); acb_set_ui(a + 1, 2); acb_set_ui(a + 2, 5);
    acb_set_ui(b + 0, 7);
    acb_set_ui(z, 4);
    acb_hypgeom_pfq(r, a, 3, b, 1, z, 0, 64);
    flint_printf("3F1(-3,2,5;7;4)    = "); acb_printn(r, 18, 0); flint_printf("\n");

    /* the same series without termination diverges */
    acb_set_d(a + 0, -2.5);
    acb_set_d(z, 0.125);
    acb_hypgeom_pfq(r, a, 3, b, 1, z, 0, 64);
    flint_printf("3F1(-2.5,2,5;7;z)  = "); acb_printn(r, 18, 0); flint_printf("\n");

    /* regularized: a lower parameter at a pole of the gamma function */
    acb_set_ui(a + 0, 1);
    acb_set_si(b + 0, -2); acb_set_ui(b + 1, 3);
    acb_set_d(z, 0.5);
    acb_hypgeom_pfq(r, a, 1, b, 2, z, 1, 64);
    flint_printf("1F2~(1;-2,3;z)     = "); acb_printn(r, 18, 0); flint_printf("\n");
    acb_hypgeom_pfq(r, a, 1, b, 2, z, 0, 64);
    flint_printf("1F2(1;-2,3;z)      = "); acb_printn(r, 18, 0); flint_printf("\n");

    _acb_vec_clear(a, 3);
    _acb_vec_clear(b, 2);
    acb_clear(z);
    acb_clear(r);
    flint_cleanup();
    return 0;
}
```
```output
z 3F2(1,1,1;2,2;z) = [0.58224052646501251 +/- 5.21e-18]
Li2(z)             = [0.582240526465012506 +/- 2.06e-19]
3F1(-3,2,5;7;4)    = [-501.857142857142857 +/- 2.91e-16]
3F1(-2.5,2,5;7;z)  = nan + nan*I
1F2~(1;-2,3;z)     = [0.00113163791418884985 +/- 2.03e-21]
1F2(1;-2,3;z)      = nan + nan*I
```

The polynomial ${}_3F_1(-3,2,5;7;4) = -3513/7$ is summed exactly up to
rounding. In the regularized ${}_1\mathbf{F}_2(1;-2,3;z)$ the terms with
$k \le 2$ vanish, because $1/\Gamma(-2+k) = 0$, and the remaining
terms sum to $\frac{(1)_3\, z^3}{3!\,\Gamma(6)}\,{}_1F_2(4; 6, 4; z) =
\frac{z^3}{120}\,{}_0F_1(6; z)$.
