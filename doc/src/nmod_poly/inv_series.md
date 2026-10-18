---
title: Power series inversion and division
functions: nmod_poly_inv_series nmod_poly_inv_series_basecase nmod_poly_inv_series_newton nmod_poly_div_series nmod_poly_div_series_basecase
---

A power series $Q = Q_0 + Q_1 x + \cdots$ over $\Z/n\Z$ is invertible if
and only if $Q_0$ is a unit, and then its inverse modulo $x^m$ is a
polynomial of length $m$. These functions compute $1/Q$ and $A/B$ modulo
$x^m$. Inversion is the basic step of division with remainder (applied to
the reversed divisor, see [](nmod_poly/modular-arithmetic)) and of most
power series functions, and it costs only a small constant number of
multiplications.

::: summary algorithms/polynomials/newton :::

## Algorithms {#algorithms}

The *basecase* algorithm solves $Q \cdot \mathit{Qinv} = 1$ coefficient
by coefficient:
$$
q_0 = Q_0^{-1}, \qquad q_i = -q_0 \sum_{j=1}^{\min(i, \mathit{Qlen}-1)}
Q_j q_{i-j},
$$
with each sum computed as one dot product (`_nmod_vec_dot_rev`, see
[](nmod_vec/dot)) and a single reduction. It costs $O(m \min(m,
\mathit{Qlen}))$ operations. A binomial $Q_0 + c x^k$ is inverted as a
geometric series. Division uses the analogous recurrence $q_i = B_0^{-1}
(A_i - \sum_{j \ge 1} B_j q_{i-j})$.

*Newton iteration* doubles the number of correct coefficients in each
step: if $g \equiv 1/Q \pmod{x^k}$, then
$$
g' = g + g\,(1 - Qg) \equiv 1/Q \pmod{x^{2k}}.
$$
The product $Qg$ is known to be $1$ in its low $k$ coefficients, so only
its middle coefficients are needed, which is a
[middle product](algorithms/polynomials/middle-product), and the
correction is a truncated product of length $k$. Summed over all steps
the cost is a small multiple of one multiplication of length $m$. Division
$A/B$ is computed in a similar way, folding the multiplication by $A$ into
the last Newton step instead of computing $1/B$ first, as
proposed by Karp and Markstein [@KM1997].

`nmod_poly_inv_series` uses the basecase algorithm for $\mathit{Qlen} \le
10$, and otherwise for $m$ below a cutoff tabulated by the bit size of
$n$ (with `fft_small`: 36 to 50 for moduli up to 7 bits, 81 to 292 for 8
to 21 bits, 388 to 407 for 22 to 28 bits and 231 to 448 above) or for
short $Q$, and Newton iteration otherwise. `nmod_poly_div_series` uses
the basecase for $\mathit{Blen} \le 5$ and analogous tables (about 80 to
600). `nmod_poly_inv_series_newton` is the same function as
`nmod_poly_inv_series`.

## Examples {#examples}

The generating function of the Fibonacci numbers is $x / (1 - x - x^2)$.
One series division gives $F_0, \ldots, F_{m-1}$ modulo $p$, and here
$F_{100000} \bmod p$ for $p = 10^9 + 7$, checked against the
matrix-free identity $F_{2k} = F_k (2 F_{k+1} - F_k)$:

```c example name=nmod_poly_fibonacci_series
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t num, den, f;
    slong m = 100001;
    ulong a, b;
    nmod_t mod;

    nmod_init(&mod, 1000000007);
    nmod_poly_init_mod(num, mod);
    nmod_poly_init_mod(den, mod);
    nmod_poly_init_mod(f, mod);

    nmod_poly_set_coeff_ui(num, 1, 1);                  /* x */
    nmod_poly_set_coeff_ui(den, 0, 1);                  /* 1 - x - x^2 */
    nmod_poly_set_coeff_ui(den, 1, mod.n - 1);
    nmod_poly_set_coeff_ui(den, 2, mod.n - 1);

    nmod_poly_div_series(f, num, den, m);

    flint_printf("F_0..F_10:");
    for (slong i = 0; i <= 10; i++)
        flint_printf(" %wu", nmod_poly_get_coeff_ui(f, i));
    flint_printf("\n");

    a = nmod_poly_get_coeff_ui(f, 50000);
    b = nmod_poly_get_coeff_ui(f, 50001);
    flint_printf("F_100000  = %wu\n", nmod_poly_get_coeff_ui(f, 100000));
    flint_printf("identity  = %wu\n",
        nmod_mul(a, nmod_sub(nmod_add(b, b, mod), a, mod), mod));

    nmod_poly_clear(num);
    nmod_poly_clear(den);
    nmod_poly_clear(f);
    return 0;
}
```
```output
F_0..F_10: 0 1 1 2 3 5 8 13 21 34 55
F_100000  = 911435502
identity  = 911435502
```

The Bernoulli numbers are defined by $x / (e^x - 1) = \sum_k B_k x^k / k!$.
Inverting the series $(e^x - 1)/x = \sum_k x^k / (k+1)!$ modulo a prime
$p$ larger than the length gives $B_k \bmod p$; the program prints $6
B_2$, $30 B_4$, $42 B_6$, $30 B_8$ and $66 B_{10}$, which are $1, -1, 1,
-1, 5$:

```c example name=nmod_poly_bernoulli_series
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_t mod;
    nmod_poly_t h, b;
    slong k, m = 12;
    ulong fac = 1, den[6] = { 0, 6, 30, 42, 30, 66 };

    nmod_init(&mod, 1000003);
    nmod_poly_init_mod(h, mod);
    nmod_poly_init_mod(b, mod);

    /* h = sum x^k / (k+1)! */
    for (k = 0; k < m; k++)
    {
        fac = nmod_mul(fac, k + 1, mod);
        nmod_poly_set_coeff_ui(h, k, nmod_inv(fac, mod));
    }

    nmod_poly_inv_series(b, h, m);

    /* B_k = k! * b_k */
    fac = 1;
    for (k = 1; k < m; k++)
    {
        fac = nmod_mul(fac, k, mod);
        if (k % 2 == 0)
        {
            ulong B = nmod_mul(fac, nmod_poly_get_coeff_ui(b, k), mod);
            ulong t = nmod_mul(B, den[k / 2], mod);
            flint_printf("%wu * B_%wd = %wd\n", den[k / 2], k,
                         (t > mod.n / 2) ? (slong) t - (slong) mod.n : (slong) t);
        }
    }

    nmod_poly_clear(h);
    nmod_poly_clear(b);
    return 0;
}
```
```output
6 * B_2 = 1
30 * B_4 = -1
42 * B_6 = 1
30 * B_8 = -1
66 * B_10 = 5
```
