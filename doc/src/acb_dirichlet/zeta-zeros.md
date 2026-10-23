---
title: Zeros of the Riemann zeta function
functions: acb_dirichlet_hardy_z_zero acb_dirichlet_zeta_zero
---

`acb_dirichlet_zeta_zeros` computes the nontrivial zeros
$\rho_n = 1/2 + i\gamma_n$ of $\zeta(s)$, numbered by increasing
imaginary part from $\rho_1 = 1/2 + 14.1347\ldots i$, and
`acb_dirichlet_hardy_z_zeros` computes their ordinates $\gamma_n$, the
zeros of the Hardy Z-function. The results are proven: each zero is
isolated in an interval that contains exactly one sign change of $Z(t)$,
and the zeros between are counted by Turing's method, so that the $n$-th
zero found is the $n$-th zero of $\zeta(s)$ in the critical strip. The
computation thereby verifies that the zeros in the Gram blocks it
examines are simple and lie on the critical line. The mathematics is
described in [](algorithms/special-functions/zeta-zeros), and real root
isolation in ball arithmetic in
[](algorithms/ball-arithmetic/real-root-isolation).

## Two methods {#methods}

**Riemann–Siegel.** `_acb_dirichlet_hardy_z_zeros_rs` locates the
Gram points $g_k$ (where $\theta(g_k) = k\pi$) around the target, uses
Gram's law for $n \le 126$, Rosser's rule for $n \le 13999526$, and
Turing's method beyond, to find intervals $(a, b)$ containing exactly
one zero, then refines each zero with `_acb_dirichlet_refine_hardy_z_zero`
(the Illinois method at moderate precision, a two-stage refinement or
Newton's method at high precision). Every evaluation of $Z(t)$ costs
$O(\sqrt{t})$ operations by the Riemann–Siegel formula, and the zeros are
refined in parallel when several threads are available.

**Platt.** `acb_dirichlet_platt_hardy_z_zeros` evaluates
$\Lambda(t) e^{\pi t/4}$, which has the same sign changes as $Z(t)$, on a
grid of thousands of points at once with fast Fourier transforms, and
interpolates between grid points with a Gaussian-windowed
Whittaker–Shannon formula with rigorous error bounds [@Pla2017]. The
fixed cost of one multi-evaluation is large, but each zero then costs
only a few interpolations. It applies to $10^4 \le n < 10^{37}$.

`acb_dirichlet_hardy_z_zeros` chooses between them with the cost model of
`_acb_dirichlet_hardy_z_zeros_use_platt`: below $n = 10^{11}$ it always
uses the Riemann–Siegel method, and at larger heights it uses Platt's
method when many zeros are requested (from about 158 zeros at
$n = 10^{11}$, 14 at $10^{13}$, 3 at $10^{15}$). A zero that Platt's method
does not deliver to the requested precision is finished with
`_acb_dirichlet_refine_hardy_z_zero_ball`, and any zero it does not find
is computed by the Riemann–Siegel method, so that the result is accurate
to about `prec` bits either way.

## Precision {#precision}

The precision `prec` is relative: $\gamma_n$ is computed with an absolute
error of about $\gamma_n 2^{-\text{prec}}$. Since
$\gamma_n \approx 2\pi n / \log n$, at least $\log_2 n$ bits are needed to
separate neighbouring zeros. The internal working precision of the
evaluations of $Z$ is chosen automatically and does not depend on
`prec` for the isolation.

## Example {#example}

The first three zeros to 128 bits, the zeros $\rho_{10^6}$ and
$\rho_{10^6+1}$, and the number of zeros below height 600270.

```c example name=acb_dirichlet_zeta_zeros
#include <flint/acb.h>
#include <flint/acb_dirichlet.h>

int main(void)
{
    acb_ptr z;
    arb_t t, N;
    fmpz_t n;
    slong k;

    z = _acb_vec_init(3);
    arb_init(t);
    arb_init(N);
    fmpz_init(n);

    fmpz_one(n);
    acb_dirichlet_zeta_zeros(z, n, 3, 128);
    for (k = 0; k < 3; k++)
    {
        acb_printn(z + k, 30, 0);
        flint_printf("\n");
    }

    fmpz_set_ui(n, 1000000);
    acb_dirichlet_zeta_zeros(z, n, 2, 64);
    for (k = 0; k < 2; k++)
    {
        acb_printn(z + k, 20, 0);
        flint_printf("\n");
    }

    arb_set_ui(t, 600270);
    acb_dirichlet_zeta_nzeros(N, t, 64);
    flint_printf("N(600270) = ");
    arb_printn(N, 20, 0);
    flint_printf("\n");

    _acb_vec_clear(z, 3);
    arb_clear(t);
    arb_clear(N);
    fmpz_clear(n);
    return 0;
}
```
```output
0.500000000000000000000000000000 + [14.1347251417346937904572519836 +/- 3.76e-29]*I
0.500000000000000000000000000000 + [21.0220396387715549926284795939 +/- 3.10e-30]*I
0.500000000000000000000000000000 + [25.0108575801456887632137909926 +/- 3.72e-29]*I
0.50000000000000000000 + [600269.6770124449555 +/- 6.27e-14]*I
0.50000000000000000000 + [600270.3010907116986 +/- 9.23e-14]*I
N(600270) = 1000000.0000000000000
```

The count confirms that the zero just below height 600270 is the
$10^6$-th.
