---
title: Evaluating Riemann theta functions
functions: acb_theta_jet acb_theta_one
---

`acb_theta_jet` is the general entry point of the module: it evaluates
one or all $2^{2g}$ theta functions $\theta_{a,b}(z, \tau)$, or their
partial derivatives with respect to $z$ up to a total order `ord`, at
several vectors $z$ for one matrix $\tau$. `acb_theta_one` (one value)
and `acb_theta_all` (all values, or their squares) are the common special
cases. The numbering of characteristics and derivatives and the
precision semantics are explained in [](acb_theta/conventions).

## Output layout {#layout}

With $n_{\text{th}} = 2^{2g}$ if `all` is nonzero and 1 otherwise, and
$n_{\text{jet}} = \binom{\mathit{ord} + g}{g}$ (`acb_theta_jet_nb`), the
output `th` has length $\mathit{nb} \cdot n_{\text{th}} \cdot n_{\text{jet}}$:
for each of the `nb` vectors $z$ in turn, for each characteristic in
increasing order, the vector of derivatives in the order of
`acb_theta_jet_tuples`. When `all` is zero, the characteristic is `ab`;
when it is nonzero, `ab` is ignored.

With `ord` = 0 and `sqr` nonzero, the squares $\theta_{a,b}(z,\tau)^2$
are computed. For `all`, this is faster than computing the values: by the
duplication formula
$\theta_{a,b}(z,\tau)^2 = \sum_{a'} (-1)^{a'^T b}\theta_{a',0}(0,2\tau) \theta_{a+a',0}(2z,2\tau)$
only the $2^g$ functions $\theta_{a,0}$ at $2\tau$ are needed, and no
square roots have to be chosen at the last step. Squares are what most
applications to modular forms need (for instance
`acb_theta_g2_even_weight`). `sqr` is ignored when `ord` is positive.

## How the evaluation works {#method}

1. **Reduction of $\tau$.** `acb_theta_reduce_tau` finds
   $m \in \mathrm{Sp}_{2g}(\Z)$ with $m \cdot \tau$ reduced
   (`acb_siegel_reduce`: LLL reduction of $Y$, reduction of the real
   part, and fundamental matrices), and replaces each $z$ by
   $(\gamma\tau + \delta)^{-T} z$.
2. **Reduction of $z$.** `acb_theta_reduce_z` subtracts from each $z$ a
   lattice vector $\tau r + s$ with $r, s$ even integer vectors, so that
   $Y^{-1}y$ and $\operatorname{Re}(z)$ are bounded; this multiplies the
   theta values by explicit exponentials.
3. **Evaluation at the reduced point.** `acb_theta_jet_notransform` uses
   the duplication formulas above, or
   $\theta_{0,0}(z,\tau) = \sum_a \theta_{a,0}(2z, 4\tau)$ for a single
   value, and calls `acb_theta_ql_jet`. This evaluates the $2^g$
   functions $\theta_{a,0}$ (or all $\theta_{a,b}$) by the quasi-linear
   algorithm of `acb_theta_ql_exact`, with finite differences for
   derivatives, or by direct summation at low precision.
4. **Transformation formula.** The values at the reduced point are
   permuted (`acb_theta_char_table`) and multiplied by the 8th roots of
   unity, the square root of $\det(\gamma\tau + \delta)$
   (`acb_siegel_kappa`) and the exponential factors; derivatives are
   transformed with `acb_theta_jet_compose` and `acb_theta_jet_mul`.

If a reduction fails (because $Y$ cannot be shown to be positive
definite, or the entries are too large for the working precision), the
output is indeterminate.

At a fixed reduced $\tau$, step 3 costs $O(M(p) \log p)$ operations at
precision $p$ for fixed $g$ and `ord` = 0, where $M(p)$ is the cost of a
multiplication; the implied constant grows like $2^{2g}$ or faster with
the genus. The algorithm is described in
[](algorithms/special-functions/theta-quasi-linear). Computing all
$2^{2g}$ values is not much slower than computing one, and several
vectors $z$ share the reduction of $\tau$ and parts of the work.

## Precision {#precision}

The argument `prec` controls the absolute error of the normalised value
$\widetilde\theta_{a,b}(z,\tau) = e^{-\pi y^T Y^{-1} y}\theta_{a,b}(z,\tau)$,
which is bounded on $\mathbb{C}^g$; see [](acb_theta/conventions#precision). The
output values themselves therefore have absolute errors of about
$e^{\pi y^T Y^{-1} y} 2^{-\mathit{prec}}$ (after reduction), and values
much smaller than the typical size, such as the odd theta constants that
vanish, have fewer correct bits relative to their size. The input balls
may be inexact: their radii are propagated with bounds on the
derivatives (`acb_theta_ql_jet_error`).

## Example {#example}

The program computes the 16 theta constants for the genus 2 matrix
$\tau = \begin{pmatrix} i & 1/2 \\ 1/2 & 2i \end{pmatrix}$, the first
derivatives of the odd function $\theta_{(1,0),(1,0)}$ at $z = 0$, and
checks the relation $\theta_{1,1} = -\theta_1$ with
[`acb_modular`](acb_modular) in genus 1. The odd theta constants vanish,
and the computed balls contain zero.

```c example name=acb_theta_jet
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/acb_modular.h>
#include <flint/acb_theta.h>

int main(void)
{
    acb_mat_t tau;
    acb_ptr z, th;
    acb_t t1, t2, t3, t4, w;
    slong k, prec = 128;

    acb_mat_init(tau, 2, 2);
    z = _acb_vec_init(2);
    th = _acb_vec_init(16 * 3);
    acb_onei(acb_mat_entry(tau, 0, 0));
    acb_set_d(acb_mat_entry(tau, 0, 1), 0.5);
    acb_set_d(acb_mat_entry(tau, 1, 0), 0.5);
    acb_set_d_d(acb_mat_entry(tau, 1, 1), 0.0, 2.0);

    acb_theta_all(th, z, tau, 0, prec);
    for (k = 0; k < 16; k++)
    {
        flint_printf("theta_%wu%wu%wu%wu = ", (k >> 3) & 1, (k >> 2) & 1,
                     (k >> 1) & 1, k & 1);
        if (acb_theta_char_is_even(k, 2))
            acb_printn(th + k, 15, ARB_STR_NO_RADIUS);
        else
            flint_printf("odd, contains zero: %d", acb_contains_zero(th + k));
        flint_printf("\n");
    }

    /* value and first derivatives of theta_{(1,0),(1,0)} at z = 0 */
    acb_theta_jet(th, z, 1, tau, 1, 10, 0, 0, prec);
    flint_printf("jet of theta_1010: ");
    for (k = 0; k < 3; k++)
    {
        acb_printn(th + k, 12, ARB_STR_NO_RADIUS);
        flint_printf(k < 2 ? ", " : "\n");
    }

    /* genus 1: theta_{1,1} = -theta_1 */
    acb_mat_clear(tau);
    acb_mat_init(tau, 1, 1);
    acb_init(t1); acb_init(t2); acb_init(t3); acb_init(t4); acb_init(w);
    acb_set_d_d(acb_mat_entry(tau, 0, 0), 0.25, 1.5);
    acb_set_d_d(z, 0.125, 0.25);
    acb_theta_all(th, z, tau, 0, prec);
    acb_modular_theta(t1, t2, t3, t4, z, acb_mat_entry(tau, 0, 0), prec);
    acb_add(w, th + 3, t1, prec);
    flint_printf("theta_11 + theta_1 = "); acb_printn(w, 5, 0); flint_printf("\n");

    acb_clear(t1); acb_clear(t2); acb_clear(t3); acb_clear(t4); acb_clear(w);
    acb_mat_clear(tau);
    _acb_vec_clear(z, 2);
    _acb_vec_clear(th, 16 * 3);
    return 0;
}
```
```output
theta_0000 = 1.08984692468251 + 0e-36*I
theta_0001 = 1.08302269779695 + 0e-36*I
theta_0010 = 0.917636847761688 + 0e-36*I
theta_0011 = 0.909521428594988 + 0e-36*I
theta_0100 = 0.415757702796908 + 0e-37*I
theta_0101 = odd, contains zero: 1
theta_0110 = 0.415757702796908 + 0e-37*I
theta_0111 = odd, contains zero: 1
theta_1000 = 0.913579138133896 + 0e-37*I
theta_1001 = 0.913579138133896 + 0e-37*I
theta_1010 = odd, contains zero: 1
theta_1011 = odd, contains zero: 1
theta_1100 = 0.267577401846795 + 0e-37*I
theta_1101 = odd, contains zero: 1
theta_1110 = odd, contains zero: 1
theta_1111 = 0e-37 - 0.267577401846795*I
jet of theta_1010: 0e-37 + 0e-37*I, -2.84869460392 + 0e-36*I, 0e-37 - 0.0214389410501*I
theta_11 + theta_1 = [+/- 6.09e-38] + [+/- 8.80e-38]*I
```

The three entries of the jet are $\theta$, $\partial\theta/\partial z_0$
and $\partial\theta/\partial z_1$, in the order of
`acb_theta_jet_tuples`: $(0,0)$, $(1,0)$, $(0,1)$.
