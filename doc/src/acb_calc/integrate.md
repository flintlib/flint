---
title: Rigorous numerical integration
functions: acb_calc_integrate acb_calc_integrate_opt_init
---

`acb_calc_integrate(res, f, param, a, b, rel_goal, abs_tol, options,
prec)` computes a ball that contains $\int_a^b f(t)\, dt$, where the
path is the segment from $a$ to $b$ in the complex plane. The enclosure
is rigorous whatever the integrand, provided that the evaluator `f`
returns correct enclosures and correctly reports non-holomorphy; the
accuracy parameters only determine how much work is spent to make it
narrow. The algorithm is that of [@Joh2018a]: adaptive bisection, with
Gauss–Legendre quadrature and error bounds from complex analysis on
segments where $f$ is holomorphic.

::: summary algorithms/ball-arithmetic/integration :::

## The algorithm {#algorithm}

The algorithm keeps a stack (or, with `use_heap`, a priority queue) of
segments, each with a crude enclosure of its integral, starting with
the whole path. For the segment $[u, v]$ at the top:

1. The *direct enclosure* $(v - u) f(B)$, with $B$ the ball containing
   the segment, evaluated with `order` 0, is accepted if its radius is
   below the current tolerance, or if $u$ and $v$ can no longer be
   distinguished at the working precision.
2. Otherwise, if the direct enclosure is finite,
   `acb_calc_integrate_gl_auto_deg` tries a Gauss–Legendre rule. It
   evaluates $f$ with `order` 1 (that is, with a holomorphy check) on
   rectangles enclosing ellipses around the segment, of increasing size,
   to obtain bounds $M$ for $|f|$, and chooses the smallest degree $n
   \le$ `deg_limit` for which the error bound
   $$\frac{64 M}{15 (\rho - 1) \rho^{2n - 1}}$$
   is below the tolerance, where $\rho$ is the sum of the semi-axes of
   the ellipse relative to the half-length of the segment. If the
   integral of $f$ is real on the segment (because the direct enclosure
   is real), the imaginary part of the result is set to zero.
3. Otherwise the segment is bisected and the direct enclosures of the
   two halves are computed; the half with the larger error is processed
   first.

The tolerance on a segment is $\max(\varepsilon_{\mathrm{abs}}, M_{\mathrm{low}}
2^{-r})$, where $M_{\mathrm{low}}$ is the largest lower bound for the
absolute value of an integral over any segment seen so far. When
`eval_limit` evaluations or `depth_limit` queued segments are reached,
the remaining segments are accepted with their direct enclosures, which
keeps the result valid but possibly wide, and
`ARB_CALC_NO_CONVERGENCE` is returned.

The error bound in step 2 means that the number of evaluations is
proportional to the number of correct bits divided by $\log \rho$: the
convergence is exponential, at a rate given by the distance of the
nearest singularity from the segment. Near a singularity or a
discontinuity, the segments become shorter (a geometric subdivision),
so that, for example, integrands with algebraic or logarithmic
endpoint singularities still converge with a number of evaluations
roughly proportional to the square of the number of bits.

## Writing integrands {#integrands}

The integrand is called with `order` 0 for the direct enclosures and
the quadrature nodes, and with `order` 1 for the ellipse bounds. With
`order` 0 it may be any function, also discontinuous; with `order` 1 it
must return a non-finite value if $f$ is not holomorphic on the input
ball. Poles are detected automatically by ball arithmetic (division by
a ball containing zero gives an infinite result), but branch cuts and
non-analytic functions such as $|z|$ or $\lfloor z \rfloor$ need the
`analytic` variants `acb_sqrt_analytic`, `acb_log_analytic`,
`acb_pow_analytic`, `acb_real_abs`, `acb_real_floor` and so on, called
with the flag `order != 0`. An integrand that ignores this rule can
make the result wrong, because the quadrature error bound then assumes
analyticity where there is none. See `acb_calc_func_t`.

## Choosing the parameters {#parameters}

- `prec` is the working precision. The algorithm does not increase it;
  if the result is not accurate enough because of cancellation in the
  integrand, the caller must repeat the computation at a higher
  precision.
- `rel_goal` is the relative accuracy goal in bits and `abs_tol` the
  absolute tolerance. The standard choice is `rel_goal` = `prec` and
  `abs_tol` $= 2^{-\mathit{prec}}$. For an integral of very small or very
  large magnitude $M$, set `abs_tol` to about $M 2^{-\mathit{rel\_goal}}$
  if $M$ can be estimated, or to zero.
- `options` can be `NULL`. Otherwise initialise it with
  `acb_calc_integrate_opt_init` and set the fields to change:
  `deg_limit` (raise it for oscillatory entire functions), `eval_limit`
  (lower it to give up sooner), `depth_limit`, `use_heap` and `verbose`.

## Examples {#examples}

The program integrates three functions at 64-bit precision: the
oscillatory entire function $\sin(x + e^x)$ on $[0, 8]$ (which has
about 950 sign changes on this interval), $\sqrt{1 - x^2}$ on $[-1, 1]$
(with algebraic branch points at both endpoints; twice the integral is
$\pi$), and the piecewise constant function $\lfloor x \rfloor$ on $[0,
10.5]$, whose integral is $0 + 1 + \cdots + 9 + 10 \cdot 0.5 = 50$.

```c example name=acb_calc_integrate_basic
#include <flint/acb.h>
#include <flint/acb_calc.h>

/* f(z) = sin(z + exp(z)), entire */
static int
f_sin_exp(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* derivatives not implemented */
    acb_exp(res, z, prec);
    acb_add(res, res, z, prec);
    acb_sin(res, res, prec);
    return 0;
}

/* f(z) = sqrt(1 - z^2), with a branch cut: must detect non-analyticity */
static int
f_circle(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();
    acb_one(res);
    acb_submul(res, z, z, prec);
    acb_sqrt_analytic(res, res, order != 0, prec);
    return 0;
}

/* f(z) = floor(z), piecewise constant */
static int
f_floor(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();
    acb_real_floor(res, z, order != 0, prec);
    return 0;
}

int main(void)
{
    acb_t s, a, b;
    mag_t tol;
    slong prec = 64;

    acb_init(s);
    acb_init(a);
    acb_init(b);
    mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -prec);

    acb_set_si(a, 0);
    acb_set_si(b, 8);
    acb_calc_integrate(s, f_sin_exp, NULL, a, b, prec, tol, NULL, prec);
    flint_printf("int_0^8 sin(x + exp(x)) dx = ");
    acb_printn(s, 15, 0);
    flint_printf("\n");

    acb_set_si(a, -1);
    acb_set_si(b, 1);
    acb_calc_integrate(s, f_circle, NULL, a, b, prec, tol, NULL, prec);
    acb_mul_2exp_si(s, s, 1);
    flint_printf("2 int_-1^1 sqrt(1 - x^2) dx = ");
    acb_printn(s, 15, 0);
    flint_printf("\n");

    acb_set_si(a, 0);
    acb_set_d(b, 10.5);
    acb_calc_integrate(s, f_floor, NULL, a, b, prec, tol, NULL, prec);
    flint_printf("int_0^10.5 floor(x) dx = ");
    acb_printn(s, 15, 0);
    flint_printf("\n");

    acb_clear(s);
    acb_clear(a);
    acb_clear(b);
    mag_clear(tol);
    return 0;
}
```
```output
int_0^8 sin(x + exp(x)) dx = [0.34740017265725 +/- 3.34e-15]
2 int_-1^1 sqrt(1 - x^2) dx = [3.14159265358979 +/- 3.25e-15] + [+/- 1.54e-19]*I
int_0^10.5 floor(x) dx = [50.0000000000000 +/- 1.10e-15]
```

The second result has a small imaginary radius: near the endpoints the
direct enclosures of $\sqrt{1 - z^2}$ on balls that cross $\pm 1$ are
not real. The printed radii are dominated by the conversion to 15
digits.
