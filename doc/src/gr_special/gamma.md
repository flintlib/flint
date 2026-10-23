---
title: The gamma function and factorials
functions: gr_gamma gr_fac
---

`gr_gamma` and the factorial functions show how one special function
interface serves very different rings. The same call computes
$\Gamma(x)$ as a rigorous real or complex ball, as a power series when
$x$ is a series, or reports a pole; `gr_fac_ui` computes $n!$ exactly in
$\Z$ and directly modulo $p$ in $\Z/p\Z$, without forming the integer.

## Rings and algorithms {#rings}

* **Real and complex balls.** `gr_gamma` calls `arb_gamma` or
  `acb_gamma`, which use the Stirling series with argument reduction or,
  for small rational and integer arguments, exact or rectangular
  splitting formulas [@Joh2023g]; see
  [](algorithms/special-functions/gamma). At the poles
  $0, -1, -2, \ldots$ it returns `GR_DOMAIN` when the argument is an exact
  integer, and `GR_UNABLE` when a ball merely contains a pole.
  `gr_rgamma` computes the entire function $1/\Gamma(x)$, which is zero
  at the poles.
* **Power series.** In `gr_series` rings over real or complex balls,
  $\Gamma(f)$ for a power series $f$ is computed with
  `arb_poly_gamma_series` or `acb_poly_gamma_series` (the Taylor
  expansion of $\Gamma$ at the constant term of $f$, composed with $f$),
  and the result carries the error term of $f$.
* **Floating point.** `nfloat`, `arf`, `dfloat` and `decimal` contexts
  compute approximations without automatic error tracking.
* **Exact rings.** $\Gamma(x)$ is not available in $\Z$ or $\Q$ as a
  method, but `gr_fac`, `gr_fac_ui` and `gr_fac_fmpz` work in every ring
  through the generic implementation `gr_generic_fac_ui`:
  `fmpz_fac_ui` in $\Z$, a table for small $n$, products of several
  consecutive integers per multiplication in rings of finite
  characteristic (with a fast algorithm for `nmod`), ball arithmetic in
  rings with a real precision, and otherwise the exact integer converted
  to the ring. `gr_fac` of a non-integer computes $\Gamma(x + 1)$.

## Example {#example}

```c example name=gr_gamma_rings
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_series.h>
#include <flint/gr_special.h>

void show_gamma(const char * s, gr_ctx_t ctx)
{
    gr_ptr x, y;
    int status;

    x = gr_heap_init(ctx);
    y = gr_heap_init(ctx);
    GR_MUST_SUCCEED(gr_set_str(x, s, ctx));
    status = gr_gamma(y, x, ctx);
    flint_printf("gamma(%s) in %{gr_ctx}: ", s, ctx);
    if (status == GR_SUCCESS)
        flint_printf("%{gr}\n", y, ctx);
    else
        flint_printf("status %d\n", status);
    gr_heap_clear(x, ctx);
    gr_heap_clear(y, ctx);
}

int main(void)
{
    gr_ctx_t ZZ, ZZp, RR, CC, RRx;
    gr_ptr f;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_nmod(ZZp, 101);
    gr_ctx_init_real_arb(RR, 64);
    gr_ctx_init_complex_acb(CC, 64);
    gr_series_ctx_init(RRx, RR, 4);

    show_gamma("1/3", RR);
    show_gamma("1+i", CC);
    show_gamma("-2", RR);
    show_gamma("1+x", RRx);

    f = gr_heap_init(ZZ);
    GR_MUST_SUCCEED(gr_fac_ui(f, 30, ZZ));
    flint_printf("30! = %{gr}\n", f, ZZ);
    gr_heap_clear(f, ZZ);
    f = gr_heap_init(ZZp);
    GR_MUST_SUCCEED(gr_fac_ui(f, 30, ZZp));
    flint_printf("30! mod 101 = %{gr}\n", f, ZZp);
    GR_MUST_SUCCEED(gr_fac_ui(f, 100, ZZp));
    flint_printf("100! mod 101 = %{gr}\n", f, ZZp);
    gr_heap_clear(f, ZZp);

    gr_ctx_clear(RRx);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(ZZp);
    gr_ctx_clear(RR);
    gr_ctx_clear(CC);
    return 0;
}
```
```output
gamma(1/3) in Real numbers (arb, prec = 64): [2.678938534707747634 +/- 9.48e-19]
gamma(1+i) in Complex numbers (acb, prec = 64): ([0.4980156681183560427 +/- 6.61e-20] + [-0.1549498283018106851 +/- 4.87e-20]*I)
gamma(-2) in Real numbers (arb, prec = 64): status 1
gamma(1+x) in Power series over Real numbers (arb, prec = 64) with precision O(x^4): 1 + [-0.5772156649015328606 +/- 9.24e-20]*x + [0.9890559953279725554 +/- 6.79e-20]*x^2 + [-0.9074790760808862890 +/- 6.24e-20]*x^3 + O(x^4)
30! = 265252859812191058636308480000000
30! mod 101 = 74
100! mod 101 = 100
```

The coefficients of $\Gamma(1 + x)$ are $1$, $-\gamma$,
$\gamma^2/2 + \pi^2/12$, ..., and $100! \equiv -1 \pmod{101}$ by
Wilson's theorem.
