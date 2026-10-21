---
title: Error functions of a real ball
functions: arb_hypgeom_erf arb_hypgeom_erfc arb_hypgeom_erfinv
---

`arb_hypgeom_erf` and `arb_hypgeom_erfc` compute the error function
$\operatorname{erf}(x) = \frac{2}{\sqrt\pi} \int_0^x e^{-t^2} dt$ and
its complement $\operatorname{erfc}(x) = 1 - \operatorname{erf}(x)$ for
a real ball $x$; `arb_hypgeom_erfinv` and `arb_hypgeom_erfcinv` compute
the inverse functions. Unlike most real functions of this module they
do not call the complex implementation: for real $x$ one of the two
convergent series has only positive terms, the asymptotic expansion
has a simple error bound, and at very high precision the incomplete
gamma function $\Gamma(1/2, x^2)$ can be computed with the bit-burst
algorithm. `acb_hypgeom_erf` calls `arb_hypgeom_erf` for real input.
The complex algorithms and the incomplete gamma functions in general
are described in [](algorithms/special-functions/error-functions).

::: summary algorithms/special-functions/error-functions :::

## The methods {#methods}

`arb_hypgeom_erf` chooses among the following.

- **Tiny or huge $|x|$.** For $|x| < 2^{-\mathit{prec}/16}$ the series
  below needs few terms. If $x^2/\log 2$ exceeds the precision (or
  $|x| > 2^{60}$), the function returns $\operatorname{sgn}(x) \pm
  \varepsilon$ with $\varepsilon = e^{-x^2}/(|x|\sqrt\pi)$, which
  bounds $|\operatorname{erfc}(|x|)|$.
- **Convergent series.** The default method is
  $$
  \operatorname{erf}(x) = \frac{2x}{\sqrt\pi} e^{-x^2} \sum_{k=0}^\infty
  \frac{(x^2)^k}{(3/2)_k},
  $$
  in which all terms are positive, so that no precision is lost to
  cancellation. The number of terms $N$ is estimated with the Lambert
  $W$ function, the tail is bounded by the tail of the exponential
  series ($(x^2)^k/(3/2)_k \le (x^2)^k/k!$), and the sum is computed by
  rectangular splitting with `_arb_hypgeom_gamma_lower_sum_rs_1`, so that
  most multiplications are by small integers.
- **Asymptotic expansion.** For $|x| > 2$, the expansion
  $$
  \operatorname{erfc}(x) = \frac{e^{-x^2}}{x\sqrt\pi}
  \left(\sum_{k=0}^{N-1} \frac{(1/2)_k}{(-x^2)^k} + R_N\right),
  \qquad |R_N| \le \frac{(1/2)_N}{x^{2N}} \le \frac{N!}{x^{2N}},
  $$
  is used when some $N$ makes the bound smaller than
  $2^{-\mathit{prec}}$; the sum is computed at the reduced precision
  that the factor $e^{-x^2}$ allows.
- **Bit-burst algorithm.** Above 30000 bits and for $150
  e^{-0.004\sqrt{\mathit{prec}}} < |x| < 0.6\sqrt{\mathit{prec}}$, the
  function computes $\operatorname{erfc}(|x|) = \Gamma(1/2,
  x^2)/\sqrt\pi$ with `arb_hypgeom_erf_bb`: $x^2$ is approached by
  numbers $z_0, z_1, \ldots$ with 64, 128, 256, ... bits, the value at
  $z_0$ comes from the asymptotic or the convergent series of the
  incomplete gamma function, and each step from $z_j$ to $z_{j+1}$
  sums a Taylor series by binary splitting. The cost is $O(M(n)
  \log^2 n)$ at $n$ bits, asymptotically less than the $O(n)$
  multiplications by small integers and $O(\sqrt n)$ full
  multiplications of rectangular splitting. See [](algorithms/special-functions/hypergeometric#bit-burst).

For inexact $x$, the series are evaluated at the midpoint and the
propagated error is bounded by $\frac{2}{\sqrt\pi} e^{-m^2} r$, where
$m$ is the lower bound of $|x|$ and $r$ the radius (and by 2); the
precision is reduced to the accuracy of $x$ plus the number of bits
that $e^{-x^2}$ gains.

`arb_hypgeom_erfc` computes $1 - \operatorname{erf}(x)$ for $x \le 1$.
For larger $x$, where the subtraction would cancel, it uses the
asymptotic expansion whenever it reaches the target accuracy (with a
single term for $x > 2^{\mathit{prec}/2+10}$), the bit-burst algorithm
above 30000 bits, and otherwise the convergent series for
$\operatorname{erf}$ at a working precision increased by $(x^2 + \log
x)/\log 2$ bits.

## Inverse functions {#inverse}

`arb_hypgeom_erfinv` solves $\operatorname{erf}(y) = x$ by Newton's
method with $1/\operatorname{erf}'(y) = \frac{\sqrt\pi}{2} e^{y^2}$.
The starting value comes from a rational approximation in double
precision [@Bla1976], from the first terms of the asymptotic expansion
of $\operatorname{erf}^{-1}$ near $\pm 1$ when $1 - |x| < 2^{-52}$, and
from the same expansion in ball arithmetic when $1 - |x| < 2^{-1000}$. Each iteration doubles the precision; when an iterate is an
interval Newton step contained in the previous enclosure, the existence
of the root in it is proved, and the loop stops once a proved enclosure
has the target accuracy. Near $|x| = 1$ the iteration works with $1 -
|x|$ and $\operatorname{erfc}$, which avoids cancellation.
`arb_hypgeom_erfcinv` passes $x$ itself in this form when $x \le 0.01$,
so that arguments such as $10^{-1000}$ lose no accuracy.

## Example {#example}

At 64 bits, $\operatorname{erf}(30)$ is $1$ up to a tiny error, and
only `arb_hypgeom_erfc` resolves the difference; `arb_hypgeom_erfcinv`
inverts it for a tiny argument:

```c example name=arb_hypgeom_erf_erfc
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>

int main(void)
{
    arb_t x, y, z;

    arb_init(x);
    arb_init(y);
    arb_init(z);

    arb_set_ui(x, 30);
    arb_hypgeom_erf(y, x, 64);
    flint_printf("erf(30)          = "); arb_printn(y, 20, 0); flint_printf("\n");
    arb_hypgeom_erfc(y, x, 64);
    flint_printf("erfc(30)         = "); arb_printn(y, 20, 0); flint_printf("\n");
    /* 1 - erf(30) cannot resolve the tiny difference at 64 bits */
    arb_hypgeom_erf(y, x, 64);
    arb_sub_ui(y, y, 1, 64);
    arb_neg(y, y);
    flint_printf("1 - erf(30)      = "); arb_printn(y, 20, 0); flint_printf("\n");

    /* inverse functions; erfcinv works directly with tiny arguments */
    arb_set_str(x, "0.5", 64);
    arb_hypgeom_erfinv(y, x, 64);
    flint_printf("erfinv(0.5)      = "); arb_printn(y, 20, 0); flint_printf("\n");
    arb_set_str(x, "1e-1000", 64);
    arb_hypgeom_erfcinv(y, x, 64);
    flint_printf("erfcinv(1e-1000) = "); arb_printn(y, 20, 0); flint_printf("\n");
    arb_hypgeom_erfc(z, y, 64);
    flint_printf("erfc of that     = "); arb_printn(z, 20, 0); flint_printf("\n");

    arb_clear(x);
    arb_clear(y);
    arb_clear(z);
    flint_cleanup();
    return 0;
}
```
```output
erf(30)          = [1.0000000000000000000 +/- 2.57e-393]
erfc(30)         = [2.5646562037561116000e-393 +/- 3.65e-413]
1 - erf(30)      = [+/- 2.57e-393]
erfinv(0.5)      = [0.4769362762044698734 +/- 6.83e-20]
erfcinv(1e-1000) = [47.93894640426451069 +/- 3.99e-18]
erfc of that     = [1.000000000000000e-1000 +/- 3.75e-1016]
```

The enclosure of $\operatorname{erf}(30)$ is the bound $1 \pm
e^{-900}/(30\sqrt\pi)$, which is correct but says nothing about
$\operatorname{erfc}(30)$ itself; computing $\operatorname{erfc}$
directly gives 64 significant bits. The string `"1e-1000"` is not
exactly representable, so the input to `arb_hypgeom_erfcinv` is a ball
of relative radius $2^{-64}$, and the last line recovers it within
that radius.
