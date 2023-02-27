.. _algorithms_constants:

Algorithms for mathematical constants
===============================================================================

Most mathematical constants are evaluated using the generic hypergeometric
summation code.

Pi
-------------------------------------------------------------------------------

`\pi` is computed using the Chudnovsky series

    .. math ::

        \frac{1}{\pi} = 12 \sum^\infty_{k=0}
        \frac{(-1)^k (6k)! (13591409 + 545140134k)}{(3k)!(k!)^3 640320^{3k + 3/2}}

which is hypergeometric and adds roughly 14 digits per term. Methods based on the
arithmetic-geometric mean seem to be slower by a factor three in practice.

A small trick
is to compute `1/\sqrt{640320}` instead of `\sqrt{640320}` at the end.

Logarithms of integers
-------------------------------------------------------------------------------

We use the formulas

.. math ::

    \log(2) = \frac{3}{4} \sum_{k=0}^{\infty} \frac{(-1)^k (k!)^2}{2^k (2k+1)!}

.. math ::

    \log(10) = 46 \operatorname{atanh}(1/31) + 34 \operatorname{atanh}(1/49) + 20 \operatorname{atanh}(1/161)


Euler's constant
-------------------------------------------------------------------------------

Euler's constant `\gamma` is computed using
the Brent-McMillan formula ([BM1980]_,  [MPFR2012]_)

.. math ::

    \gamma = \frac{S_0(2n) - K_0(2n)}{I_0(2n)} - \log(n)

in which `n` is a free parameter and

.. math ::

    S_0(x) = \sum_{k=0}^{\infty} \frac{H_k}{(k!)^2} \left(\frac{x}{2}\right)^{2k}, \quad
    I_0(x) = \sum_{k=0}^{\infty} \frac{1}{(k!)^2} \left(\frac{x}{2}\right)^{2k}

.. math ::

    2x I_0(x) K_0(x) \sim \sum_{k=0}^{\infty} \frac{[(2k)!]^3}{(k!)^4 8^{2k} x^{2k}}.

All series are evaluated using binary splitting.
The first two series are evaluated simultaneously, with the summation
taken up to `k = N - 1` inclusive where `N \ge \alpha n + 1` and
`\alpha \approx 4.9706257595442318644`
satisfies `\alpha (\log \alpha - 1) = 3`. The third series is taken
up to `k = 2n-1` inclusive. With these parameters, it is shown in
[BJ2013]_ that the error is bounded by `24e^{-8n}`.

Catalan's constant
-------------------------------------------------------------------------------

Catalan's constant is computed using the hypergeometric series

.. math ::

    C = \frac{1}{64} \sum_{k=1}^{\infty} \frac{256^k (580k^2-184k+15)}{k^3(2k-1){6k\choose 3k}{6k\choose 4k}{4k\choose 2k}}

given in [PP2010]_.

Khinchin's constant
-------------------------------------------------------------------------------

Khinchin's constant `K_0` is computed using the formula

.. math ::

    \log K_0 = \frac{1}{\log 2} \left[
    \sum_{k=2}^{N-1} \log \left(\frac{k-1}{k} \right) \log \left(\frac{k+1}{k} \right)
    + \sum_{n=1}^\infty 
    \frac {\zeta (2n,N)}{n} \sum_{k=1}^{2n-1} \frac{(-1)^{k+1}}{k}
    \right]

where `N \ge 2` is a free parameter that can be used for tuning [BBC1997]_.
If the infinite series is truncated after `n = M`, the remainder
is smaller in absolute value than

.. math ::

    \sum_{n=M+1}^{\infty} \zeta(2n, N) = 
    \sum_{n=M+1}^{\infty} \sum_{k=0}^{\infty} (k+N)^{-2n} \le
    \sum_{n=M+1}^{\infty} \left( N^{-2n} + \int_0^{\infty} (t+N)^{-2n} dt \right)

    = \sum_{n=M+1}^{\infty} \frac{1}{N^{2n}} \left(1 + \frac{N}{2n-1}\right)
    \le \sum_{n=M+1}^{\infty} \frac{N+1}{N^{2n}} = \frac{1}{N^{2M} (N-1)}
    \le \frac{1}{N^{2M}}.

Thus, for an error of at most `2^{-p}` in the series,
it is sufficient to choose `M \ge p / (2 \log_2 N)`.

Glaisher's constant
-------------------------------------------------------------------------------

Glaisher's constant `A = \exp(1/12 - \zeta'(-1))` is computed directly
from this formula. We don't use the reflection formula for the zeta function,
as the arithmetic in Euler-Maclaurin summation is faster at `s = -1`
than at `s = 2`.

Apery's constant
-------------------------------------------------------------------------------

Apery's constant `\zeta(3)` is computed using the hypergeometric series

.. math ::

    \zeta(3) = \frac{1}{64} \sum_{k=0}^\infty
        (-1)^k (205k^2 + 250k + 77) \frac{(k!)^{10}}{[(2k+1)!]^5}.

