.. _algorithms_constants:

Algorithms for mathematical constants
===============================================================================

Most mathematical constants are evaluated using the generic hypergeometric
summation code.

Pi
-------------------------------------------------------------------------------

`\pi` is computed using the Chudnovsky series

    .. math::

        \frac{1}{\pi} = 12 \sum^\infty_{k=0}
        \frac{(-1)^k (6k)! (13591409 + 545140134k)}{(3k)!(k!)^3 640320^{3k + 3/2}}

which is hypergeometric and adds roughly 14 digits per term. Methods based on the
arithmetic-geometric mean seem to be slower by a factor three in practice.

A small trick
is to compute `1/\sqrt{640320}` instead of `\sqrt{640320}` at the end.

Logarithms of integers
-------------------------------------------------------------------------------

We use the formulas

.. math::

    \log(2) = \frac{3}{4} \sum_{k=0}^{\infty} \frac{(-1)^k (k!)^2}{2^k (2k+1)!}

.. math::

    \log(10) = 46 \operatorname{atanh}(1/31) + 34 \operatorname{atanh}(1/49) + 20 \operatorname{atanh}(1/161)


Euler's constant
-------------------------------------------------------------------------------

Euler's constant `\gamma` is computed using
the Brent-McMillan formula ([BM1980]_,  [MPFR2012]_)

.. math::

    \gamma = \frac{S_0(2n) - K_0(2n)}{I_0(2n)} - \log(n)

in which `n` is a free parameter and

.. math::

    S_0(x) = \sum_{k=0}^{\infty} \frac{H_k}{(k!)^2} \left(\frac{x}{2}\right)^{2k}, \quad
    I_0(x) = \sum_{k=0}^{\infty} \frac{1}{(k!)^2} \left(\frac{x}{2}\right)^{2k}

.. math::

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

.. math::

    C = \frac{1}{768} \sum_{k=1}^{\infty} \frac{(-4096)^k P(k)}
        {k^3 (2k-1)(3k-1)(3k-2)(6k-1)(6k-5) {5k \choose k} {10k \choose 5k} {12k \choose 6k}}

where

.. math::

    \begin{matrix}
        P(k) & = -43203456k^6 + 92809152k^5 - 76613904k^4 \\
             & + 30494304k^3 - 6004944k^2 + 536620^k - 17325,
    \end{matrix}

discovered by Zuniga [Zun2023]_.
It was previously computed using a series given in [PP2010]_.

Apery's constant
-------------------------------------------------------------------------------

Apery's constant `\zeta(3)` is computed using the hypergeometric series

.. math::

    \zeta(3) = \frac{1}{48} \sum_{k=1}^{\infty} \frac{(-1)^{k-1} P(k)}{k^5 (2k-1)^3(3k-1)(3k-2)(4k-1)(4k-3)(6k-1)(6k-5){5k \choose k}{5k \choose 2k}{9k \choose 4k}{10k \choose 5k}{12k \choose 6k}}

where

.. math::

    \begin{matrix}
        P(k) & = 1565994397644288k^{11} - 6719460725627136k^{10} + 12632254526031264k^9 \\
             & - 13684352515879536k^8 + 9451223531851808k^7 - 4348596587040104k^6 \\
             & + 1352700034136826k^5 - 282805786014979k^4 + 38721705264979k^3 \\
             & - 3292502315430k^2 + 156286859400k - 3143448000,
    \end{matrix}

discovered by Zuniga [Zun2023]_.

Khinchin's constant
-------------------------------------------------------------------------------

Khinchin's constant `K_0` is computed using the formula

.. math::

    \log K_0 = \frac{1}{\log 2} \left[
    \sum_{k=2}^{N-1} \log \left(\frac{k-1}{k} \right) \log \left(\frac{k+1}{k} \right)
    + \sum_{n=1}^\infty 
    \frac {\zeta (2n,N)}{n} \sum_{k=1}^{2n-1} \frac{(-1)^{k+1}}{k}
    \right]

where `N \ge 2` is a free parameter that can be used for tuning [BBC1997]_.
If the infinite series is truncated after `n = M`, the remainder
is smaller in absolute value than

.. math::

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

Reciprocal Fibonacci constant
-------------------------------------------------------------------------------

We use Gosper's series ([Gos1974]_, corrected in [Arn2012]_)

.. math::

    \sum_{n=1}^{\infty} \frac{1}{F_n} = \sum_{n=0}^{\infty}
        \frac{(-1)^{n(n-1)/2} (F_{4n+3} + (-1)^n F_{2n+2})}{F_{2n+1} F_{2n+2} L_1 L_3 \cdots L_{2n+1}}

where `L_n = 2F_{n-1} + F_n` denotes a Lucas number.
The truncation error after `N \ge 1` terms is bounded by `(1 / \phi)^{N^2}`.
The series is not of hypergeometric type, but we can evaluate it
in quasilinar time using binary splitting; factoring out a
multiplicative recurrence for `L_1 L_3 \cdots` allows computing the series
as a product of `O(\sqrt{p})` matrices with `O(\sqrt{p})`-bit entries.
