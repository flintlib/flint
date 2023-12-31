.. _algorithms_gamma:

Algorithms for the gamma function
===============================================================================

The Stirling series
-------------------------------------------------------------------------------

In general, the gamma function is computed via the Stirling series

.. math::

    \log \Gamma(z) = \left(z-\frac{1}{2}\right)\log z - z +
          \frac{\ln {2 \pi}}{2}
            + \sum_{k=1}^{n-1}  \frac{B_{2k}}{2k(2k-1)z^{2k-1}}
          + R(n,z)

where ([Olv1997]_ pp. 293-295) the remainder term is exactly

.. math::

    R_n(z) = \int_0^{\infty} \frac{B_{2n} - {\tilde B}_{2n}(x)}{2n(x+z)^{2n}} dx.

To evaluate the gamma function of a power series argument, we substitute
`z \to z + t \in \mathbb{C}[[t]]`.

Using the bound for `|x+z|` given by [Olv1997]_ and the fact
that the numerator of the integrand is bounded in
absolute value by `2 |B_{2n}|`, the remainder can be shown
to satisfy the bound

.. math::

    |[t^k] R_n(z+t)| \le 2 |B_{2n}|
        \frac{\Gamma(2n+k-1)}{\Gamma(k+1) \Gamma(2n+1)}
        \; |z| \; \left(\frac{b}{|z|}\right)^{2n+k}

where `b = 1/\cos(\operatorname{arg}(z)/2)`.
Note that by trigonometric identities, assuming that `z = x+yi`, we
have `b = \sqrt{1 + u^2}` where

.. math::

    u = \frac{y}{\sqrt{x^2 + y^2} + x} = \frac{\sqrt{x^2 + y^2} - x}{y}.

To use the Stirling series at `p`-bit precision,
we select parameters `r`, `n` such that the
remainder `R(n,z)` approximately is bounded by `2^{-p}`.
If `|z|` is too small for the Stirling series
to give sufficient accuracy directly, we first translate to `z + r`
using the formula `\Gamma(z) = \Gamma(z+r) / 
(z (z+1) (z+2) \cdots (z+r-1))`.

To obtain a remainder smaller than `2^{-p}`, we must choose an `r` such
that, in the real case, `z + r > \beta p`, where
`\beta > \log(2) / (2 \pi) \approx 0.11`.
In practice, a slightly larger factor `\beta \approx 0.2` more closely
balances `n` and `r`. A much larger `\beta` (e.g. `\beta = 1`) could be
used to reduce the number of Bernoulli numbers that have to be
precomputed, at the expense of slower repeated evaluation.

Rational arguments
-------------------------------------------------------------------------------

We use efficient methods to compute `y = \Gamma(p/q)` where `q` is
one of `1, 2, 3, 4, 6` and `p` is a small integer.

The cases `\Gamma(1) = 1` and `\Gamma(1/2) = \sqrt \pi` are trivial.
We reduce all remaining cases to `\Gamma(1/3)` or `\Gamma(1/4)`
using the following relations:

.. math::

    \Gamma(2/3) = \frac{2 \pi}{3^{1/2} \Gamma(1/3)}, \quad \quad
    \Gamma(3/4) = \frac{2^{1/2} \pi}{\Gamma(1/4)},

.. math::

    \Gamma(1/6) = \frac{\Gamma(1/3)^2}{(\pi/3)^{1/2} 2^{1/3}}, \quad \quad
    \Gamma(5/6) = \frac{2 \pi (\pi/3)^{1/2} 2^{1/3}}{\Gamma(1/3)^2}.

We compute `\Gamma(1/3)` and `\Gamma(1/4)` rapidly to high precision using

.. math::

    \Gamma(1/3) = \left( \frac{12 \pi^4}{\sqrt{10}}
        \sum_{k=0}^{\infty}
        \frac{(6k)!(-1)^k}{(k!)^3 (3k)! 3^k 160^{3k}} \right)^{1/6}, \quad \quad
    \Gamma(1/4) = \sqrt{\frac{(2\pi)^{3/2}}{\operatorname{agm}(1, \sqrt 2)}}.

An alternative formula which could be used for `\Gamma(1/3)` is

.. math::

    \Gamma(1/3) = \frac{2^{4/9} \pi^{2/3}}{3^{1/12} \left( \operatorname{agm}\left(1,\frac{1}{2} \sqrt{2+\sqrt{3}}\right)\right)^{1/3}},

but this appears to be slightly slower in practice.

