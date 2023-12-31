.. _algorithms_polylogarithms:

Algorithms for polylogarithms
===============================================================================

The polylogarithm is defined for `s, z \in \mathbb{C}` with `|z| < 1` by

.. math::

    \operatorname{Li}_s(z) = \sum_{k=1}^{\infty} \frac{z^k}{k^s}

and for `|z| \ge 1` by analytic continuation, except for the singular
point `z = 1`.

Computation for small z
-------------------------------------------------------------------------------

The power sum converges rapidly when `|z| \ll 1`.
To compute the series expansion with respect to `s`, we substitute
`s \to s + x \in \mathbb{C}[[x]]` and obtain

.. math::

    \operatorname{Li}_{s+x}(z) = \sum_{d=0}^{\infty} x^d
        \frac{(-1)^d}{d!} \sum_{k=1}^{\infty} T(k)

where

.. math::

        T(k) = \frac{z^k \log^d(k)}{k^s}.

The remainder term `\left| \sum_{k=N}^{\infty} T(k) \right|` is bounded
via the following strategy, implemented in :func:`mag_polylog_tail`.

Denote the terms by `T(k)`. We pick a nonincreasing function `U(k)` such that

.. math::

    \frac{T(k+1)}{T(k)} = z \left(\frac{k}{k+1}\right)^s
        \left( \frac{\log(k+1)}{\log(k)} \right)^d \le U(k).

Then, as soon as `U(N) < 1`,

.. math::

    \sum_{k=N}^{\infty} T(k)
        \le T(N) \sum_{k=0}^{\infty} U(N)^k = \frac{T(N)}{1 - U(N)}.

In particular, we take

.. math::

    U(k) = z \; B(k, \max(0, -s)) \; B(k \log(k), d)

where `B(m,n) = (1 + 1/m)^n`. This follows from the bounds

.. math::

    \left(\frac{k}{k+1}\right)^{s}
    \le \begin{cases}
        1                    & \text{if }         s \ge 0 \\
        (1 + 1/k)^{-s}  & \text{if }         s < 0.
        \end{cases}

and

.. math::

    \left( \frac{\log(k+1)}{\log(k)} \right)^d \le
        \left(1 + \frac{1}{k \log(k)}\right)^d.

Expansion for general z
-------------------------------------------------------------------------------

For general complex `s, z`, we write the polylogarithm as a sum of
two Hurwitz zeta functions

.. math::

    \operatorname{Li}_s(z) = \frac{\Gamma(v)}{(2\pi)^v}
        \left[
            i^v
            \zeta \left(v, \frac{1}{2} + \frac{\log(-z)}{2\pi i}\right)
            + i^{-v}
            \zeta \left(v, \frac{1}{2} - \frac{\log(-z)}{2\pi i}\right)
        \right]

in which `s = 1-v`.
With the principal branch of `\log(-z)`, we obtain the conventional
analytic continuation of the polylogarithm with a branch
cut on `z \in (1,+\infty)`.

To compute the series expansion with respect to `v`, we substitute
`v \to v + x \in \mathbb{C}[[x]]` in this formula
(at the end of the computation, we map `x \to -x` to
obtain the power series for `\operatorname{Li}_{s+x}(z)`).
The right hand side becomes

.. math::

    \Gamma(v+x) [E_1 Z_1 + E_2 Z_2]

where `E_1 = (i/(2 \pi))^{v+x}`, `Z_1 = \zeta(v+x,\ldots)`,
`E_2 = (1/(2 \pi i))^{v+x}`, `Z_2 = \zeta(v+x,\ldots)`.

When `v = 1`, the `Z_1` and `Z_2` terms become Laurent series with
a leading `1/x` term. In this case,
we compute the deflated series `\tilde Z_1, \tilde Z_2 = \zeta(x,\ldots) - 1/x`.
Then

.. math::

    E_1 Z_1 + E_2 Z_2 = (E_1 + E_2)/x + E_1 \tilde Z_1 + E_2 \tilde Z_2.

Note that `(E_1 + E_2) / x` is a power series, since the constant term in
`E_1 + E_2` is zero when `v = 1`. So we simply compute one extra derivative
of both `E_1` and `E_2`, and shift them one step.
When `v = 0, -1, -2, \ldots`, the `\Gamma(v+x)` prefactor has a pole.
In this case, we proceed analogously and formally multiply
`x \, \Gamma(v+x)` with `[E_1 Z_1 + E_2 Z_2] / x`.

Note that the formal cancellation only works when the order `s` (or `v`)
is an exact integer: it is not currently possible to use this method when
`s` is a small ball containing any of `0, 1, 2, \ldots` (then the
result becomes indeterminate).

The Hurwitz zeta method becomes inefficient when `|z| \to 0` (it
gives an indeterminate
result when `z = 0`). This is not a problem since we just use the defining series
for the polylogarithm in that region.
It also becomes inefficient when `|z| \to \infty`, for which an asymptotic
expansion would better.

