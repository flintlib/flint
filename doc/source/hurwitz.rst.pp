.. _algorithms_hurwitz:

Algorithms for the Hurwitz zeta function
===============================================================================

Euler-Maclaurin summation
-------------------------------------------------------------------------------

The Euler-Maclaurin formula allows evaluating the Hurwitz zeta function and
its derivatives for general complex input. The algorithm is described
in [Joh2013]_.

Parameter Taylor series
-------------------------------------------------------------------------------

To evaluate `\zeta(s,a)` for several nearby parameter values, the following
Taylor expansion is useful:

.. math::

    \zeta(s,a+x) = \sum_{k=0}^{\infty} (-x)^k \frac{(s)_k}{k!} \zeta(s+k,a)

We assume that `a \ge 1` is real and that `\sigma = \operatorname{re}(s)`
with `K + \sigma > 1`. The tail is bounded by

.. math::

    \sum_{k=K}^{\infty} |x|^k \frac{|(s)_k|}{k!} \zeta(\sigma+k,a) \le
    \sum_{k=K}^{\infty}
        |x|^k \frac{|(s)_k|}{k!} \left[
            \frac{1}{a^{\sigma+k}} + \frac{1}{(\sigma+k-1) a^{\sigma+k-1}} \right].

Denote the term on the right by `T(k)`. Then

.. math::

    \left|\frac{T(k+1)}{T(k)}\right| =
            \frac{|x|}{a}
            \frac{(k+\sigma-1)}{(k+\sigma)}
            \frac{(k+\sigma+a)}{(k+\sigma+a-1)}
            \frac{|k+s|}{(k+1)}
        \le
            \frac{|x|}{a}
            \left(1 + \frac{1}{K+\sigma+a-1}\right)
            \left(1 + \frac{|s-1|}{K+1}\right) = C

and if `C < 1`,

.. math::

    \sum_{k=K}^{\infty} T(k) \le \frac{T(K)}{1-C}.

