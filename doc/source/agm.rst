.. _algorithms_agm:

Algorithms for the arithmetic-geometric mean
===============================================================================

With complex variables, it is convenient to work with the univariate
function `M(z) = \operatorname{agm}(1,z)`. The general case is given by
`\operatorname{agm}(a,b) = a M(1,b/a)`.

Functional equation
------------------------------------------------------------------------------

If the real part of *z* initially is not completely nonnegative, we
apply the functional equation `M(z) = (z+1) M(u) / 2`
where `u = \sqrt{z} / (z+1)`.

Note that *u* has nonnegative real part, absent rounding error.
It is not a problem for correctness if rounding makes the interval
contain negative points, as this just inflates the final result.

For the derivative, the functional equation becomes
`M'(z) = [M(u) - (z-1) M'(u) / ((1+z) \sqrt{z})] / 2`.

AGM iteration
------------------------------------------------------------------------------

Once *z* is in the right half plane, we can apply the AGM iteration
(`2a_{n+1} = a_n + b_n, b_{n+1}^2 = a_n b_n`) directly.
The correct square root is given by `\sqrt{a} \sqrt{b}`,
which is computed as `\sqrt{ab}, i \sqrt{-ab}, -i \sqrt{-ab}, \sqrt{a} \sqrt{b}`
respectively if both *a* and *b* have positive real part, nonnegative
imaginary part, nonpositive imaginary part, or otherwise.

The iteration should be terminated when `a_n` and `b_n` are close enough.
For positive real variables, we can simply take lower and upper bounds
to get a correct enclosure at this point. For complex variables, it is shown
in [Dup2006]_, p. 87 that, for *z* with nonnegative real part,
`|M(z) - a_n| \le |a_n - b_n|`, giving a convenient error bound.

Rather than running the AGM iteration until `a_n` and `b_n` agree to
`p` bits, it is slightly more efficient to iterate until
they agree to about `p/10` bits and finish with a series expansion.
With `z = (a-b)/(a+b)`, we have

.. math::

    \operatorname{agm}(a,b) = \frac{(a+b) \pi}{4 K(z^2)},

valid at least when `|z| < 1` and `a, b` have nonnegative real part,
and

.. math::

    \frac{\pi}{4 K(z^2)} = \tfrac{1}{2} - \tfrac{1}{8} z^2 - \tfrac{5}{128} z^4 - \tfrac{11}{512} z^6 - \tfrac{469}{32768} z^8 + \ldots

where the tail is bounded by `\sum_{k=10}^{\infty} |z|^k/64`.

First derivative
------------------------------------------------------------------------------

Assuming that *z* is exact and that `|\arg(z)| \le 3 \pi / 4`,
we compute `(M(z), M'(z))` simultaneously using a finite difference.

The basic inequality we need is `|M(z)| \le \max(1, |z|)`, which is
an immediate consequence of the AGM iteration.

By Cauchy's integral formula, `|M^{(k)}(z) / k!| \le C D^k` where
`C = \max(1, |z| + r)` and `D = 1/r`, for any `0 < r < |z|` (we
choose *r* to be of the order `|z| / 4`). Taylor expansion now gives

.. math::

    \left|\frac{M(z+h) - M(z)}{h} - M'(z)\right| \le \frac{C D^2 h}{1 - D h}

    \left|\frac{M(z+h) - M(z-h)}{2h} - M'(z)\right| \le \frac{C D^3 h^2}{1 - D h}

    \left|\frac{M(z+h) + M(z-h)}{2} - M(z)\right| \le \frac{C D^2 h^2}{1 - D h}

assuming that *h* is chosen so that it satisfies `h D < 1`.

The forward finite difference would require two function evaluations
at doubled precision. We use the central difference as it
only requires 1.5 times the precision.

When *z* is not exact, we evaluate at the midpoint as above
and bound the propagated error using derivatives.
Again by Cauchy's integral formula, we have

.. math::

    |M'(z+\varepsilon)| \le \frac{\max(1, |z|+|\varepsilon|+r)}{r}

    |M''(z+\varepsilon)| \le \frac{2 \max(1, |z|+|\varepsilon|+r)}{r^2}

assuming that the circle centered on *z* with radius `|\varepsilon| + r`
does not cross the negative half axis. We choose *r* of order `|z| / 2`
and verify that all assumptions hold.

Higher derivatives
-------------------------------------------------------------------------------

The function `W(z) = 1 / M(z)` is D-finite. The coefficients of
`W(z+x) = \sum_{k=0}^{\infty} c_k x^k` satisfy

.. math::

    -2 z (z^2-1) c_2 = (3z^2-1) c_1 + z c_0,

.. math::

    -(k+2)(k+3) z (z^2-1) c_{k+3} = (k+2)^2 (3z^2-1) c_{k+2} + (3k(k+3)+7)z c_{k+1} + (k+1)^2 c_{k}

in general, and

.. math::

    -(k+2)^2 c_{k+2} = (3k(k+3)+7) c_{k+1} + (k+1)^2 c_{k}

when `z = 1`.

