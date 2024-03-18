.. _general_formulas:

General formulas and bounds
===============================================================================

This section collects some results from real and complex
analysis that are useful when deriving error bounds.
Beware of typos.

Error propagation
-------------------------------------------------------------------------------

We want to bound the error when `f(x+a)` is approximated by `f(x)`.
Specifically, the goal is to bound `f(x + a) - f(x)` in terms of `r`
for the set of values `a` with `|a| \le r`.
Most bounds will be monotone increasing with `|a|` (assuming that `x` is
fixed), so for brevity we simply express the bounds in terms of `|a|`.

**Theorem (generic first-order bound)**:

.. math::

    |f(x+a) - f(x)| \le \min(2 C_0, C_1 |a|)

where

.. math::

    C_0 = \sup_{|t| \le |a|} |f(x+t)|, \quad C_1 = \sup_{|t| \le |a|} |f'(x+t)|.

The statement is valid with either `a, t \in \mathbb{R}` or `a, t \in \mathbb{C}`.

**Theorem (product)**: For `x, y \in \mathbb{C}` and `a, b \in \mathbb{C}`,

.. math::

    \left| (x+a)(y+b) - x y \right| \le |xb| + |ya| + |ab|.

**Theorem (quotient)**: For `x, y \in \mathbb{C}` and `a, b \in \mathbb{C}`
with `|b| < |y|`,

.. math::

        \left| \frac{x}{y} - \frac{x+a}{y+b} \right|
        \le \frac{|xb|+|ya|}{|y| (|y|-|b|)}.

**Theorem (square root)**: For `x, a \in \mathbb{R}` with `0 \le |a| \le x`,

.. math::

    \left| \sqrt{x+a} - \sqrt{x} \right|
        \le \sqrt{x} \left(1 - \sqrt{1-\frac{|a|}{x}}\right)
        \le \frac{\sqrt{x}}{2} \left(\frac{|a|}{x} + \frac{|a|^2}{x^2}\right)

where the first inequality is an equality if `a \le 0`.
(When `x = a = 0`, the limiting value is 0.)

**Theorem (reciprocal square root)**: For `x, a \in \mathbb{R}` with `0 \le |a| < x`,

.. math::

    \left| \frac{1}{\sqrt{x+a}} - \frac{1}{\sqrt{x}} \right|
        \le \frac{|a|}{2 (x-|a|)^{3/2}}.

**Theorem (k-th root)**: For `k > 1` and `x, a \in \mathbb{R}` with `0 \le |a| \le x`,

.. math::

    \left| (x+a)^{1/k} - x^{1/k} \right|
        \le x^{1/k} \min\left(1, \frac{1}{k} \, \log\left(1+\frac{|a|}{x-|a|}\right)\right).

*Proof*: The error is largest when `a = -r` is negative, and

.. math::

    x^{1/k} - (x-r)^{1/k} &= x^{1/k} [1 - (1-r/x)^{1/k}]

    &= x^{1/k} [1 - \exp(\log(1-r/x)/k)] \le x^{1/k} \min(1, -\log(1-r/x)/k)

    &= x^{1/k} \min(1, \log(1+r/(x-r))/k).

**Theorem (sine, cosine)**: For `x, a \in \mathbb{R}`, `|\sin(x+a) - \sin(x)| \le \min(2, |a|)`.

**Theorem (logarithm)**: For `x, a \in \mathbb{R}` with `0 \le |a| < x`,

.. math::

    |\log(x+a) - \log(x)| \le \log\left(1 + \frac{|a|}{x-|a|}\right),

with equality if `a \le 0`.

**Theorem (exponential)**: For `x, a \in \mathbb{R}`, 
`|e^{x+a} - e^x| = e^x (e^a-1) \le e^x (e^{|a|}-1)`, with equality if `a \ge 0`.

**Theorem (inverse tangent)**: For `x, a \in \mathbb{R}`,

.. math::

    |\operatorname{atan}(x+a) - \operatorname{atan}(x)| \le \min(\pi, C_1 |a|).

where

.. math::

    C_1 = \sup_{|t| \le |a|} \frac{1}{1 + (x+t)^2}.

If `|a| < |x|`, then `C_1 = (1 + (|x| - |a|)^2)^{-1}` gives a monotone bound.

An exact bound: if `|a| < |x|` or `|x(x+a)| < 1`, then

.. math::

    |\operatorname{atan}(x+a) - \operatorname{atan}(x)| =
        \operatorname{atan}\left(\frac{|a|}{1 + x(x+a)}\right).

In the last formula, a case distinction has to be made depending on the
signs of *x* and *a*.

Sums and series
-------------------------------------------------------------------------------

**Theorem (geometric bound)**: If `|c_k| \le C` and `|z| \le D < 1`, then

.. math::

    \left| \sum_{k=N}^{\infty} c_k z^k \right| \le \frac{C D^N}{1 - D}.

**Theorem (integral bound)**: If `f(x)` is nonnegative and
monotone decreasing, then

.. math::

    \int_N^{\infty} f(x) \le \sum_{k=N}^{\infty} f(k) \le f(N) + \int_N^{\infty} f(x) dx.

Complex analytic functions
-------------------------------------------------------------------------------

**Theorem (Cauchy's integral formula)**:
If `f(z) = \sum_{k=0}^{\infty} c_k z^k` is analytic (on an open
subset of `\mathbb{C}` containing the disk `D = \{ z : |z| \le R \}`
in its interior, where `R > 0`), then

.. math::

    c_k = \frac{1}{2\pi i} \int_{|z|=R} \frac{f(z)}{z^{k+1}}\, dz.

**Corollary (derivative bound)**:

.. math::

    |c_k| \le \frac{C}{R^k}, \quad C = \max_{|z|=R} |f(z)|.

**Corollary (Taylor series tail)**:
If `0 \le r < R` and `|z| \le r`, then

.. math::

    \left|\sum_{k=N}^{\infty} c_k z^k\right| \le
        \frac{C D^N}{1-D}, \quad D = \left|\frac{r}{R}\right|.

Euler-Maclaurin formula
-------------------------------------------------------------------------------

**Theorem (Euler-Maclaurin)**:
If `f(t)` is `2M`-times differentiable, then

.. math::

    \sum_{k=L}^U f(k) = S + I + T + R

.. math::

    S = \sum_{k=L}^{N-1} f(k), \quad I = \int_N^U f(t) dt,

.. math::

    T = \frac{1}{2} \left( f(N) + f(U) \right) + 
        \sum_{k=1}^M \frac{B_{2k}}{(2k)!}
        \left(f^{(2k-1)}(U) - f^{(2k-1)}(N)\right),

.. math::

    R = -\int_N^U \frac{B_{2M}(t - \lfloor t \rfloor)}{(2M)!} f^{(2M)}(t) dt.

**Lemma (Bernoulli polynomials)**: `|B_n(t - \lfloor t \rfloor)| \le 4 n! / (2 \pi)^n`.

**Theorem (remainder bound)**:

.. math::

    |R| \le \frac{4}{(2\pi)^{2M}} \int_N^U \left| f^{(2M)}(t) \right| dt.

**Theorem (parameter derivatives)**:
If `f(t) = f(t,x) = \sum_{k=0}^{\infty} a_k(t) x^k` and
`R = R(x) = \sum_{k=0}^{\infty} c_k x^k`
are analytic functions of `x`, then

.. math::

    |c_k| \le \frac{4}{(2\pi)^{2M}} \int_N^U |a_k^{(2M)}(t)| dt.

