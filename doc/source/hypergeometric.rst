.. _algorithms_hypergeometric:

Algorithms for hypergeometric functions
===============================================================================

The algorithms used to compute hypergeometric functions are
described in [Joh2016]_. Here, we state the most important error bounds.

.. _algorithms_hypergeometric_convergent:

Convergent series
-------------------------------------------------------------------------------

Let

.. math::

    T(k) = \frac{\prod_{i=0}^{p-1} (a_i)_k}{\prod_{i=0}^{q-1} (b_i)_k} z^k.

We compute a factor *C* such that

.. math::

    \left|\sum_{k=n}^{\infty} T(k)\right| \le C |T(n)|.

We check that `\operatorname{Re}(b+n) > 0` for all lower
parameters *b*. If this does not hold, *C* is set to infinity.
Otherwise, we cancel out pairs of parameters
`a` and `b` against each other. We have

.. math::

    \left|\frac{a+k}{b+k}\right| = \left|1 + \frac{a-b}{b+k}\right| \le 1 + \frac{|a-b|}{|b+n|}

and

.. math::

    \left|\frac{1}{b+k}\right| \le \frac{1}{|b+n|}

for all `k \ge n`. This gives us a constant *D* such that
`T(k+1) \le D T(k)` for all `k \ge n`.
If `D \ge 1`, we set *C* to infinity. Otherwise, we take
`C = \sum_{k=0}^{\infty} D^k = (1-D)^{-1}`.

Convergent series of power series
-------------------------------------------------------------------------------

The same principle is used to get tail bounds for
with `a_i, b_i, z \in \mathbb{C}[[x]]`,
or more precisely, bounds for each coefficient in
`\sum_{k=N}^{\infty} T(k) \in \mathbb{C}[[x]] / \langle x^n \rangle`
given `a_i, b_i, z \in \mathbb{C}[[x]] / \langle x^n \rangle`.
First, we fix some notation, assuming that `A` and `B` are power series:

* `A_{[k]}` denotes the coefficient of `x^k` in `A`, and `A_{[m:n]}` denotes the power series `\sum_{k=m}^{n-1} A_{[k]} x^k`.
* `|A|` denotes `\sum_{k=0}^{\infty} |A_{[k]}| x^k` (this can be viewed as an element of `\mathbb{R}_{\ge 0}[[x]]`).
* `A \le B` signifies that `|A|_{[k]} \le |B|_{[k]}` holds for all `k`.
* We define `\mathcal{R}(B) = |B_{[0]}| - |B_{[1:\infty]}|`.

Using the formulas

.. math::

    (A B)_{[k]} = \sum_{j=0}^k A_{[j]} B_{[k-j]}, \quad (1 / B)_{[k]} = \frac{1}{B_{[0]}} \sum_{j=1}^k -B_{[j]} (1/B)_{[k-j]},

it is easy prove the following bounds for the coefficients
of sums, products and quotients of formal power series:

.. math::

    |A + B| \le |A| + |B|,
    \quad |A B|  \le |A| |B|,
    \quad |A / B| \le |A| / \mathcal{R}(B).

If `p \le q` and `\operatorname{Re}({b_i}_{[0]}+N) > 0` for all `b_i`, then we may take

.. math::

    D = |z| \, \prod_{i=1}^p \left(1 + \frac{|a_i-b_i|}{\mathcal{R}(b_i+N)}\right) \prod_{i=p+1}^{q} \frac{1}{\mathcal{R}(b_i + N)}.

If `D_{[0]} < 1`,then  `(1 - D)^{-1} |T(n)|` gives the error bound.

Note when adding and multiplying power series with (complex) interval coefficients,
we can use point-valued upper bounds for the absolute values instead
of performing interval arithmetic throughout.
For `\mathcal{R}(B)`, we must then pick a lower bound for `|B_{[0]}|` and upper bounds for
the coefficients of `|B_{[1:\infty]}|`.

.. _algorithms_hypergeometric_asymptotic_confluent:

Asymptotic series for the confluent hypergeometric function
-------------------------------------------------------------------------------

Let `U(a,b,z)` denote the confluent hypergeometric function of the second
kind with the principal branch cut, and
let `U^{*} = z^a U(a,b,z)`.
For all `z \ne 0` and `b \notin \mathbb{Z}` (but valid for all `b` as a limit),
we have (DLMF 13.2.42)

.. math::

    U(a,b,z)
        = \frac{\Gamma(1-b)}{\Gamma(a-b+1)} M(a,b,z)
        + \frac{\Gamma(b-1)}{\Gamma(a)} z^{1-b} M(a-b+1,2-b,z).

Moreover, for all `z \ne 0` we have

.. math::

    \frac{{}_1F_1(a,b,z)}{\Gamma(b)}
        = \frac{(-z)^{-a}}{\Gamma(b-a)} U^{*}(a,b,z)
        + \frac{z^{a-b} e^z}{\Gamma(a)} U^{*}(b-a,b,-z)

which is equivalent to DLMF 13.2.41 (but simpler in form).

We have the asymptotic expansion

.. math::

    U^{*}(a,b,z) \sim {}_2F_0(a, a-b+1, -1/z)

where `{}_2F_0(a,b,z)` denotes a formal hypergeometric series, i.e.

.. math::

    U^{*}(a,b,z) = \sum_{k=0}^{n-1} \frac{(a)_k (a-b+1)_k}{k! (-z)^k} + \varepsilon_n(z).

The error term `\varepsilon_n(z)` is bounded according to DLMF 13.7.
A case distinction is made depending on whether `z` lies in one
of three regions which we index by `R`.
Our formula for the error bound increases with the value of `R`, so we
can always choose the larger out of two indices if `z` lies in
the union of two regions.

Let `r = |b-2a|`.
If `\operatorname{Re}(z) \ge r`, set `R = 1`.
Otherwise, if `\operatorname{Im}(z) \ge r` or `\operatorname{Re}(z) \ge 0 \land |z| \ge r`, set `R = 2`.
Otherwise, if `|z| \ge 2r`, set `R = 3`.
Otherwise, the bound is infinite.
If the bound is finite, we have

.. math::

    |\varepsilon_n(z)| \le 2 \alpha C_n \left|\frac{(a)_n (a-b+1)_n}{n! z^n} \right| \exp(2 \alpha \rho C_1 / |z|)

in terms of the following auxiliary quantities

.. math::

    \sigma = |(b-2a)/z|

.. math::

    C_n = \begin{cases}
    1                              & \text{if } R = 1 \\
    \chi(n)                        & \text{if } R = 2 \\
    (\chi(n) + \sigma \nu^2 n) \nu^n & \text{if } R = 3
    \end{cases}

.. math::

    \nu = \left(\tfrac{1}{2} + \tfrac{1}{2}\sqrt{1-4\sigma^2}\right)^{-1/2} \le 1 + 2 \sigma^2

.. math::

    \chi(n) = \sqrt{\pi} \Gamma(\tfrac{1}{2}n+1) / \Gamma(\tfrac{1}{2} n + \tfrac{1}{2})

.. math::

    \sigma' = \begin{cases}
    \sigma & \text{if } R \ne 3 \\
    \nu \sigma & \text{if } R = 3
    \end{cases}

.. math::

    \alpha = (1 - \sigma')^{-1}

.. math::

    \rho = \tfrac{1}{2} |2a^2-2ab+b| + \sigma' (1+ \tfrac{1}{4} \sigma') (1-\sigma')^{-2}

.. _algorithms_hypergeometric_asymptotic_airy:

Asymptotic series for Airy functions
-------------------------------------------------------------------------------

Error bounds are based on Olver (DLMF section 9.7).
For `\arg(z) < \pi` and `\zeta = (2/3) z^{3/2}`, we have

.. math::

    \operatorname{Ai}(z) = \frac{e^{-\zeta}}{2 \sqrt{\pi} z^{1/4}} \left[S_n(\zeta) + R_n(z)\right], \quad
    \operatorname{Ai}'(z) = -\frac{z^{1/4} e^{-\zeta}}{2 \sqrt{\pi}} \left[(S'_n(\zeta) + R'_n(z)\right]

.. math::

    S_n(\zeta) = \sum_{k=0}^{n-1} (-1)^k \frac{u(k)}{\zeta^k}, \quad
    S'_n(\zeta) = \sum_{k=0}^{n-1} (-1)^k \frac{v(k)}{\zeta^k}

.. math::

    u(k) = \frac{(1/6)_k (5/6)_k}{2^k k!}, \quad
    v(k) = \frac{6k+1}{1-6k} u(k).

Assuming that *n* is positive, the error terms are bounded by

.. math::

    |R_n(z)|  \le C |u(n)| |\zeta|^{-n}, \quad |R'_n(z)| \le C |v(n)| |\zeta|^{-n}

where

.. math::

    C = \begin{cases}
        2 \exp(7 / (36 |\zeta|)) & |\arg(z)| \le \pi/3 \\
        2 \chi(n) \exp(7 \pi / (72 |\zeta|)) & \pi/3 \le |\arg(z)| \le 2\pi/3 \\
        4 \chi(n) \exp(7 \pi / (36 |\operatorname{re}(\zeta)|)) |\cos(\arg(\zeta))|^{-n} & 2\pi/3 \le |\arg(z)| < \pi.
        \end{cases}

For computing Bi when *z* is roughly in the positive half-plane, we use the
connection formulas

.. math::

    \operatorname{Bi}(z) = -i (2 w^{+1} \operatorname{Ai}(z w^{-2}) - \operatorname{Ai}(z))

    \operatorname{Bi}(z) = +i (2 w^{-1} \operatorname{Ai}(z w^{+2}) - \operatorname{Ai}(z))

where `w = \exp(\pi i/3)`. Combining roots of unity gives

.. math::

    \operatorname{Bi}(z) = \frac{1}{2 \sqrt{\pi} z^{1/4}} [2X + iY]

.. math::

    \operatorname{Bi}(z) = \frac{1}{2 \sqrt{\pi} z^{1/4}} [2X - iY]

.. math::

    X = \exp(+\zeta) [S_n(-\zeta) + R_n(z w^{\mp 2})], \quad Y = \exp(-\zeta) [S_n(\zeta) + R_n(z)]

where the upper formula is valid
for `-\pi/3 < \arg(z) < \pi` and the lower formula is valid for `-\pi < \arg(z) < \pi/3`.
We proceed analogously for the derivative of Bi.

In the negative half-plane, we use the connection formulas

.. math::

    \operatorname{Ai}(z) = e^{+\pi i/3} \operatorname{Ai}(z_1)  +  e^{-\pi i/3} \operatorname{Ai}(z_2)

.. math::

    \operatorname{Bi}(z) = e^{-\pi i/6} \operatorname{Ai}(z_1)  +  e^{+\pi i/6} \operatorname{Ai}(z_2)

where `z_1 = -z e^{+\pi i/3}`, `z_2 = -z e^{-\pi i/3}`.
Provided that `|\arg(-z)| < 2 \pi / 3`, we have
`|\arg(z_1)|, |\arg(z_2)| < \pi`, and thus the asymptotic expansion
for Ai can be used. As before, we collect roots of unity to obtain

.. math::

    \operatorname{Ai}(z) = A_1 [S_n(i \zeta)  + R_n(z_1)]
                         + A_2 [S_n(-i \zeta) + R_n(z_2)]

.. math::

    \operatorname{Bi}(z) = A_3 [S_n(i \zeta)  + R_n(z_1)]
                         + A_4 [S_n(-i \zeta) + R_n(z_2)]

where `\zeta = (2/3) (-z)^{3/2}` and

.. math::

    A_1 = \frac{\exp(-i (\zeta - \pi/4))}{2 \sqrt{\pi} (-z)^{1/4}}, \quad
    A_2 = \frac{\exp(+i (\zeta - \pi/4))}{2 \sqrt{\pi} (-z)^{1/4}}, \quad
    A_3 = -i A_1, \quad
    A_4 = +i A_2.

The differentiated formulas are analogous.

Corner case of the Gauss hypergeometric function
-------------------------------------------------------------------------------

In the corner case where `z` is near `\exp(\pm \pi i / 3)`, none of the
linear fractional transformations is effective.
In this case, we use Taylor series to analytically continue the solution
of the hypergeometric differential equation from the origin.
The function `f(z) = {}_2F_1(a,b,c,z_0+z)` satisfies

.. math::

    f''(z) = -\frac{((z_0+z)(a+b+1)-c)}{(z_0+z)(z_0-1+z)} f'(z) - \frac{a b}{(z_0+z)(z_0-1+z)} f(z).

Knowing `f(0), f'(0)`, we can compute the consecutive derivatives
recursively, and evaluating the truncated Taylor series allows us to
compute `f(z), f'(z)` to high accuracy
for sufficiently small `z`.
Some experimentation showed that two continuation steps

.. math::

    0 \quad \to \quad 0.375 \pm 0.625i \quad \to \quad 0.5 \pm 0.8125i \quad \to \quad z

gives good performance.
Error bounds for the truncated Taylor series are obtained
using the Cauchy-Kovalevskaya majorant method,
following the outline in [Hoe2001]_.
The differential equation is majorized by

.. math::

    g''(z) = \frac{N+1}{2} \left( \frac{\nu}{1-\nu z} \right) g'(z)
    + \frac{(N+1)N}{2} \left( \frac{\nu}{1-\nu z} \right)^2 g(z)

provided that `N` and `\nu \ge \max(1/|z_0|, 1/|z_0-1|)`
are chosen sufficiently large. It follows that we can compute explicit
numbers `A, N, \nu` such that the simple solution `g(z) = A (1-\nu z)^{-N}`
of the differential equation provides the bound

.. math::

    |f_{[k]}| \le g_{[k]} = A {{N+k} \choose k} \nu^k.

