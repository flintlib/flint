.. _hypgeom:

**hypgeom.h** -- support for hypergeometric series
===============================================================================

This module provides functions for high-precision evaluation of series
of the form

.. math::

    \sum_{k=0}^{n-1} \frac{A(k)}{B(k)} \prod_{j=1}^k \frac{P(k)}{Q(k)} z^k

where `A, B, P, Q` are polynomials. The present version only supports
`A, B, P, Q \in \mathbb{Z}[k]` (represented using the
FLINT *fmpz_poly_t* type). This module also provides functions
for high-precision evaluation of infinite series (`n \to \infty`),
with automatic, rigorous error bounding.

Note that we can standardize to `A = B = 1` by
setting `\tilde P(k) = P(k) A(k) B(k-1), \tilde Q(k) = Q(k) A(k-1) B(k)`.
However, separating out `A` and `B` is convenient and improves
efficiency during evaluation.


Strategy for error bounding
-------------------------------------------------------------------------------

We wish to evaluate `S(z) = \sum_{k=0}^{\infty} T(k) z^k` where `T(k)`
satisfies `T(0) = 1` and

.. math::

    T(k) = R(k) T(k-1) = \left( \frac{P(k)}{Q(k)} \right) T(k-1)

for given polynomials

.. math::

    P(k) = a_p k^p + a_{p-1} k^{p-1} + \ldots a_0

    Q(k) = b_q k^q + b_{q-1} k^{q-1} + \ldots b_0.

For convergence, we require `p < q`, or `p = q` with `|z| |a_p| < |b_q|`.
We also assume that `P(k)` and `Q(k)` have no roots among the positive
integers (if there are positive integer roots, the sum is either finite
or undefined). With these conditions satisfied, our goal is to find a
parameter `n \ge 0` such that

.. math::

    \left\lvert \sum_{k=n}^{\infty} T(k) z^k \right\rvert \le 2^{-d}.

We can rewrite the hypergeometric term ratio as

.. math::

    z R(k) = z \frac{P(k)}{Q(k)} =
        z \left( \frac{a_p}{b_q} \right) \frac{1}{k^{q-p}} F(k)

where

.. math::

    F(k) = \frac{
    1 + \tilde a_{1} / k + \tilde a_{2} / k^2 + \ldots + \tilde a_q / k^p
    }{
    1 + \tilde b_{1} / k + \tilde b_{2} / k^2 + \ldots + \tilde b_q / k^q
    } = 1 + O(1/k)

and where `\tilde a_i = a_{p-i} / a_p`, `\tilde b_i = b_{q-i} / b_q`.
Next, we define

.. math::

    C = \max_{1 \le i \le p} |\tilde a_i|^{(1/i)},
    \quad D = \max_{1 \le i \le q} |\tilde b_i|^{(1/i)}.

Now, if `k > C`, the magnitude of the numerator of `F(k)` is
bounded from above by

.. math::

    1 + \sum_{i=1}^p \left(\frac{C}{k}\right)^i \le 1 + \frac{C}{k-C}

and if `k > 2D`, the magnitude of the denominator of `F(k)` is bounded
from below by

.. math::

    1 - \sum_{i=1}^q \left(\frac{D}{k}\right)^i \ge 1 + \frac{D}{D-k}.

Putting the inequalities together gives the following bound,
valid for `k > K = \max(C, 2D)`:

.. math::

    |F(k)| \le \frac{k (k-D)}{(k-C)(k-2D)} = \left(1 + \frac{C}{k-C} \right)
    \left(1 + \frac{D}{k-2D} \right).

Let `r = q-p` and `\tilde z = |z a_p / b_q|`. Assuming
`k > \max(C, 2D, {\tilde z}^{1/r})`, we have

.. math::

    |z R(k)| \le G(k) = \frac{\tilde z F(k)}{k^r}

where `G(k)` is monotonically decreasing. Now we just need to find an
`n` such that `G(n) < 1` and for which `|T(n)| / (1 - G(n)) \le 2^{-d}`.
This can be done by computing a floating-point guess for `n` then
trying successively larger values.

This strategy leaves room for some improvement. For example, if
`\tilde b_1` is positive and large, the bound `B` becomes very pessimistic
(a larger positive `\tilde b_1` causes faster convergence,
not slower convergence).


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: hypgeom_struct

.. type:: hypgeom_t

    Stores polynomials *A*, *B*, *P*, *Q* and precomputed bounds,
    representing a fixed hypergeometric series.


Memory management
-------------------------------------------------------------------------------

.. function:: void hypgeom_init(hypgeom_t hyp)

.. function:: void hypgeom_clear(hypgeom_t hyp)


Error bounding
-------------------------------------------------------------------------------

.. function:: slong hypgeom_estimate_terms(const mag_t z, int r, slong d)

    Computes an approximation of the largest `n` such
    that `|z|^n/(n!)^r = 2^{-d}`, giving a first-order estimate of the
    number of terms needed to approximate the sum of a hypergeometric
    series of weight `r \ge 0` and argument `z` to an absolute
    precision of `d \ge 0` bits. If `r = 0`, the direct solution of the
    equation is given by `n = (\log(1-z) - d \log 2) / \log z`.
    If `r > 0`, using `\log n! \approx n \log n - n` gives an equation
    that can be solved in terms of the Lambert *W*-function as
    `n = (d \log 2) / (r\,W\!(t))` where
    `t = (d \log 2) / (e r z^{1/r})`.

    The evaluation is done using double precision arithmetic.
    The function aborts if the computed value of `n` is greater
    than or equal to LONG_MAX / 2.

.. function:: slong hypgeom_bound(mag_t error, int r, slong C, slong D, slong K, const mag_t TK, const mag_t z, slong prec)

    Computes a truncation parameter sufficient to achieve *prec* bits
    of absolute accuracy, according to the strategy described above.
    The input consists of `r`, `C`, `D`, `K`, precomputed bound for `T(K)`,
    and `\tilde z = z (a_p / b_q)`, such that for `k > K`, the hypergeometric
    term ratio is bounded by

    .. math::

        \frac{\tilde z}{k^r} \frac{k(k-D)}{(k-C)(k-2D)}.

    Given this information, we compute a `\varepsilon` and an
    integer `n` such that
    `\left| \sum_{k=n}^{\infty} T(k) \right| \le \varepsilon \le 2^{-\mathrm{prec}}`.
    The output variable *error* is set to the value of `\varepsilon`,
    and `n` is returned.

.. function:: void hypgeom_precompute(hypgeom_t hyp)

    Precomputes the bounds data `C`, `D`, `K` and an upper bound for `T(K)`.

Summation
-------------------------------------------------------------------------------

.. function:: void arb_hypgeom_sum(arb_t P, arb_t Q, const hypgeom_t hyp, slong n, slong prec)

    Computes `P, Q` such that `P / Q = \sum_{k=0}^{n-1} T(k)` where `T(k)`
    is defined by *hyp*,
    using binary splitting and a working precision of *prec* bits.

.. function:: void arb_hypgeom_infsum(arb_t P, arb_t Q, hypgeom_t hyp, slong tol, slong prec)

    Computes `P, Q` such that `P / Q = \sum_{k=0}^{\infty} T(k)` where `T(k)`
    is defined by *hyp*, using binary splitting and
    working precision of *prec* bits.
    The number of terms is chosen automatically to bound the
    truncation error by at most `2^{-\mathrm{tol}}`.
    The bound for the truncation error is included in the output
    as part of *P*.

