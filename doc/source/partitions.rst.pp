.. _partitions:

**partitions.h** -- computation of the partition function
===============================================================================

This module implements the asymptotically fast algorithm
for evaluating the integer partition function `p(n)`
described in [Joh2012]_.
The idea is to evaluate a truncation of the Hardy-Ramanujan-Rademacher series
using tight precision estimates, and symbolically factoring the
occurring exponential sums.

An implementation based on floating-point arithmetic can
also be found in FLINT. That version relies on some numerical subroutines
that have not been proved correct.

The implementation provided here uses ball arithmetic throughout to guarantee
a correct error bound for the numerical approximation of `p(n)`.
Optionally, hardware double arithmetic can be used for low-precision
terms. This gives a significant speedup for small (e.g. `n < 10^6`).

.. function:: void partitions_rademacher_bound(arf_t b, const fmpz_t n, ulong N)

    Sets `b` to an upper bound for

    .. math::

        M(n,N) = \frac{44 \pi^2}{225 \sqrt 3} N^{-1/2}
                  + \frac{\pi \sqrt{2}}{75} \left( \frac{N}{n-1} \right)^{1/2}
                \sinh\left(\frac{\pi}{N} \sqrt{\frac{2n}{3}}\right).

    This formula gives an upper bound for the truncation error in the
    Hardy-Ramanujan-Rademacher formula when the series is taken up
    to the term `t(n,N)` inclusive.

.. function:: void partitions_hrr_sum_arb(arb_t x, const fmpz_t n, slong N0, slong N, int use_doubles)

    Evaluates the partial sum `\sum_{k=N_0}^N t(n,k)` of the
    Hardy-Ramanujan-Rademacher series.

    If *use_doubles* is nonzero, doubles and the system's standard library math
    functions are used to evaluate the smallest terms. This significantly
    speeds up evaluation for small `n` (e.g. `n < 10^6`), and gives a small speed
    improvement for larger `n`, but the result is not guaranteed to be correct.
    In practice, the error is estimated very conservatively, and unless
    the system's standard library is broken, use of doubles can be considered
    safe. Setting *use_doubles* to zero gives a fully guaranteed
    bound.

.. function:: void partitions_fmpz_fmpz(fmpz_t p, const fmpz_t n, int use_doubles)

    Computes the partition function `p(n)` using the Hardy-Ramanujan-Rademacher
    formula. This function computes a numerical ball containing `p(n)`
    and verifies that the ball contains a unique integer.

    If *n* is sufficiently large and a number of threads greater than 1
    has been selected with :func:`flint_set_num_threads()`, the computation
    time will be reduced by using two threads.

    See :func:`partitions_hrr_sum_arb` for an explanation of the
    *use_doubles* option.

.. function:: void partitions_fmpz_ui(fmpz_t p, ulong n)

    Computes the partition function `p(n)` using the Hardy-Ramanujan-Rademacher
    formula. This function computes a numerical ball containing `p(n)`
    and verifies that the ball contains a unique integer.

.. function:: void partitions_fmpz_ui_using_doubles(fmpz_t p, ulong n)

    Computes the partition function `p(n)`, enabling the use of doubles
    internally. This significantly speeds up evaluation for small `n`
    (e.g. `n < 10^6`), but the error bounds are not certified
    (see remarks for :func:`partitions_hrr_sum_arb`).

.. function:: void partitions_leading_fmpz(arb_t res, const fmpz_t n, slong prec)

    Sets *res* to the leading term in the Hardy-Ramanujan series
    for `p(n)` (without Rademacher's correction of this term, which is
    vanishingly small when `n` is large), that is,
    `\sqrt{12} (1-1/t) e^t / (24n-1)` where `t = \pi \sqrt{24n-1} / 6`.

