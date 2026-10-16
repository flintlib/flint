.. _partitions:

**partitions.h** -- computation of the partition function
===============================================================================

This module computes the integer partition function `p(n)`: by a lookup
table and the pentagonal recurrence for small `n`, and otherwise by the
asymptotically fast algorithm described in [Joh2012]_, which evaluates a
truncation of the Hardy-Ramanujan-Rademacher series using tight precision
estimates and symbolically factoring the occurring exponential sums.  The
evaluation is done by :func:`mp_real_partitions_hrr`, in ball arithmetic
throughout, so the result is proved correct; it uses hardware
double-based ball arithmetic (dfloat) for the low-precision terms and
several threads when they are available.  The previous implementation
based on arb balls (``partitions_hrr_sum_arb``) is kept as an example in
``examples/partitions.c``.

.. function:: void partitions_rademacher_bound(arf_t b, const fmpz_t n, ulong N)

    Sets `b` to an upper bound for

    .. math::

        M(n,N) = \frac{44 \pi^2}{225 \sqrt 3} N^{-1/2}
                  + \frac{\pi \sqrt{2}}{75} \left( \frac{N}{n-1} \right)^{1/2}
                \sinh\left(\frac{\pi}{N} \sqrt{\frac{2n}{3}}\right).

    This formula gives an upper bound for the truncation error in the
    Hardy-Ramanujan-Rademacher formula when the series is taken up
    to the term `t(n,N)` inclusive.

.. function:: void partitions_fmpz_fmpz(fmpz_t p, const fmpz_t n, int use_doubles)
              void partitions_fmpz_ui(fmpz_t p, ulong n)

    Computes the partition function `p(n)` (zero for negative `n`): from
    a table for `n < 128`, by the pentagonal recurrence while `p(n)` fits
    a word (`n < 417` on 64-bit machines), and otherwise by
    :func:`mp_real_partitions_hrr`, which computes a ball containing
    `p(n)` and verifies that it contains a unique integer.  The number of
    threads selected with :func:`flint_set_num_threads` is used for large
    `n` (from about `10^8`).  The *use_doubles* argument is ignored (the
    hardware doubles are always used, with proved error bounds).  The
    result must fit an :type:`fmpz_t` (`n` up to about `1.4 \cdot
    10^{21}`); use :func:`mp_real_partitions_hrr` beyond.

.. function:: void partitions_fmpz_ui_using_doubles(fmpz_t p, ulong n)

    Deprecated: the same as :func:`partitions_fmpz_ui`.

.. function:: void partitions_leading_fmpz(arb_t res, const fmpz_t n, slong prec)

    Sets *res* to the leading term in the Hardy-Ramanujan series
    for `p(n)` (without Rademacher's correction of this term, which is
    vanishingly small when `n` is large), that is,
    `\sqrt{12} (1-1/t) e^t / (24n-1)` where `t = \pi \sqrt{24n-1} / 6`.

