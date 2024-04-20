.. _mpn-extras:

**mpn_extras.h** -- support functions for limb arrays
===============================================================================

Macros
--------------------------------------------------------------------------------

.. macro:: MPN_NORM(a, an)

    Normalise ``(a, an)`` so that either ``an`` is zero or 
    ``a[an - 1]`` is nonzero.

.. macro:: MPN_SWAP(a, an, b, bn)

    Swap ``(a, an)`` and ``(b, bn)``, i.e. swap pointers and sizes.


Utility functions
--------------------------------------------------------------------------------


.. function:: void flint_mpn_debug(mp_srcptr x, mp_size_t xsize)

    Prints debug information about ``(x, xsize)`` to ``stdout``. 
    In particular, this will print binary representations of all the limbs.

.. function:: int flint_mpn_zero_p(mp_srcptr x, mp_size_t xsize)

    Returns `1` if all limbs of ``(x, xsize)`` are zero, otherwise `0`.

.. function:: mp_limb_t flint_mpn_sumdiff_n(mp_ptr s, mp_ptr d, mp_srcptr x, mp_srcptr y, mp_size_t n)

    Simultaneously computes the sum ``s`` and difference ``d`` of ``(x, n)`` and ``(y, n)``,
    returning carry multiplied by two plus borrow.

Multiplication
--------------------------------------------------------------------------------

.. function:: mp_limb_t flint_mpn_mul(mp_ptr z, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn)

    Sets ``(z, xn+yn)`` to the product of ``(x, xn)`` and ``(y, yn)``
    and returns the top limb of the result.
    We require `xn \ge yn \ge 1`
    and that ``z`` is not aliased with either input operand.
    This function is intended for all operand sizes. It will automatically
    select an appropriate algorithm out of the following:

    * A hardcoded multiplication function for small sizes.
    * Karatsuba or Toom-Cook multiplication for intermediate sizes.
    * FFT multiplication for huge sizes.
    * A GMP fallback for cases where we do currently not have optimized code.

.. function:: void flint_mpn_mul_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_size_t n)

    Sets ``z`` to the product of ``(x, n)`` and ``(y, n)``.
    We require `n \ge 1`
    and that ``z`` is not aliased with either input operand.
    The algorithm selection is similar to :func:`flint_mpn_mul`.

.. function:: void flint_mpn_sqr(mp_ptr z, mp_srcptr x, mp_size_t n)

    Sets ``z`` to the square of ``(x, n)``.
    We require `n \ge 1`
    and that ``z`` is not aliased with the input operand.
    The algorithm selection is similar to :func:`flint_mpn_sqr`.

.. function:: mp_size_t flint_mpn_fmms1(mp_ptr y, mp_limb_t a1, mp_srcptr x1, mp_limb_t a2, mp_srcptr x2, mp_size_t n)

    Given not-necessarily-normalized `x_1` and `x_2` of length `n > 0` and output `y` of length `n`, try to compute `y = a_1\cdot x_1 - a_2\cdot x_2`.
    Return the normalized length of `y` if `y \ge 0` and `y` fits into `n` limbs. Otherwise, return `-1`.
    `y` may alias `x_1` but is not allowed to alias `x_2`.

.. function:: void flint_mpn_mul_toom22(mp_ptr pp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn, mp_ptr scratch)

    Toom-22 (Karatsuba) multiplication. The *scratch* space must have room for
    `2 \text{an} + k` limbs where `k` is the number of limbs. If *NULL* is passed,
    space will be allocated internally.

Truncating multiplication
--------------------------------------------------------------------------------

Given two `n`-limb integers, a *high product* (or *mulhigh*) is an approximation
of the leading `n` limbs of the full `2n`-limb product.
In the basecase regime, a high product can be computed in roughly half the
time of the full product, and in some fraction `0.5 < c < 1` of the time
in the Toom-Cook regime. This speedup vanishes asymptotically in the FFT
regime. Contrary to polynomial high products or integer low products, integer
high products are not uniquely defined due to carry propagation.
We make the following definitions:

* *Rough mulhigh* accumulates at least `n + 1` limbs of partial products,
  outputting `n` limbs where the `n - 1` most significant limbs are essentially
  correct and the `n`-th most significant limb may have an error of `O(n)` ulp.
  This is the version of mulhigh used in [HZ2011]_.
* *Precise mulhigh* accumulates at least `n + 2` limbs of partial products,
  outputting `n + 1` limbs where the `n` most significant limbs are essentially
  correct and the `(n+1)`-th most significant limb may have an error of `O(n)` ulp.
* *Exact mulhigh* is the exact truncation of the full product. This cannot be
  computed faster than the full product in the worst case, but it can be
  computed faster on average by performing a precise mulhigh, inspecting
  the low output limb, and correcting with a low product when necessary.

In all cases, a high product is either equal to or smaller than the high part
of the full product.

More generally, we can define `n`-limb high products of `m`-limb and
`p`-limb integers where `m + p > n`, but this is not currently implemented.

.. function:: void _flint_mpn_mulhigh_n_mulders_recursive(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
              void _flint_mpn_sqrhigh_mulders_recursive(mp_ptr res, mp_srcptr u, mp_size_t n)

    Rough mulhigh implemented using Mulders' recursive algorithm as described in [HZ2011]_.
    Puts in *res[n], ..., res[2n-1]* an approximation of the `n` high limbs of *{u, n}* times *{v, n}*.
    The error is less than *n* ulps of *res[n]*. Assumes `2n` limbs are allocated at *res*;
    the low limbs will be used as scratch space.
    The *sqrhigh* version implements squaring.

.. function:: mp_limb_t _flint_mpn_mulhigh_basecase(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
              mp_limb_t _flint_mpn_mulhigh_n_mulders(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
              mp_limb_t _flint_mpn_mulhigh_n_mul(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
              mp_limb_t flint_mpn_mulhigh_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)

    Precise mulhigh. Puts in *res[0], ..., res[n-1]* an approximation of the `n` high limbs of
    *{u, n}* times *{v, n}*. and returns the `(n+1)`-th most significant limb.
    The error is at most *n + 2* ulp in the returned limb.

    * The *basecase* version implements the `O(n^2)` schoolbook algorithm.
      On x86-64 machines with ADX, the basecase version currently assumes
      that `n \ge 6`.
    * The *mulders* version computes a rough mulhigh with one extra limb of precision
      in temporary scratch space using :func:`_flint_mpn_mulhigh_n_mulders_recursive`
      and then copies the high limbs to the output.
    * The *mul* version computes a full product in temporary scratch space and
      copies the high limbs to the output. The output is actually the exact
      mulhigh.
    * The default version looks up a hardcoded basecase multiplication routine
      in a table for small *n*, and otherwise calls the *basecase*, *mulders*
      or *mul* implementations.

.. function:: mp_limb_t _flint_mpn_sqrhigh_basecase(mp_ptr res, mp_srcptr u, mp_size_t n)
              mp_limb_t _flint_mpn_sqrhigh_mulders(mp_ptr res, mp_srcptr u, mp_size_t n)
              mp_limb_t _flint_mpn_sqrhigh_sqr(mp_ptr res, mp_srcptr u, mp_size_t n)
              mp_limb_t flint_mpn_sqrhigh(mp_ptr res, mp_srcptr u, mp_size_t n)

    Squaring counterparts of :func:`flint_mpn_mulhigh_n`.

    On x86-64 machines with ADX, the basecase version currently assumes
    that `n \ge 8`.

Divisibility
--------------------------------------------------------------------------------


.. function:: int flint_mpn_divisible_1_odd(mp_srcptr x, mp_size_t xsize, mp_limb_t d)

    Expression determining whether ``(x, xsize)`` is divisible by the
    ``mp_limb_t d`` which is assumed to be odd-valued and at least `3`.

    This function is implemented as a macro.

.. function:: mp_size_t flint_mpn_remove_2exp(mp_ptr x, mp_size_t xsize, flint_bitcnt_t * bits)

    Divides ``(x, xsize)`` by `2^n` where `n` is the number of trailing 
    zero bits in `x`. The new size of `x` is returned, and `n` is stored in 
    the bits argument. `x` may not be zero.

.. function:: mp_size_t flint_mpn_remove_power_ascending(mp_ptr x, mp_size_t xsize, mp_ptr p, mp_size_t psize, ulong * exp)

    Divides ``(x, xsize)`` by the largest power `n` of ``(p, psize)`` 
    that is an exact divisor of `x`. The new size of `x` is returned, and 
    `n` is stored in the ``exp`` argument. `x` may not be zero, and `p` 
    must be greater than `2`.

    This function works by testing divisibility by ascending squares
    `p, p^2, p^4, p^8, \dotsc`, making it efficient for removing potentially
    large powers. Because of its high overhead, it should not be used as
    the first stage of trial division.

.. function:: int flint_mpn_factor_trial(mp_srcptr x, mp_size_t xsize, slong start, slong stop)

    Searches for a factor of ``(x, xsize)`` among the primes in positions 
    ``start, ..., stop-1`` of ``flint_primes``. Returns `i` if 
    ``flint_primes[i]`` is a factor, otherwise returns `0` if no factor 
    is found. It is assumed that ``start >= 1``.

.. function:: int flint_mpn_factor_trial_tree(slong * factors, mp_srcptr x, mp_size_t xsize, slong num_primes)

    Searches for a factor of ``(x, xsize)`` among the primes in positions
    approximately in the range ``0, ..., num_primes - 1`` of ``flint_primes``.
    
    Returns the number of prime factors found and fills ``factors`` with their
    indices in ``flint_primes``. It is assumed that ``num_primes`` is in the
    range ``0, ..., 3512``.

    If the input fits in a small ``fmpz`` the number is fully factored instead.

    The algorithm used is a tree based gcd with a product of primes, the tree
    for which is cached globally (it is threadsafe).

Division
--------------------------------------------------------------------------------


.. function:: int flint_mpn_divides(mp_ptr q, mp_srcptr array1, mp_size_t limbs1, mp_srcptr arrayg, mp_size_t limbsg, mp_ptr temp)

    If ``(arrayg, limbsg)`` divides ``(array1, limbs1)`` then
    ``(q, limbs1 - limbsg + 1)`` is set to the quotient and 1 is 
    returned, otherwise 0 is returned. The temporary space ``temp``
    must have space for ``limbsg`` limbs.

    Assumes ``limbs1 >= limbsg > 0``.

.. function:: mp_limb_t flint_mpn_preinv1(mp_limb_t d, mp_limb_t d2)

    Computes a precomputed inverse from the leading two limbs of the
    divisor ``b, n`` to be used with the ``preinv1`` functions.
    We require the most significant bit of ``b, n`` to be 1.

.. function:: mp_limb_t flint_mpn_divrem_preinv1(mp_ptr q, mp_ptr a, mp_size_t m, mp_srcptr b, mp_size_t n, mp_limb_t dinv)

    Divide ``a, m`` by ``b, n``, returning the high limb of the 
    quotient (which will either be 0 or 1), storing the remainder in-place 
    in ``a, n`` and the rest of the quotient in ``q, m - n``.
    We require the most significant bit of ``b, n`` to be 1.
    ``dinv`` must be computed from ``b[n - 1]``, ``b[n - 2]`` by 
    ``flint_mpn_preinv1``. We also require ``m >= n >= 2``.

.. function:: void flint_mpn_mulmod_preinv1(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_size_t n, mp_srcptr d, mp_limb_t dinv, ulong norm)

    Given a normalised integer `d` with precomputed inverse ``dinv`` 
    provided by ``flint_mpn_preinv1``, computes `ab \pmod{d}` and
    stores the result in `r`. Each of `a`, `b` and `r` is expected to 
    have `n` limbs of space, with zero padding if necessary. 

    The value ``norm`` is provided for convenience. If `a`, `b` and
    `d` have been shifted left by ``norm`` bits so that `d` is
    normalised, then `r` will be shifted right by ``norm`` bits
    so that it has the same shift as all the inputs.

    We require `a` and `b` to be reduced modulo `n` before calling the
    function. 

.. function:: void flint_mpn_preinvn(mp_ptr dinv, mp_srcptr d, mp_size_t n)

    Compute an `n` limb precomputed inverse ``dinv`` of the `n` limb
    integer `d`.

    We require that `d` is normalised, i.e. with the most significant
    bit of the most significant limb set.

.. function:: void flint_mpn_mod_preinvn(mp_ptr r, mp_srcptr a, mp_size_t m, mp_srcptr d, mp_size_t n, mp_srcptr dinv)

    Given a normalised integer `d` of `n` limbs, with precomputed inverse
    ``dinv`` provided by ``flint_mpn_preinvn`` and integer `a` of `m`
    limbs, computes `a \pmod{d}` and stores the result in-place in the lower
    `n` limbs of `a`. The remaining limbs of `a` are destroyed.

    We require `m \geq n`. No aliasing of `a` with any of the other operands
    is permitted.

    Note that this function is not always as fast as ordinary division.

.. function:: mp_limb_t flint_mpn_divrem_preinvn(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t m, mp_srcptr d, mp_size_t n, mp_srcptr dinv)

    Given a normalised integer `d` with precomputed inverse ``dinv`` 
    provided by ``flint_mpn_preinvn``, computes the quotient of `a` by `d` 
    and stores the result in `q` and the remainder in the lower `n` limbs of
    `a`. The remaining limbs of `a` are destroyed.

    The value `q` is expected to have space for `m - n` limbs and we require
    `m \ge n`. No aliasing is permitted between `q` and `a` or between these
    and any of the other operands. 

    Note that this function is not always as fast as ordinary division.

.. function:: void flint_mpn_mulmod_preinvn(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm)

    Given a normalised integer `d` with precomputed inverse ``dinv`` 
    provided by ``flint_mpn_preinvn``, computes `ab \pmod{d}` and
    stores the result in `r`. Each of `a`, `b` and `r` is expected to 
    have `n` limbs of space, with zero padding if necessary. 

    The value ``norm`` is provided for convenience. If `a`, `b` and
    `d` have been shifted left by ``norm`` bits so that `d` is
    normalised, then `r` will be shifted right by ``norm`` bits
    so that it has the same shift as all the inputs.

    We require `a` and `b` to be reduced modulo `n` before calling the
    function. 

    Note that this function is not always as fast as ordinary division.


GCD
--------------------------------------------------------------------------------


.. function:: mp_size_t flint_mpn_gcd_full2(mp_ptr arrayg, mp_srcptr array1, mp_size_t limbs1, mp_srcptr array2, mp_size_t limbs2, mp_ptr temp)

    Sets ``(arrayg, retvalue)`` to the gcd of ``(array1, limbs1)`` and
        ``(array2, limbs2)``.

    The only assumption is that neither ``limbs1`` nor ``limbs2`` is
    zero.

    The function must be supplied with ``limbs1 + limbs2`` limbs of temporary
    space, or ``NULL`` must be passed to ``temp`` if the function should
    allocate its own space.

.. function:: mp_size_t flint_mpn_gcd_full(mp_ptr arrayg, mp_srcptr array1, mp_size_t limbs1, mp_srcptr array2, mp_size_t limbs2)

    Sets ``(arrayg, retvalue)`` to the gcd of ``(array1, limbs1)`` and
    ``(array2, limbs2)``. 

    The only assumption is that neither ``limbs1`` nor ``limbs2`` is
    zero.


Random Number Generation
--------------------------------------------------------------------------------


.. function:: void flint_mpn_rrandom(mp_ptr rp, flint_rand_t state, mp_size_t n)

    Generates a random number with ``n`` limbs and stores 
    it on ``rp``. The number it generates will tend to have
    long strings of zeros and ones in the binary representation.

    Useful for testing functions and algorithms, since this kind of random
    numbers have proven to be more likely to trigger corner-case bugs.

.. function:: void flint_mpn_urandomb(mp_ptr rp, flint_rand_t state, flint_bitcnt_t n)

    Generates a uniform random number of ``n`` bits and stores 
    it on ``rp``.

