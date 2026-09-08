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

.. function:: char * flint_mpn_get_str(char * res, int base, mp_srcptr x, mp_size_t xn, int negative)

    Returns the string representation of ``(x, xn)`` (or its negation if
    ``negative`` is set to 1) in base *base* which must be a base
    supported by GMP. If ``res`` is ``NULL``, a new string will be allocated;
    otherwise, the given pointer ``res`` will be used and is assumed
    to have sufficient space to represent the full output, one extra
    digit, minus sign (if negative), and null terminator.

.. function:: int flint_mpn_zero_p(mp_srcptr x, mp_size_t xsize)

    Returns `1` if all limbs of ``(x, xsize)`` are zero, otherwise `0`.

.. function:: int flint_mpn_equal_p(mp_srcptr x, mp_srcptr y, mp_size_t xsize)

    Returns `1` if all limbs of ``(x, xsize)`` and ``(y, xsize)`` are equal, otherwise `0`.

Addition and subtraction
--------------------------------------------------------------------------------

.. function:: mp_limb_t flint_mpn_sumdiff_n(mp_ptr s, mp_ptr d, mp_srcptr x, mp_srcptr y, mp_size_t n)

    Simultaneously computes the sum ``s`` and difference ``d`` of ``(x, n)`` and ``(y, n)``,
    returning carry multiplied by two plus borrow.

.. function:: void flint_mpn_negmod_n(mp_ptr res, mp_srcptr x, mp_srcptr m, mp_size_t n)
              void flint_mpn_addmod_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m, mp_size_t n)
              void flint_mpn_submod_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m, mp_size_t n)
              void flint_mpn_addmod_n_m(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t yn, mp_srcptr m, mp_size_t n)
              void flint_mpn_submod_n_m(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t yn, mp_srcptr m, mp_size_t n)

    Arithmetic modulo ``(m, n)``. These functions assume that
    ``(x, n)`` and ``(y, n)`` are already reduced modulo ``(m, n)``.
    The ``n_m`` variants accept ``(y, yn)`` with ``yn <= n``,
    where ``(y, yn)`` is already reduced modulo ``(m, n)``.

.. function:: void flint_mpn_negmod_2(mp_ptr res, mp_srcptr x, mp_srcptr m)
              void flint_mpn_addmod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
              void _flint_mpn_addmod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
              void flint_mpn_submod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)

    Modular arithmetic specialized for two limbs.
    The ``_flint_mpn_addmod_2`` version assumes that the most significant
    bit of ``m[1]`` is not set.

.. function:: int flint_mpn_signed_sub_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t n)

    Sets ``res`` to `|x - y|`, returning 0 if the result equals `x - y`
    and returning 1 if the result equals `y - x`.


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

.. function:: void _flint_mpn_mullow_n_mulders_recursive(mp_ptr rp, mp_srcptr u, mp_srcptr v, mp_size_t n)
              mp_limb_t flint_mpn_mullow_basecase(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
              mp_limb_t _flint_mpn_mullow_n_mulders(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
              mp_limb_t _flint_mpn_mullow_n_mul(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
              mp_limb_t _flint_mpn_mullow_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
              mp_limb_t flint_mpn_mullow_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)

    Compute the low `n` limbs of the product.

    The `(n + 1)`-th limb is also computed and returned.
    Warning: this extra limb of output may be removed in the future.

.. function:: void flint_mpn_mul_or_mullow_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)

    Write the low `n + 1` limbs of the product `uv` to ``res``.
    The output is assumed to have space for `2n` limbs so that the high
    limbs can be used as scratch space or to write the whole product
    when this is the fastest method.

    Warning: the one extra limb of output may be removed in the future.

.. function:: void flint_mpn_mul_or_mulhigh_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)

    Write the high `n + 1` limbs of the product `uv` to ``res + (n - 1)``
    (with possible error of a few ulps as for :func:`flint_mpn_mulhigh_n`).
    The low `n - 1` limbs of the output may be used as scratch space or
    to write the whole product when this is the fastest method.

Middle product
--------------------------------------------------------------------------------

The *windowed middle product* extracts a chosen limb window of a product.  For
`\mathrm{an} \ge 1`, `\mathrm{bn} \ge 1` and `0 \le \mathrm{zlo} < \mathrm{zhi}
\le \mathrm{an} + \mathrm{bn}`, it writes `\mathrm{zhi} - \mathrm{zlo}` limbs to
``z`` approximating limbs `[\mathrm{zlo}, \mathrm{zhi})` of `a b`.  It is a
*lower approximation*: partial products `a[p] b[q]` with `p + q < \mathrm{zlo}`
are dropped, so the computed value never exceeds the exact window, and the
deficit (a single carry from below `\mathrm{zlo}`) is bounded by
`\min(\mathrm{an}, \mathrm{bn}, \mathrm{zlo}) \cdot 2^{64}`.  With
`\mathrm{zlo} = 0` the window is exact.

.. function:: void flint_mpn_mulmid(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)

    Compute the window `[\mathrm{zlo}, \mathrm{zhi})` of `a b`, dispatching to
    whichever of the routines below is expected to be fastest for the given
    shape.  Individual backends may return the exact window or a tighter
    approximation than the classical drop; all satisfy the contract above.

.. function:: void flint_mpn_mulmid_classical(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)

    Row-based schoolbook implementation, `O((\mathrm{zhi} - \mathrm{zlo}) \cdot
    \min(\mathrm{an}, \mathrm{bn}))`.

.. function:: void flint_mpn_mulmid_via_mul(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)
              void flint_mpn_mulmid_via_mullow_n(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)
              void flint_mpn_mulmid_via_mulhigh_n(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)
              void flint_mpn_mulmid_via_n_padded(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)
              void flint_mpn_mulmid_fft_small(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)

    Reductions of a general window to, respectively, a full product
    (:func:`flint_mpn_mul`), a balanced low product (:func:`flint_mpn_mullow_n`),
    a balanced high product (:func:`flint_mpn_mulhigh_n`), a balanced middle
    product (:func:`flint_mpn_mulmid_n`) and the small-prime FFT.  Each is valid
    for arbitrary input by padding internally, but is only economical in its own
    regime; :func:`flint_mpn_mulmid` chooses between them.

.. function:: void flint_mpn_mulmid_n(mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)

    Exact balanced middle product of `\{ap, 2n-1\}` and `\{bp, n\}`, writing
    `n + 2` limbs: the high `n` limbs are exact and the low two are guard limbs
    (they lack the carry into the band from below).  This is a wrapper around
    GMP's ``mpn_mulmid_n`` and is only defined when
    ``FLINT_HAVE_NATIVE_mpn_mulmid_n`` is set (that is, when the build may call
    GMP internals and GMP exports the symbol; see ``configure``).  When it is
    unavailable, :func:`flint_mpn_mulmid_via_n_padded` is likewise unavailable
    and :func:`flint_mpn_mulmid` uses its other methods.

Dot products and polynomial multiplication with fixed limb sizes
--------------------------------------------------------------------------------

The following routines operate on arrays of contiguous, homogeneous
integers of fixed limb size, e.g. representing vectors or polynomials
whose coefficients have been packed into a uniform number of limbs.
Results are computed modulo `2^{64 s}` for a fixed number of output
limbs `s`, which makes them exact when the true results (as unsigned
integers, or as two's complement integers in the signed versions) are
known to fit in `s` limbs. The dot product kernels are generated by
``dev/gen_mpn_dot_rev.py``, using a separate short carry chain for each
limb offset to maximize instruction-level parallelism while keeping the
number of live registers small.

.. type:: flint_mpn_dot_func_t
          flint_mpn_dot_strided_func_t

    Function pointer types ``void (*)(nn_ptr, nn_srcptr, nn_srcptr, slong)``
    and ``void (*)(nn_ptr, nn_srcptr, slong, nn_srcptr, slong, slong)``
    for the dot product kernels below.

.. function:: void _mpn_dot_n1xn2_s(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len)
              void _mpn_dot_rev_n1xn2_s(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len)
              void _mpn_dot_strided_n1xn2_s(nn_ptr res, nn_srcptr a, slong astride, nn_srcptr b, slong bstride, slong len)
              void _mpn_dot_n1xn2_s_signed(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len)
              void _mpn_dot_rev_n1xn2_s_signed(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len)
              void _mpn_dot_strided_n1xn2_s_signed(nn_ptr res, nn_srcptr a, slong astride, nn_srcptr b, slong bstride, slong len)

    Sets ``{res, s}`` to `\sum_{k=0}^{len-1} a_k b_k \bmod 2^{64 s}`
    where the `a_k` are ``n1``-limb and the `b_k` are ``n2``-limb unsigned
    (or in the signed versions, two's complement) integers, for all
    `1 \le n_2 \le n_1 \le 4` and
    `s \in \{n_1 + n_2 - 1, n_1 + n_2, n_1 + n_2 + 1\}` with `s \le 8`.
    In the forward versions the entries are contiguous, `a_k = a[k]` and
    `b_k = b[k]`; in the reversed versions `b_k = b[len - 1 - k]`; in the
    strided versions the entries are ``astride`` respectively ``bstride``
    limbs apart (the strides may be negative). For example,
    :func:`_mpn_dot_2x2_5`, :func:`_mpn_dot_rev_2x1_3_signed` and
    :func:`_mpn_dot_strided_4x4_8` are kernels of this kind.

    For `n_1 \le 2`, dedicated forward and reversed kernels are generated;
    for `n_1 = 3, 4` these are inline wrappers around the strided kernels,
    where the stride overhead is negligible. For `n_1 \le 3` the products
    are evaluated with independent carry chains per limb offset; for
    `n_1 = 4` the assembly multiplication routines are used with the
    products accumulated in registers. In the signed versions, the top
    limb corrections are done without branches.

.. var:: const flint_mpn_dot_strided_func_t flint_mpn_dot_strided_tab[2][FLINT_MPN_DOT_TAB_N + 1][FLINT_MPN_DOT_TAB_N + 1][3]
         const flint_mpn_dot_func_t flint_mpn_dot_tab[2][FLINT_MPN_DOT_DEDICATED_TAB_N + 1][FLINT_MPN_DOT_DEDICATED_TAB_N + 1][3]
         const flint_mpn_dot_func_t flint_mpn_dot_rev_tab[2][FLINT_MPN_DOT_DEDICATED_TAB_N + 1][FLINT_MPN_DOT_DEDICATED_TAB_N + 1][3]

    Tables of the above kernels, indexed by ``[sgn][n1][n2][s - (n1 + n2 - 1)]``
    with ``NULL`` entries where no kernel exists (in particular when
    ``n1 < n2``). ``FLINT_MPN_DOT_TAB_N`` is 4 and
    ``FLINT_MPN_DOT_DEDICATED_TAB_N`` is 2.

.. function:: void _flint_mpn_dot_rev_generic(nn_ptr res, nn_srcptr a, slong n1, nn_srcptr b, slong n2, slong len, slong s, nn_ptr scratch)
              void _flint_mpn_dot_rev_generic_signed(nn_ptr res, nn_srcptr a, slong n1, nn_srcptr b, slong n2, slong len, slong s, nn_ptr scratch)

    Generic versions using :func:`flint_mpn_mul` and :func:`mpn_add_n`,
    valid for any `n_1 \ge n_2 \ge 1` and
    `s \in \{n_1 + n_2 - 1, n_1 + n_2, n_1 + n_2 + 1\}`. The scratch space
    must have at least `n_1 + n_2` limbs. The signed version handles
    the sign corrections with branches and is only competitive when
    the signs are predictable.

.. function:: void _flint_mpn_poly_mulmid_classical(nn_ptr res, nn_srcptr f, slong flen, slong n1, nn_srcptr g, slong glen, slong n2, slong nlo, slong nhi, slong s, int sgn)

    Sets ``{res, (nhi - nlo) s}`` to the coefficients `[nlo, nhi)` of the
    product of ``{f, flen n1}`` and ``{g, glen n2}``, viewed as polynomials
    with ``n1``-limb respectively ``n2``-limb coefficients (unsigned if
    ``sgn`` is 0, two's complement if ``sgn`` is 1), computed modulo
    `2^{64 s}` by classical multiplication, where `n_1 \ge n_2` and
    `s \in \{n_1 + n_2 - 1, n_1 + n_2, n_1 + n_2 + 1\}`. The results are
    exact provided that the requested coefficients fit in `s` limbs.
    Squaring is detected when ``f == g`` and ``flen == glen``.

    Each output coefficient is a dot product. For `n_1 \le 3` the routine
    dispatches to loop instances with the dot product code inlined
    (see below); for `n_1 = 3, 4` and short lengths, to instances using
    the assembly multiplication routines with register accumulation;
    and otherwise to loops using the dot product kernels or
    :func:`flint_mpn_mul` with :func:`mpn_add_n`.

.. function:: void _flint_mpn_poly_mulmid_n1xn2_s(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlo, slong nhi)
              void _flint_mpn_poly_mulmid_n1xn2_s_signed(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlo, slong nhi)
              void _flint_mpn_poly_sqrmid_nxn_s(nn_ptr res, nn_srcptr f, slong flen, slong nlo, slong nhi)
              void _flint_mpn_poly_sqrmid_nxn_s_signed(nn_ptr res, nn_srcptr f, slong flen, slong nlo, slong nhi)
              void _flint_mpn_poly_mulmid_n1xn2_s_short(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlo, slong nhi)
              void _flint_mpn_poly_sqrmid_nxn_s_short(nn_ptr res, nn_srcptr f, slong flen, slong nlo, slong nhi)

    Instances of :func:`_flint_mpn_poly_mulmid_classical` (and of squaring)
    for fixed limb sizes, with the coefficient loop and the dot product
    code fused so that no function calls are made per term or per
    coefficient, and for squaring, so that the doubling of the cross
    terms and the addition of the square term are done in registers.
    They exist for `n_1 \le 3` (for example :func:`_flint_mpn_poly_mulmid_2x2_5`
    and :func:`_flint_mpn_poly_sqrmid_3x3_7_signed`). For `n_1 \le 2` each term's
    full product is computed with a minimal number of additions and added
    to a single accumulator; for `n_1 = 3` one carry chain per limb offset
    is kept across the dot product as in the kernels, with the operands
    read from memory by the multiplication instructions to relieve
    register pressure. The ``_short`` versions, for `n_1 = 4` with
    `s \le 8`, use the assembly multiplication routines and accumulate
    in registers. Generated by ``dev/gen_mpn_dot_rev.py``.

.. function:: void _flint_mpn_poly_mul_classical(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlimbs, slong slimbs, int sgn)
              void _flint_mpn_poly_mul_karatsuba(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlimbs, slong slimbs, slong cutoff, int norm, int sgn)

    Sets ``{res, (flen + glen - 1) slimbs}`` to the full product of
    ``{f, flen nlimbs}`` and ``{g, glen nlimbs}``, viewed as polynomials
    with ``nlimbs``-limb coefficients (unsigned if ``sgn`` is 0,
    two's complement if ``sgn`` is 1), using classical respectively
    Karatsuba multiplication. The output coefficients are computed
    modulo `2^{64 \cdot slimbs}` where ``slimbs`` is one of
    `2 \cdot nlimbs - 1`, `2 \cdot nlimbs`, `2 \cdot nlimbs + 1`.
    Squaring is detected when ``f == g`` and ``flen == glen``.

    The results are exact provided that all coefficients of the product
    fit in ``slimbs`` limbs. For Karatsuba, this must also hold for all
    intermediate coefficients (the products of sums of the parts,
    which grow by up to two bits per level of recursion), and
    ``norm`` must bound the number of leading unused bits (excluding the
    sign bit when ``sgn`` is 1) in the top limb of the input coefficients.
    When ``norm`` reaches zero in the recursion, the coefficients of
    the sums are extended to ``nlimbs + 1`` limbs, so ``slimbs`` must
    also be at least `2 \cdot nlimbs + 1` in that case. Karatsuba recursion
    switches to the classical algorithm when the shorter length is
    below ``cutoff``; the classical instances are resolved once per
    level of recursion, so that the base cases are called directly.
    Unbalanced lengths are handled by splitting the longer operand into
    pieces of the length of the shorter one at every level, so that
    the cost is proportional to the length ratio.

.. function:: void _flint_mpn_poly_mulmid_karatsuba(nn_ptr res, nn_srcptr f, slong flen, slong nlimbs1, slong norm1, nn_srcptr g, slong glen, slong nlimbs2, slong norm2, slong nlo, slong nhi, slong slimbs, slong cutoff, int sgn)

    Sets ``{res, (nhi - nlo) slimbs}`` to the coefficients `[nlo, nhi)`
    of the product of ``{f, flen nlimbs1}`` and ``{g, glen nlimbs2}``,
    with the same conventions as :func:`_flint_mpn_poly_mul_karatsuba` except
    that the operands may have different numbers of limbs, with
    ``norm1`` and ``norm2`` their respective headroom; ``slimbs`` must be
    at least ``nlimbs1 + nlimbs2 - 1``, and at least
    ``nlimbs1 + nlimbs2 + 1`` if either operand may get extended in the
    recursion (when its headroom is smaller than the recursion depth
    plus one). Squaring is detected when ``f == g``, ``flen == glen``
    and ``nlimbs1 == nlimbs2``, also for the short products.
    The coefficients of the inputs that cannot contribute to the window
    are discarded first, and the algorithm is chosen according to the
    shape:

    - a full product (as above) if the window is the whole product;
    - Mulders' short product for a low product (`nlo = 0`): a full
      product of the low parts, split at about `0.7 n`, plus two recursive
      short products (one, when squaring); a high product (`nhi` equal to
      the product length) is computed as the reversal of a short product
      of the reversed inputs;
    - the Karatsuba middle product of Hanrot, Quercia and Zimmermann for
      windows in the "middle" region where all coefficients of the shorter
      operand contribute (for example the `n` coefficients `[n-1, 2n-1)` of a
      `(2n-1) \times n` product), computed in blocks of the length of
      the shorter operand, with the edges computed classically;
    - otherwise a full product of the trimmed inputs, clipped to the
      window, or classical multiplication for small windows.

    The short products are used from about 4 times ``cutoff``, below
    which the classical windowed algorithm is faster. The Karatsuba
    middle product uses two's complement arithmetic internally (it
    involves differences of parts of the longer operand), so ``slimbs``
    must include a sign bit even for unsigned inputs; moreover, it
    consumes one bit of headroom on every branch of the recursion, so
    when ``norm`` is smaller than the recursion depth the coefficients
    get extended early (about doubling the cost) and a larger cutoff is
    used.

.. macro:: FLINT_MPN_POLY_MUL_UNSIGNED
           FLINT_MPN_POLY_MUL_SIGNED
           FLINT_MPN_POLY_MUL_SIGNMAG
           FLINT_MPN_POLY_MUL_BIAS

    Representations of the packed coefficients accepted by
    :func:`_flint_mpn_poly_mulmid`: nonnegative integers (zero-extended);
    two's complement integers; sign-magnitude, i.e. a sign limb (0 or 1)
    followed by the magnitude, so that a coefficient of `n` magnitude
    limbs occupies `n + 1` limbs (only available with the classical
    algorithm); or biased, i.e. `x + 2^b` where `b` is the bit bound of
    the operand, which the multiplication treats as unsigned, leaving
    it to the caller to correct the outputs using sliding window sums
    of the inputs (this costs `O(\mathtt{slimbs})` per output coefficient
    but nothing per term, so it wins for long dot products and clearly
    with Karatsuba).

.. type:: mpn_poly_mul_params_struct
          mpn_poly_mul_params_t

    Packing parameters for :func:`_flint_mpn_poly_mulmid`: the representation
    ``method`` (one of the above), the numbers of limbs ``nlimbs1`` and
    ``nlimbs2`` per coefficient of the two operands (the magnitude limbs
    in the sign-magnitude case), their headroom ``norm1`` and ``norm2``
    (leading unused bits, not counting the sign bit), and the number of
    output limbs ``slimbs``.

.. function:: void _flint_mpn_poly_mulmid_params(mpn_poly_mul_params_t P, slong len1, slong bits1, slong len2, slong bits2, slong nlo, slong nhi, int squaring)

    Chooses packing parameters that are correct and (as far as the tuning
    goes) optimal for computing the coefficients `[nlo, nhi)` of the
    product of two polynomials of lengths ``len1`` and ``len2`` with
    :func:`_flint_mpn_poly_mulmid`, given bounds for the coefficients in the
    format of :func:`_fmpz_vec_max_bits`: the coefficients of the
    operands have at most ``|bits1|`` and ``|bits2|`` bits, and a
    negative bound indicates that negative coefficients may be present.
    ``squaring`` should be set if the operands are the same.

    For nonnegative coefficients the representation is always unsigned,
    with the minimal number of limbs. Otherwise the choice between
    two's complement, biased and sign-magnitude coefficients depends on
    the limb counts, the lengths and on whether Karatsuba will be used:
    the two's complement kernels are somewhat slower than the unsigned
    ones and are only available in the classical algorithm for small
    limb counts, the bias method has an overhead per output coefficient,
    and when the sign bit would cost a limb (coefficients of exactly
    `64 k` bits), the classical sign-magnitude kernels remain the
    fastest option up to many times the Karatsuba cutoff, since the
    growth in the Karatsuba recursion forces the extra limb in any
    representation.

.. function:: slong _flint_mpn_poly_mulmid_slimbs(slong minlen, slong nlimbs1, slong norm1, slong nlimbs2, slong norm2, int method, slong cutoff)

    Returns the number of output limbs required by :func:`_flint_mpn_poly_mulmid`
    for the given packed operands (see :type:`mpn_poly_mul_params_t`)
    when the dot products have at most ``minlen`` terms and the Karatsuba
    cutoff is ``cutoff``: the bound for the classical algorithm, or, if
    Karatsuba is used (``minlen`` at least ``cutoff``, except for
    sign-magnitude coefficients), the bound including a sign bit and the
    growth in the recursion, and ``nlimbs1 + nlimbs2 + 1`` if an operand
    may get extended in the recursion. Passing ``WORD_MAX`` as the
    cutoff gives the classical bound.

.. function:: slong _flint_mpn_poly_mulmid_cutoff(slong nlimbs1, slong nlimbs2, int squaring)

    The Karatsuba cutoff used by :func:`_flint_mpn_poly_mulmid` for coefficients
    of the given numbers of limbs: the tuned value
    :func:`_flint_mpn_poly_karatsuba_cutoff`, unless overridden with the
    global ``_flint_mpn_poly_mulmid_force_cutoff`` (for tests and profiling;
    the representation chosen by :func:`_flint_mpn_poly_mulmid_params` can
    likewise be overridden with ``_flint_mpn_poly_mulmid_force_method``).

.. function:: void _flint_mpn_poly_mulmid(nn_ptr res, nn_srcptr f, slong flen, slong nlimbs1, slong norm1, nn_srcptr g, slong glen, slong nlimbs2, slong norm2, slong nlo, slong nhi, slong slimbs, int method)
              void _flint_mpn_poly_mul(nn_ptr res, nn_srcptr f, slong flen, slong nlimbs1, slong norm1, nn_srcptr g, slong glen, slong nlimbs2, slong norm2, slong slimbs, int method)

    Sets ``{res, (nhi - nlo) slimbs}`` to the coefficients `[nlo, nhi)`
    of the product (respectively the full product) of ``{f, flen}`` and
    ``{g, glen}``, packed as described by ``method`` (a biased operand is
    treated as unsigned), choosing automatically between the classical
    algorithm (:func:`_flint_mpn_poly_mulmid_classical`, always used for
    sign-magnitude coefficients) and :func:`_flint_mpn_poly_mulmid_karatsuba`
    with the cutoff :func:`_flint_mpn_poly_mulmid_cutoff`, after discarding the
    coefficients of the inputs that cannot contribute to the window.
    ``slimbs`` must be at least the value given by
    :func:`_flint_mpn_poly_mulmid_slimbs` (as computed by
    :func:`_flint_mpn_poly_mulmid_params`). Squaring is detected when
    ``f == g``, ``flen == glen`` and ``nlimbs1 == nlimbs2``.

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

.. function:: void flint_mpn_signed_div2(mp_ptr res, mp_srcptr x, mp_size_t n)

    Sets ``res`` to ``(x, n)`` divided by two, where ``x`` is viewed
    as a signed integer in two's complement form.

.. function:: int flint_mpn_divides(mp_ptr q, mp_srcptr array1, mp_size_t limbs1, mp_srcptr arrayg, mp_size_t limbsg, mp_ptr temp)

    If ``(arrayg, limbsg)`` divides ``(array1, limbs1)`` then
    ``(q, limbs1 - limbsg + 1)`` is set to the quotient and 1 is 
    returned, otherwise 0 is returned. The temporary space ``temp``
    must have space for ``limbsg`` limbs.

    Assumes ``limbs1 >= limbsg > 0``.

.. function:: void flint_mpn_tdiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              void flint_mpn_tdiv_q(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              void flint_mpn_tdiv_r(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

    Truncating (Euclidean) division: sets `(q, an - bn + 1)` to
    `\lfloor a / b \rfloor` and `(r, bn)` to `a - qb`. Requires
    `an \ge bn \ge 1` and `b_{bn-1} \ne 0`; the top limb of `a` may be zero.
    No aliasing between the output and input arrays is permitted.

    For small and medium operands this wraps GMP's ``mpn_tdiv_qr``.
    When both `bn` and `an - bn + 1` exceed ``FLINT_MPN_TDIV_QR_NEWTON_CUTOFF``
    the quotient is computed by Karp-Markstein Newton division
    (:func:`fixed_div_newton`) using FLINT's multiplication, and unbalanced
    divisions (`an \ge 3 bn` for divisors above
    ``FLINT_MPN_TDIV_QR_UNBALANCED3_CUTOFF`` limbs, `an \ge 4 bn` above
    ``FLINT_MPN_TDIV_QR_UNBALANCED4_CUTOFF`` limbs) are done as a sequence
    of `2bn \times bn` block divisions sharing one approximate inverse of `b`
    (:func:`fixed_inv_newton`). Short divisors with long dividends
    (`4 \le bn < 64` and `an \ge 32 bn`, or `32 \le bn < 64` and
    `an \ge 4 bn`) go through :func:`flint_mpn_preinvn` and
    :func:`flint_mpn_divrem_preinvn` after normalising the divisor, which
    beats GMP by up to 40% in that regime.

.. function:: void _flint_mpn_tdiv_qr_newton(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              void _flint_mpn_tdiv_qr_unbalanced(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              void _flint_mpn_tdiv_qr_preinv(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_srcptr binv, mp_size_t binvn)
              void _flint_mpn_tdiv_qr_preinvn(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              void _flint_mpn_tdiv_qr_gmp(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              void _flint_mpn_tdiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

    The individual algorithms behind :func:`flint_mpn_tdiv_qr`; in all of
    them `r` may be ``NULL``. The Newton version views `a` and `b` as fixed-point
    numbers and computes `n + 2` fraction limbs of `a/b`, where `n = an - bn + 1`;
    since the error is below `4 B^{-2}` at the integer scale, the integer part is
    certified whenever the first fraction limb lies in `[2, B-2]`, in which
    case the remainder follows from a low product; otherwise the candidate
    quotient is corrected by `O(1)` steps using a full product.
    The preinv version requires `(binv, binvn + 2)` to be the output of
    ``fixed_inv_newton(binv, b, bn, binvn)`` with `binvn \ge n + 2` and
    `an \ge n + 2`; the unbalanced version requires `an > 2 bn` and
    `bn \ge 3`.

.. function:: void flint_mpn_cdiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              void flint_mpn_cdiv_q(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              void flint_mpn_cdiv_r(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

    Ceiling division: sets `(q, an - bn + 2)` to `\lceil a / b \rceil` and
    `(r, bn)` to `qb - a \in [0, b)`. Unlike truncating division, the
    quotient can overflow `an - bn + 1` limbs (e.g. `a = B^2 - 1`, `b = B`),
    so one extra output limb, which is always 0 or 1, is written.
    Requirements as for :func:`flint_mpn_tdiv_qr`.

.. function:: int flint_mpn_ndiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              void flint_mpn_ndiv_q(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              int flint_mpn_ndiv_r(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

    Division with rounding to nearest, ties to even: sets `(q, an - bn + 2)`
    to the integer nearest `a/b` and `(r, bn)` to `|a - qb| \le b/2`, and
    returns the sign (`-1`, `0` or `1`) of `a - qb`. As for ceiling division,
    one extra quotient limb (0 or 1) is written.
    Requirements as for :func:`flint_mpn_tdiv_qr`.

.. function:: int flint_mpn_div(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

    If `b` divides `a`, sets `(q, an - bn + 1)` to the quotient and returns 1;
    otherwise returns 0, in which case the contents of `q` are undefined.
    Requires `an \ge bn \ge 1` and `b_{bn-1} \ne 0`. This is a version of
    :func:`flint_mpn_divides` using :func:`flint_mpn_tdiv_qr`.

.. function:: void flint_mpn_divexact(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
              void _flint_mpn_divexact_hensel(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

    Exact division: sets `(q, an - bn + 1)` to `a / b`, assuming that `b`
    divides `a`. Requires `an \ge bn \ge 1` and `b_{bn-1} \ne 0`.

    Single-limb divisors are handled by ``mpn_divexact_1``. When GMP's
    internal ``mpn_divexact`` is available (``FLINT_HAVE_NATIVE_mpn_divexact``)
    it is used when `\min(bn, an - bn + 1)` is below
    ``FLINT_MPN_DIVEXACT_NEWTON_CUTOFF``, unless the quotient is at least four
    times longer than a divisor of at least
    ``FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF`` limbs. Otherwise
    the division is done 2-adically: writing `b = 2^v B^k b'` with `b'` odd,
    the quotient is `(a / (2^v B^k)) \, b'^{-1} \bmod B^{an - bn + 1}`, computed
    with :func:`flint_mpn_bdiv_q`. (A bidirectional variant computing the
    high half of the quotient by Euclidean division was found to be slower
    and is not used.)

.. type:: flint_mpn_divexact_preinv_struct
          flint_mpn_divexact_preinv_t

    Precomputed data for repeated exact division by a fixed divisor: its
    odd part `b'` (with `b = 2^v B^k b'`) and the 2-adic inverse
    `b'^{-1} \bmod B^{bn'+1}`.

.. function:: void flint_mpn_divexact_preinv_init(flint_mpn_divexact_preinv_t pre, mp_srcptr b, mp_size_t bn)
              void flint_mpn_divexact_preinv_clear(flint_mpn_divexact_preinv_t pre)
              void flint_mpn_divexact_preinv(mp_ptr q, mp_srcptr a, mp_size_t an, const flint_mpn_divexact_preinv_t pre)

    Exact division by a fixed divisor `(b, bn)` with a precomputed inverse:
    sets `(q, an - bn + 1)` to `a / b`, assuming `b` divides `a`. Quotients
    of at most `bn' + 1` limbs cost a single :func:`flint_mpn_mullow_n`;
    longer ones use the block Hensel division with the stored inverse (one
    low and one high `bn' \times bn'` product per block). Only the low
    limbs of `a` that can influence the quotient are read (and shifted when
    `b` is even). Compared with a fresh ``mpn_divexact`` for every division,
    which recomputes an inverse each time, this is 1.5-2.5 times faster for
    balanced operands and 10-25% faster for long quotients; divisors of
    2-4 limbs with long quotients are handed to GMP. Used by
    :func:`_fmpz_vec_scalar_divexact_fmpz` and
    :func:`fmpz_mat_scalar_divexact_fmpz` when at least two entries are as
    large as the divisor.

.. function:: int flint_mpn_divisible(mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

    Returns 1 if `b` divides `a` and 0 otherwise. Requires `bn \ge 1` and
    `b_{bn-1} \ne 0`; `a` may have zero top limbs and `an` may be zero.

    After the trivial cases and the 2-adic part (`b = 2^v B^k b'` with `b'`
    odd), single-limb divisors use the Hensel remainder
    (``mpn_modexact_1_odd``), two-limb divisors of dividends of at most four
    limbs use inline Hensel steps with a limb inverse, and long dividends
    are first screened by trial division: the residues of `b` and `a`
    modulo `2^{48} - 1` (nine small prime factors, computed at a quarter of
    a nanosecond per limb by :func:`flint_mpn_mod_2exp48m1`) reveal a prime
    dividing `b` but not `a` for about 60% of random pairs, rejecting them
    in `O(an)` time. The remaining cases go to GMP's ``mpn_divisible_p``
    below ``FLINT_MPN_DIVISIBLE_GMP_CUTOFF`` divisor limbs when available,
    and otherwise to the Hensel division with remainder
    :func:`flint_mpn_bdiv_qr` with `an - bn + 1` quotient limbs, `b'`
    dividing `a` iff the remainder vanishes. Used by :func:`fmpz_divisible`.

.. function:: mp_limb_t flint_mpn_mod_2exp48m1(mp_srcptr a, mp_size_t n)

    Returns a value below `2^{50}` congruent to `a` modulo `2^{48} - 1`
    (the same quantity as GMP's ``mpn_mod_34lsub1``, which is called for
    inputs shorter than 48 limbs when available). Longer inputs use an AVX2
    kernel when compiled with AVX2 support, running at about half the time
    of GMP's assembly, or a portable add/adc carry-chain loop. Only
    available for 64-bit limbs.

.. function:: mp_size_t flint_mpn_pow_bound_limbs(mp_srcptr x, mp_size_t xn, ulong e)

    Returns an upper bound on the number of limbs needed to hold `x^e`
    with :func:`flint_mpn_pow`, including one limb of slack for the
    intermediate products. Requires `x_{xn-1} \ne 0`. For small results the
    bound `e \cdot \mathrm{bits}(x)` is used; otherwise `e \log_2 x`
    is evaluated in double precision from the top two limbs and rounded
    up, which is essentially tight also for small bases with large
    exponents.

.. function:: mp_size_t flint_mpn_pow(mp_ptr res, mp_srcptr x, mp_size_t xn, ulong e)

    Sets `res` to `x^e` by left-to-right binary exponentiation and returns
    the number of limbs of the result. `res` must have room for
    :func:`flint_mpn_pow_bound_limbs` limbs. Requires `xn \ge 1` and
    `x_{xn-1} \ne 0`; `e = 0` gives 1.

    Low zero limbs of `x` are stripped and restored as an offset in the
    output; trailing zero bits are shifted out only when the base is short
    (at most 8 limbs) and restored with one shift of the result. While the
    running power fits in one or two limbs, fixed-size limb arithmetic is
    used. Afterwards the running power alternates between `res` and a
    scratch buffer without copying: the parity of the remaining
    out-of-place operations is chosen so that the final product lands in
    `res`, and multiplications by a single-limb base are done in place with
    ``mpn_mul_1``. Compared with ``mpz_pow_ui`` this is at parity for
    medium sizes and about twice as fast for large bases (where
    :func:`flint_mpn_sqr` uses ``fft_small``). Used by :func:`fmpz_pow_ui`.

.. function:: void flint_mpn_inv(mp_ptr q, mp_srcptr x, mp_size_t xn, mp_size_t n)

    Correctly truncated reciprocal: sets `(q, n - xn + 2)` to
    `\lfloor B^n / x \rfloor`. Requires `x_{xn-1} \ne 0` and `n \ge xn`. The
    top output limb is nonzero only when `x` is a power of `B`.
    For large operands this uses :func:`fixed_inv_newton` with three extra
    fraction limbs, with the same certification and correction scheme as
    :func:`flint_mpn_tdiv_qr`, and avoids forming the numerator `B^n`
    altogether in the common case.

mpz-like interface
--------------------------------------------------------------------------------

.. function:: void flint_mpz_tdiv_qr(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_tdiv_q(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_tdiv_r(mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_fdiv_qr(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_fdiv_q(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_fdiv_r(mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_cdiv_qr(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_cdiv_q(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_cdiv_r(mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_mod(mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_divexact(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
              void flint_mpz_sqrtrem(mpz_ptr s, mpz_ptr r, mpz_srcptr a)
              void flint_mpz_sqrt(mpz_ptr s, mpz_srcptr a)

    Drop-in replacements for the GMP functions of the same names (identical
    semantics for signs, rounding modes and aliasing), with the limb-level
    work done by :func:`flint_mpn_tdiv_qr`, :func:`flint_mpn_divexact` and
    :func:`flint_mpn_sqrtrem`. Divisors (respectively square root
    arguments) of at most two limbs are passed straight to GMP, which has
    dedicated fast paths for them; otherwise the flint_mpn routines choose
    between GMP's mpn layer and Newton iteration. These are used by the
    ``fmpz`` division and square root functions for large operands.

Hensel (2-adic) division and square root
--------------------------------------------------------------------------------

The following functions perform arithmetic modulo `B^n`, developing results
from the least significant limb upwards. They are ports of the
``radix_invmod_bn``, ``radix_divmod_bn``, ``radix_rsqrtmod_bn`` and
``radix_sqrtmod_bn`` algorithms of the :ref:`radix <radix>` module to the machine
word radix `B = 2^{64}`, where the reductions modulo the digit radix
disappear and only the `p = 2` branches of the square root algorithms
survive. Output arrays may not alias input arrays unless stated otherwise.

.. function:: void _flint_mpn_mulhigh_known_low(mp_ptr out, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn, mp_srcptr kl, mp_size_t kl_len, mp_size_t klo, mp_size_t khi, mp_ptr scratch)

    Sets `(out, khi - klo)` to limbs `[klo, khi)` of the product `xy`, given
    that its low `klo` limbs are known and equal to `(kl, kl_len)` (zero
    above ``kl_len``). Since the windowed middle product is a lower
    approximation whose deficit is below `B^2`, the high limbs can be
    recovered from a middle product with three guard limbs, the known low
    limbs telling whether a carry was lost. ``scratch`` must have room for
    `khi` limbs.

.. function:: void flint_mpn_binv(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)

    Sets `(res, n)` to `x^{-1} \bmod B^n`. Requires `x` to be odd.
    Newton iteration on middle products; up to
    ``FLINT_MPN_BINV_MULLOW_CUTOFF`` limbs full low products
    (:func:`flint_mpn_mullow_n`) are used instead, and the cases `n \le 3`
    are inlined.

.. function:: void flint_mpn_bdiv_qr_1(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_limb_t b, mp_size_t n)
              void flint_mpn_bdiv_qr_classical(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t n)
              void flint_mpn_bdiv_qr_karp_markstein(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t n)
              void flint_mpn_bdiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t n)
              void flint_mpn_bdiv_q(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t n)

    Hensel division: given `(a, an)` and odd `(b, bn)`, develop `n` limbs of
    the 2-adic quotient `q` with `qb = a \bmod B^n`. If `r` is not ``NULL``
    (``bdiv_q`` passes ``NULL``), the remainder `(a - qb) / B^n \bmod B^{bn}`
    is written to `(r, bn)`, so that `a = qb + B^n r \pmod{B^{n + bn}}`; the
    division is exact with an `n`-limb quotient iff `r = 0` and `a` has no
    nonzero limbs at positions `\ge n + bn`. `q` may alias `a`.

    Unlike GMP's ``mpn_bdiv_qr``, the number of quotient limbs `n` is an
    explicit parameter independent of `an` and `bn`.

    The classical version runs from the bottom up in blocks of `bn` limbs,
    each costing one low and one high `bn \times bn` product, and is used by
    the dispatcher when `bn < n/2`; the Karp-Markstein version uses a
    half-precision inverse (:func:`flint_mpn_binv`) and one refinement step.

.. function:: int flint_mpn_brsqrt(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)
              int flint_mpn_bsqrt(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)

    If `x \equiv 1 \pmod 8`, sets `(res, n)` to the 2-adic reciprocal
    square root `y` (with `x y^2 \equiv 1 \pmod{B^n}`) respectively square root
    `s` (with `s^2 \equiv x \pmod{B^n}`) that is congruent to 1 modulo 4,
    and returns 1. Otherwise returns 0 without writing anything (an odd
    `x` is a square modulo `B^n` iff `x \equiv 1 \pmod 8`).

    The reciprocal square root is computed by a full-recompute Newton
    iteration in bit precision (each step maps `c` correct bits to
    `2c - 2`), using known-low high products above the basecase range and
    full low products below ``FLINT_MPN_BRSQRT_MULLOW_CUTOFF`` limbs; the
    square root uses one Karp-Markstein refinement on top of a reciprocal
    square root with `\lceil n/2 \rceil + 1` limbs.

Square root
--------------------------------------------------------------------------------

.. function:: mp_size_t flint_mpn_sqrtrem(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
              void _flint_mpn_sqrtrem_newton(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
              void _flint_mpn_sqrtrem_gmp(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)

    Sets `(s, \lceil an/2 \rceil)` to `\lfloor \sqrt{a} \rfloor`. If `r` is
    not ``NULL``, `(r, \lceil an/2 \rceil + 1)` is set to the remainder
    `a - s^2 \le 2s` (zero-padded) and the number of limbs of the remainder
    is returned; `r` must have room for `\max(an, 2)` limbs, since GMP's
    ``mpn_sqrtrem`` writes the remainder in place using `an` limbs. If `r` is ``NULL`` the remainder is not formed and the
    return value is 0 if `a` is a perfect square and 1 otherwise (GMP's
    convention). Requires `an \ge 1` and `a_{an-1} \ne 0`.

    Two-limb inputs use the hardware double-precision square root for the
    initial approximation, followed by one Newton step with a 128/64-bit
    division when `a \ge 2^{100}` and a final adjustment; this is about
    twice as fast as GMP below `2^{100}` and 1.3-1.4 times faster above.
    Below ``FLINT_MPN_SQRTREM_NEWTON_CUTOFF`` input limbs the function
    otherwise wraps GMP's ``mpn_sqrtrem`` (which needs `an` limbs of
    remainder space; the remainder is copied). Above it, `a` is viewed as a
    fixed-point number in `[B^{-2}, 1)` and :func:`fixed_sqrt_newton` is used
    with three extra fraction limbs, so that the truncated root is certified
    whenever the first fraction limb lies in `[2, B-2]`; otherwise it is
    corrected by `O(1)` steps. With a ``NULL`` remainder the Newton path
    tests `a = s^2` on the low two limbs before comparing in full.

.. function:: int flint_mpn_sqrt(mp_ptr s, mp_srcptr a, mp_size_t an)

    Checked exact square root, the analogue of :func:`flint_mpn_div`: if
    `a` is a perfect square, sets `(s, \lceil an/2 \rceil)` to its square
    root and returns 1; otherwise returns 0, in which case the contents of
    `s` are undefined. Nonsquares are usually rejected before any root is
    computed by quadratic residue tests modulo 256 and modulo the factors
    9, 5, 7, 13, 17, 97, 241 and 257 of `2^{48} - 1`, which are all obtained
    from one reduction modulo `2^{48} - 1` (about one random nonsquare in
    900 survives). That reduction uses an AVX2 kernel when available, which
    runs in about half the time of GMP's ``mpn_mod_34lsub1`` from 128 limbs
    on and is otherwise equivalent to it; GMP's routine is used for short
    inputs when it is available (``FLINT_HAVE_NATIVE_mpn_mod_34lsub1``), and
    a portable add/adc carry-chain version exists as well. On 32-bit builds
    the residues come from two ``mpn_mod_1`` calls instead.

.. function:: int flint_mpn_is_square(mp_srcptr a, mp_size_t an)

    Returns 1 if `a` is a perfect square and 0 otherwise, using the same
    screening as :func:`flint_mpn_sqrt` and only allocating space for the
    root when the screening passes. Used by :func:`fmpz_is_square`; random
    nonsquares are rejected in about the time of GMP's
    ``mpz_perfect_square_p`` (faster from a few thousand bits), and large
    squares are certified with the Newton square root.

Division and modular arithmetic with precomputed inverses
--------------------------------------------------------------------------------

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

.. function:: mp_limb_t flint_mpn_divrem_1_preinv(mp_ptr q, mp_srcptr a, mp_size_t n, mp_limb_t d, mp_limb_t dinv, unsigned int norm)

    Divide ``a, n`` by the limb ``d``, writing the quotient to ``q, n``
    and returning the remainder. Requires ``n`` and ``d`` to be positive.
    Allows ``a`` and ``q`` to be aliased. Requires a single-limb inverse ``dinv``
    precomputed by :func:`n_preinvert_limb` and the
    number of leading zero bits of ``d`` as ``norm``.

    This is equivalent to ``mpn_divrem_1(q, 0, a, n, d)`` but faster for small
    ``n``. Typically ``mpn_divrem_1`` will be faster for large ``n`` as it
    has dedicated assembly code on many architectures whereas
    ``flint_mpn_divrem_1_preinv`` currently does not.

.. function:: mp_limb_t flint_mpn_divrem_2_1_preinv_norm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv)
              mp_limb_t flint_mpn_divrem_2_1_preinv_unnorm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv, unsigned int norm)
              mp_limb_t flint_mpn_divrem_3_1_preinv_norm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv)
              mp_limb_t flint_mpn_divrem_3_1_preinv_unnorm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv, unsigned int norm)

    Versions of :func:`flint_mpn_divrem_1_preinv` specialized for length 2 and 3.
    The ``_norm`` functions require a normalised divisor while the ``_unnorm``
    functions require an unnormalised divisor with positive ``norm``.

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

    We require `a` and `b` to be reduced modulo `d` before calling the
    function. 

.. function:: void flint_mpn_mulmod_preinvn_2(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_srcptr d, mp_srcptr dinv, ulong norm)

    Version of :func:`flint_mpn_mulmod_preinv1` specialized for two limbs.
    The behavior is not exactly the same: `a` and `b` are assumed to
    be unshifted, and the output is unshifted.

.. function:: void flint_mpn_fmmamod_preinvn(mp_ptr r, mp_srcptr a1, mp_srcptr b1, mp_srcptr a2, mp_srcptr b2, mp_size_t n, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)
              void flint_mpn_fmmamod_preinvn_2(mp_ptr r, mp_srcptr a1, mp_srcptr b1, mp_srcptr a2, mp_srcptr b2, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)

    Given ``dnormed`` containing a normalised integer `d 2^{norm}` with precomputed inverse ``dinv``
    provided by ``flint_mpn_preinvn``, computes `a_1 b_1 + a_2 b_2 \pmod{d}`. We require
    all operands to be reduced modulo `d`.

Preconditioned modular multiplication
--------------------------------------------------------------------------------

Currently two algorithms are implemented for preconditioned multiplication:
Shoup multiplication and the matrix algorithm. An FFT variant may be added in the future.

.. function:: int flint_mpn_mulmod_want_precond(mp_size_t n, slong num, ulong norm)

    Assuming a precision of `n` limbs and that one wants to perform `num`
    multiplications with a fixed (preconditioned) operand with norm ``norm``,
    return one of the following constants indicating
    which algorithm is better (accounting for the cost of pretransforming
    the operand).

    * ``MPN_MULMOD_PRECOND_NONE`` - should use :func:`flint_mpn_mulmod_preinvn` (no precomputation)

    * ``MPN_MULMOD_PRECOND_SHOUP`` - should use :func:`flint_mpn_mulmod_precond_shoup`

    * ``MPN_MULMOD_PRECOND_MATRIX`` - should use :func:`flint_mpn_mulmod_precond_matrix`

.. function:: void flint_mpn_mulmod_precond_shoup(mp_ptr res, mp_srcptr a, mp_srcptr apre, mp_srcptr b, mp_size_t n, mp_srcptr d, ulong norm)

    Compute `ab \pmod{d}` given precomputed data for ``apre``
    generated with :func:`flint_mpn_mulmod_precond_shoup_precompute`.
    We require that `b` is reduced modulo `d`.

.. function:: void flint_mpn_mulmod_precond_shoup_precompute(mp_ptr apre, mp_srcptr a, mp_size_t n, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)

    Given `0 \le a < d`, precompute data for
    multiplication by `a` modulo `d` using Shoup's method.
    The modulus is given as ``dnormed`` containing `d 2^{norm}` together
    with precomputed inverse ``dinv``.
    The destination ``apre`` must have space for `n` limbs.

.. function:: void flint_mpn_mulmod_precond_matrix(mp_ptr rp, mp_srcptr apre, mp_srcptr b, mp_size_t n, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)

    Given ``dnormed`` containing a normalised integer `d 2^{norm}` with precomputed inverse ``dinv``
    provided by :func:`flint_mpn_preinvn`, computes `ab \pmod{d}`. We require
    `b` to be reduced modulo `d`.
    The user provides the operand `a` via the ``apre`` argument in the
    pretransformed representation returned by :func:`flint_mpn_mulmod_precond_matrix_precompute`.
    The complexity of this function is `O(n^2)`. Requires `n \ge 2`.

.. function:: void flint_mpn_mulmod_precond_matrix_precompute(mp_ptr apre, mp_srcptr a, mp_size_t n, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)

    Given ``dnormed`` containing a normalised integer `d 2^{norm}` with precomputed inverse ``dinv``
    and an integer `a` which is reduced modulo `d`,
    write to ``apre`` a pretransformed representation of `a`
    for use with :func:`flint_mpn_mulmod_precond_matrix`.
    Currently, the output consists of `n \times n` limbs storing
    `a 2^{norm} \beta^i \mod {d 2^{norm}}` for `0 \le i < n` where `\beta` is the limb
    radix, plus one junk limb.

.. function:: mp_size_t flint_mpn_mulmod_precond_matrix_alloc(mp_size_t n)

    The *alloc* function returns the number of limbs of space required for
    :func:`flint_mpn_mulmod_precond_matrix_precompute`
    given a modulus with `n` limbs.

.. function:: void flint_mpn_fmmamod_precond_matrix(mp_ptr rp, mp_srcptr a1pre, mp_srcptr b1, mp_srcptr a2pre, mp_srcptr b2, mp_size_t n, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)

    Analogous to :func:`flint_mpn_mulmod_precond_matrix`, but computes `a_1 b_1 + a_2 b_2` modulo `d`.


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

.. function:: void flint_mpn_urandomb(mp_ptr rp, flint_rand_t state, flint_bitcnt_t n)

    Generates a uniform random number of ``n`` bits and stores
    it on ``rp``.

.. function:: void flint_mpn_urandomm(mp_ptr rp, flint_rand_t state, mp_srcptr xp, mp_size_t xn)

    Generates a uniform random number between 0 inclusive and ``(xp, xn)``
    exclusive`[0, x)` and stores it on ``rp``. The most significant limb of
    ``xp`` is required to be nonzero. This function will write ``xn`` limbs to
    ``rp`` even if the largest possible value has one fewer limb.

.. function:: void flint_mpn_rrandom(mp_ptr rp, flint_rand_t state, mp_size_t n)

    Generates a random number with ``n`` limbs and stores
    it on ``rp``. The number it generates will tend to have
    long strings of zeros and ones in the binary representation.

    Useful for testing functions and algorithms, since this kind of random
    numbers have proven to be more likely to trigger corner-case bugs.

.. function:: void flint_mpn_rrandomb(mp_ptr rp, flint_rand_t state, flint_bitcnt_t nbits)

    Generates a random number with ``nbits`` bits and stores
    it on ``rp``. The number it generates will tend to have
    long strings of zeros and ones in the binary representation.



Complex multiplication
--------------------------------------------------------------------------------

Multiplication of Gaussian integers represented as pairs of limb arrays
with separate sign bits (0 meaning nonnegative). No aliasing is permitted
between output and input arrays.

.. function:: void flint_mpn_mul_complex(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
              void flint_mpn_sqr_complex(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)

    Sets `zr + zi i = (ar + ai i)(br + bi i)` (respectively the square of
    `ar + ai i`). Each part takes an independent length, at least 1 limb
    and not necessarily normalized. A *signed length* is written for each
    output: the magnitude occupies ``|*zr_len|`` limbs and a negative
    value means the result is negative; nothing above ``|*zr_len|`` limbs
    is written. The outputs must have room for
    ``max(arn, ain) + max(brn, bin) + 1`` limbs (``2 max(arn, ain) + 1``
    for the square). The algorithm is selected from the shape: schoolbook
    when a part is much shorter than its partner, Karatsuba when the
    parts are internally balanced, and a transformed (fft_small) method
    for large balanced operands.

.. function:: void flint_mpn_mul_complex_classical(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
              void flint_mpn_mul_complex_karatsuba(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
              int flint_mpn_mul_complex_fft_small(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
              void flint_mpn_sqr_complex_classical(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)
              void flint_mpn_sqr_complex_karatsuba(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)
              int flint_mpn_sqr_complex_fft_small(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)

    The individual algorithms behind the general functions, exposed for
    comparison and tuning. All accept any shape. The *fft_small* variants
    return 0, leaving the outputs untouched, when the method is
    unavailable or refuses the operands.

.. function:: void flint_mpn_mulhigh_n_complex(nn_ptr zr, int * zr_sgn, nn_ptr zi, int * zi_sgn, nn_srcptr ar, int ar_sgn, nn_srcptr ai, int ai_sgn, nn_srcptr br, int br_sgn, nn_srcptr bi, int bi_sgn, mp_size_t n)
              void flint_mpn_sqrhigh_n_complex(nn_ptr zr, int * zr_sgn, nn_ptr zi, int * zi_sgn, nn_srcptr ar, int ar_sgn, nn_srcptr ai, int ai_sgn, mp_size_t n)

    High products: all parts share the length `n`, and each output
    receives exactly `n + 1` limbs, zero padded, plus a sign -- the limbs
    `[n, 2n]` of the exact result. Relative to the exact value the error
    is below `2 + 3(n + 4)/2^{64}` ulp of the lowest returned limb for
    the product and below `2 + 2(n + 4)/2^{64}` for the square -- each
    underlying :func:`flint_mpn_mulhigh_n`, read as `n` limbs, errs by
    `(-1 - \varepsilon, +\varepsilon)` ulp against the exact value with
    `\varepsilon = (n + 4)/2^{64}`, and each output combines at most
    three -- so below 3 ulp for any practical `n`. The transformed path
    stays within `(-1.5, +0.5)` ulp.

.. var:: slong flint_mpn_mul_complex_fft_cutoff
         slong flint_mpn_sqr_complex_fft_cutoff

    Sizes in limbs from which the general functions use the transformed
    path. Machine dependent; writable for tuning.

Multimodular reduction and Chinese remaindering
--------------------------------------------------------------------------------

The following functions perform simultaneous reduction of a multiprecision
integer modulo a fixed list of single-limb moduli (multi mod) and the
inverse reconstruction (multi CRT). All the work is done with
mpn arithmetic and precomputed data stored in a :type:`flint_mpn_crt_t`.
The :type:`fmpz_comb_t` in the ``fmpz`` module is a thin wrapper around
this structure.

The moduli must be nonzero and (for Chinese remaindering) pairwise
coprime, but need not be prime.

.. type:: flint_mpn_crt_struct

.. type:: flint_mpn_crt_t

    Precomputed data for multi mod / multi CRT with respect to a list of
    single-limb moduli `m_0, \ldots, m_{n-1}` with product `P`.

    Internally, consecutive tiny moduli are batched into single-limb
    products (leaves), which are grouped into chunks whose products are
    the leaves of a balanced subproduct tree. Each level of the tree is
    stored as a contiguously packed array, and each node above a
    threshold is accompanied by a precomputed inverse for fast division.
    Modular reduction descends the tree using :func:`flint_mpn_mod_preinvn`
    and finishes with a basecase using dot products with precomputed powers
    of `2^{\mathtt{FLINT\_BITS}}` modulo the leaves; Chinese remaindering
    starts with a basecase using :func:`mpn_addmul_1` with precomputed
    multipliers (with the fractional cofactors of the tree folded in),
    ascends the tree without intermediate reductions, and performs a single
    reduction modulo `P` at the end. Values that are much smaller than `P`
    are detected early (by reconstructing modulo the product of the first
    few moduli and verifying against all residues) and returned without
    traversing the tree. When the whole product is small, the
    fixed-length templates from ``crt_helpers.h`` are used.

    Precomputation costs `O(M(N) \log N)` operations for `N` total bits,
    using a remainder tree for the cofactors (no large modular inverses),
    and `O(N \log N)` memory.

    The following fields are public: ``num_primes``, ``primes``
    (a copy of the moduli), ``prod`` and ``prod_len`` (the product `P` as an
    mpn integer), and ``tmp_limbs`` (the size of the workspace required
    by the conversion functions).

.. macro:: FLINT_MPN_CRT_MOD
           FLINT_MPN_CRT_CRT

    Flags selecting the operations to precompute for.

.. function:: void flint_mpn_crt_init(flint_mpn_crt_t C, nn_srcptr primes, slong num_primes)
              void flint_mpn_crt_init2(flint_mpn_crt_t C, nn_srcptr primes, slong num_primes, int flags)

    Initialises *C* for the given moduli. The version with *flags*
    (a bitwise or of ``FLINT_MPN_CRT_MOD`` and ``FLINT_MPN_CRT_CRT``) only
    performs the precomputations needed for the selected operations,
    which saves time and memory when only one direction is needed;
    :func:`flint_mpn_crt_init` selects both.
    Throws an exception if the moduli are zero, or if Chinese
    remaindering is selected and the moduli are not pairwise coprime.

.. function:: void flint_mpn_crt_init_tuned(flint_mpn_crt_t C, nn_srcptr primes, slong num_primes, int flags, slong crt_chunk_bits, slong mod_base_bits, slong preinv_cutoff)

    Like :func:`flint_mpn_crt_init2`, but with explicit tuning parameters:
    the approximate number of bits in a CRT basecase chunk, the number of
    bits at which modular reduction switches to the basecase, and the
    number of limbs above which precomputed inverses are used for
    division. This is mainly intended for profiling.

.. function:: void flint_mpn_crt_clear(flint_mpn_crt_t C)

    Frees the memory allocated by *C*.

.. function:: void flint_mpn_multi_mod(nn_ptr out, nn_srcptr x, slong xn, const flint_mpn_crt_t C, nn_ptr tmp)

    Reduces the nonnegative integer `(x, xn)` modulo all the moduli,
    writing the residues to ``out[0], ..., out[num_primes - 1]``.
    The input may have any size, including zero limbs and non-normalised
    top limbs. The workspace *tmp* must have space for ``C->tmp_limbs``
    limbs; alternatively, *tmp* may be *NULL*, in which case the
    workspace is allocated internally.

.. function:: int flint_mpn_multi_crt(nn_ptr out, nn_srcptr res, const flint_mpn_crt_t C, int sign, nn_ptr tmp)

    Reconstructs the integer `x` from the residues
    ``res[0], ..., res[num_primes - 1]``, each of which must be reduced
    modulo the corresponding modulus. The result is written as
    ``C->prod_len`` limbs (zero padded) to *out*.

    If *sign* is zero, the result is the unique `x` with `0 \le x < P`,
    and the return value is zero.
    If *sign* is nonzero, the result is the unique `x` with `-P < 2x \le P`;
    its absolute value is written to *out* and the return value indicates
    whether `x` is negative.

    The workspace *tmp* is as for :func:`flint_mpn_multi_mod`.

.. function:: void flint_mpn_multi_mod_vec(nn_ptr out, slong out_stride, nn_srcptr x, slong xn, slong len, const flint_mpn_crt_t C, nn_ptr tmp)
              void flint_mpn_multi_crt_vec(nn_ptr out, slong out_stride, int * negative, nn_srcptr res, slong res_stride, slong len, const flint_mpn_crt_t C, int sign, nn_ptr tmp)

    Vector versions. In the mod version, *x* is a packed array of *len*
    nonnegative integers of *xn* limbs each and the residue of entry `i`
    modulo modulus `l` is written to ``out[l * out_stride + i]``
    (so *out_stride* must be at least *len*). In the CRT version, the
    residue of entry `i` modulo modulus `l` is read from
    ``res[l * res_stride + i]``, entry `i` of the output is written to
    ``out + i * out_stride`` (with ``out_stride >= C->prod_len``), and if
    *sign* is nonzero the sign flags are written to ``negative[i]``
    (*negative* may be *NULL* when *sign* is zero).

    These are equivalent to looping over the entries but faster when the
    entries are small compared to the product, since the precomputed
    tables are then traversed once per block of entries instead of once
    per entry.

.. function:: void flint_mpn_multi_mod_once(nn_ptr out, nn_srcptr x, slong xn, nn_srcptr primes, slong num_primes)
              int flint_mpn_multi_crt_once(nn_ptr out, slong * outn, nn_ptr prod, nn_srcptr res, nn_srcptr primes, slong num_primes, int sign)

    One-shot versions of :func:`flint_mpn_multi_mod` and
    :func:`flint_mpn_multi_crt` which do not require a precomputed
    structure. The mod version creates one internally (with mod-only
    data). The CRT version instead traverses the subproduct tree depth
    first, computing subproducts and cofactors on the fly and freeing
    them as soon as possible, so that the memory usage is a small
    constant multiple of the size of the product `P` rather than
    `O(N \log N)`; this matters for very large reconstructions (e.g.
    Bernoulli numbers with billions of bits). It is also somewhat faster
    than initialising a full structure for a single use.
    In the CRT version, *out* must have space
    for *num_primes* limbs; the actual length of the product of the
    moduli is written to *outn*, and if *prod* is not *NULL*, the product
    itself is written there (*prod* must also have space for *num_primes*
    limbs). The return value is as for :func:`flint_mpn_multi_crt`.

    Large instances are parallelised over the nodes of the subproduct
    tree (both in precomputation and in the conversions) when several
    threads are available.

.. function:: ulong flint_mpn_crt_mod_leaf(nn_srcptr a, slong an, const flint_mpn_crt_t C, slong j)

    Reduces the nonnegative integer `(a, an)` modulo the `j`-th leaf
    modulus of *C* (the product of one or more consecutive batched
    moduli; leaf `j` is modulus `j` when all moduli exceed 32 bits).
    Requires ``an <= C->mod_pow_limbs``. This exposes the basecase of
    :func:`flint_mpn_multi_mod` and is useful with a single modulus, as
    a faster replacement for :func:`mpn_mod_1` when many reductions with
    the same modulus are needed.
