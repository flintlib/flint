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

Montgomery reduction and modular exponentiation
--------------------------------------------------------------------------------

.. function:: mp_size_t flint_mpn_mulmod_bnm1_next_size(mp_size_t n)
              mp_size_t flint_mpn_mulmod_bnm1_itch(mp_size_t rn)
              void flint_mpn_mulmod_bnm1(nn_ptr rp, mp_size_t rn, nn_srcptr ap, mp_size_t an, nn_srcptr bp, mp_size_t bn, nn_ptr tp)

    Sets `r` to a representative below `B^{rn}` of `a b \bmod
    (B^{rn} - 1)`, for ``rn`` obtained from ``next_size`` and operands
    of at most ``rn`` limbs, using ``itch(rn)`` limbs of scratch space.
    No aliasing is permitted.

    For even ``rn`` the computation splits by CRT into products modulo
    `B^{h} - 1` (recursively) and `B^{h} + 1`
    (:func:`flint_mpn_mulmod_2expp1_basecase`) with `h = rn/2`,
    stopping at odd or small sizes with a plain multiplication; the
    cost is about `0.6` multiplications of size ``rn``.

.. function:: void _flint_mpn_binvert(nn_ptr v, nn_srcptr m, mp_size_t n)

    Sets `v = m^{-1} \bmod B^n` for odd `m` of `n` limbs, by Hensel
    lifting. No aliasing is permitted.

.. function:: void _flint_mpn_redc_n(nn_ptr t, nn_srcptr X, nn_srcptr m, nn_srcptr minv, mp_size_t n, nn_ptr scratch)

    Montgomery reduction: sets `t = X B^{-n} \bmod m`, canonical in
    `[0, m)`, for odd `m` of `n` limbs and `0 \le X < m B^n` (`2n` limbs),
    where ``minv`` is `-m^{-1} \bmod B^n` (the negation of the output of
    :func:`_flint_mpn_binvert`). Requires `2n` limbs of scratch space.
    `t` must not alias `X`, `m` or ``minv``; `X` is preserved.

    With `q = X_{lo} \cdot minv \bmod B^n` one has `X + q m = t B^n`
    exactly. The high product `\lfloor q m / B^n \rfloor` is obtained
    from :func:`flint_mpn_mulhigh_n`, whose rounding is resolved in
    constant time by comparing the returned guard limb against the top
    limb of the exactly known low half `B^n - X_{lo}`; see the source
    for the argument, [Mul2000]_ and [HanZim2004b]_ for the analysis
    of short products, and [DEGR2024]_ for truncated multiplication
    in the Montgomery setting.

.. function:: void flint_mpn_powm_preinvn(nn_ptr r, nn_srcptr b, nn_srcptr e, mp_size_t en, nn_srcptr m, mp_size_t mn, nn_srcptr dinv, flint_bitcnt_t norm)

    As :func:`flint_mpn_powm`, taking the precomputed inverse ``dinv``
    of `m 2^{norm}` from :func:`flint_mpn_preinvn` with
    ``norm = clz(m[mn-1])``, as stored for instance in an
    ``fmpz_preinvn_t``. Both entry points share a single dispatcher, so
    this one differs only in having an inverse to hand: exponents small
    enough that the division-based basecase wins are run on the supplied
    inverse with no per-call precomputation, and `2`, `3` and `4` skip
    the sliding-window bookkeeping entirely. The cutoff tracks the point
    where :func:`flint_mpn_powm` hands a small modulus to ``mpz_powm``,
    since past that there is nothing for the inverse to save. Larger
    exponents fall through to the Montgomery machinery, whose setup then
    amortizes.

.. function:: void flint_mpn_powm(nn_ptr r, nn_srcptr b, nn_srcptr e, mp_size_t en, nn_srcptr m, mp_size_t mn)

    Sets `r = b^e \bmod m`. The modulus `m` has ``mn`` limbs and must be
    normalised (nonzero top limb); `b` must be reduced modulo `m` and
    occupies ``mn`` limbs; the exponent `e` has ``en`` limbs (``en`` may
    be zero, denoting `e = 0`, in which case `r = 1 \bmod m`). No
    aliasing of `r` with the inputs is permitted.

    Even moduli are handled by CRT between the odd part and the power
    of two; the odd part is computed by the *redc* stage below.
    Single-limb moduli with a single-limb exponent use plain ``nmod``
    arithmetic.

.. function:: void _flint_mpn_powm_basecase(nn_ptr r, nn_srcptr b, nn_srcptr e, mp_size_t en, nn_srcptr m, mp_size_t mn)
              void _flint_mpn_powm_redc(nn_ptr r, nn_srcptr b, nn_srcptr e, mp_size_t en, nn_srcptr m, mp_size_t mn)

    Stage implementations of :func:`flint_mpn_powm`, requiring `m > 1`
    normalised, `b < m` and `e` nonzero with nonzero top limb. Both use
    sliding window exponentiation with single-limb bases handled by
    scalar window multipliers and `O(n)` reductions.

    The *basecase* version works on residues shifted left by the norm
    of `m` via :func:`flint_mpn_mulmod_preinvn` and accepts any `m > 1`.

    The *redc* version requires `m` odd and uses Montgomery
    multiplication. Below a threshold the reduction is
    :func:`_flint_mpn_redc_n`. Above it (with ``FLINT_HAVE_FFT_SMALL``),
    the reduction follows the structure of GMP's ``redc_n``: since the
    low half of `q m` is exactly `B^{mn} - X_{lo}`, the high half
    `H < m` is recovered from the wraparound product
    `q m \bmod (B^{nnc} - 1)`, computed as a cyclic convolution against
    a cached transform of `m` (see
    :func:`fft_small_plan_init_mpn_cyclic`) followed by a limb
    rotation. At still larger sizes the quotient
    `q = X_{lo} \cdot (-m^{-1}) \bmod B^{mn}` also runs in the
    transform domain against a cached transform of the inverse. The
    Montgomery radix stays `B^{mn}` throughout, so the multiplications
    themselves are unpadded calls to ``flint_mpn_sqr``,
    ``flint_mpn_mul_n`` and ``flint_mpn_mullow_n``.

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
