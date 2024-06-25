.. _nmod-vec:

**nmod_vec.h** -- vectors over integers mod n (word-size n)
===============================================================================

Memory management
--------------------------------------------------------------------------------


.. function:: nn_ptr _nmod_vec_init(slong len)

    Returns a vector of the given length. The entries are not necessarily
    zero.

.. function:: void _nmod_vec_clear(nn_ptr vec)

    Frees the memory used by the given vector.


Random functions
--------------------------------------------------------------------------------


.. function:: void _nmod_vec_randtest(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod)

    Sets ``vec`` to a random vector of the given length with entries
    reduced modulo ``mod.n``.


Basic manipulation and comparison
--------------------------------------------------------------------------------


.. function:: void _nmod_vec_set(nn_ptr res, nn_srcptr vec, slong len)

    Copies ``len`` entries from the vector ``vec`` to ``res``.

.. function:: void _nmod_vec_zero(nn_ptr vec, slong len)

    Zeros the given vector of the given length.

.. function:: void _nmod_vec_swap(nn_ptr a, nn_ptr b, slong length)

    Swaps the vectors ``a`` and ``b`` of length `n` by actually
    swapping the entries.

.. function:: void _nmod_vec_reduce(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)

    Reduces the entries of ``(vec, len)`` modulo ``mod.n`` and set
    ``res`` to the result.

.. function:: flint_bitcnt_t _nmod_vec_max_bits(nn_srcptr vec, slong len)

    Returns the maximum number of bits of any entry in the vector.

.. function:: int _nmod_vec_equal(nn_srcptr vec, nn_srcptr vec2, slong len)

    Returns `1` if ``(vec, len)`` is equal to ``(vec2, len)``,
    otherwise returns `0`.


Printing
--------------------------------------------------------------------------------


.. function:: void _nmod_vec_print_pretty(nn_srcptr vec, slong len, nmod_t mod)

    Pretty-prints ``vec`` to ``stdout``. A header is printed followed by the
    vector enclosed in brackets. Each entry is right-aligned to the width of
    the modulus written in decimal, and the entries are separated by spaces.
    For example::

        <length-12 integer vector mod 197>
        [ 33 181 107  61  32  11  80 138  34 171  86 156]

.. function:: int _nmod_vec_fprint_pretty(FILE * file, nn_srcptr vec, slong len, nmod_t mod)

    Same as ``_nmod_vec_print_pretty`` but printing to ``file``.

.. function:: int _nmod_vec_print(nn_srcptr vec, slong len, nmod_t mod)

    Currently, same as ``_nmod_vec_print_pretty``.

.. function:: int _nmod_vec_fprint(FILE * f, nn_srcptr vec, slong len, nmod_t mod)

    Currently, same as ``_nmod_vec_fprint_pretty``.


Arithmetic operations
--------------------------------------------------------------------------------


.. function:: void _nmod_vec_add(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)

    Sets ``(res, len)`` to the sum of ``(vec1, len)``
    and ``(vec2, len)``.

.. function:: void _nmod_vec_sub(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)

    Sets ``(res, len)`` to the difference of ``(vec1, len)``
    and ``(vec2, len)``.

.. function:: void _nmod_vec_neg(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)

    Sets ``(res, len)`` to the negation of ``(vec, len)``.

.. function:: void _nmod_vec_scalar_mul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)

    Sets ``(res, len)`` to ``(vec, len)`` multiplied by `c`. The element
    `c` and all elements of ``vec`` are assumed to be less than ``mod.n``.

.. function:: void _nmod_vec_scalar_mul_nmod_shoup(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)

    Sets ``(res, len)`` to ``(vec, len)`` multiplied by `c` using
    :func:`n_mulmod_shoup`. `mod.n` should be less than `2^{\mathtt{FLINT\_BITS} - 1}`. `c`
    and all elements of ``vec`` should be less than ``mod.n``.

.. function:: void _nmod_vec_scalar_addmul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)

    Adds ``(vec, len)`` times `c` to the vector ``(res, len)``. The element
    `c` and all elements of ``vec`` are assumed to be less than ``mod.n``.


Dot products
--------------------------------------------------------------------------------

Dot products functions and macros rely on several implementations, depending on
the length of this dot product and on the underlying modulus. What
implementations will be called is determined via ``_nmod_vec_dot_params``,
which returns a ``dot_params_t`` element which can then be used as input to the
dot product routines.

The efficiency of the different approaches range roughly as follows, from
faster to slower, on 64 bit machines. In all cases, modular reduction is only
performed at the very end of the computation.

- moduli up to `1515531528` (about `2^{30.5}`): implemented via single limb
  integer multiplication, using explicit vectorization if supported (current
  support is for AVX2);

- moduli that are a power of `2` up to `2^{32}`: same efficiency as the above
  case;

- moduli that are a power of `2` between `2^{33}` and `2^{63}`: efficiency
  between that of the above case and that of the below one (depending on the
  machine and on automatic vectorization);

- other moduli up to `2^{32}`: implemented via single limb integer
  multiplication combined with accumulation in two limbs;

- moduli more than `2^{32}`, unreduced dot product fits in two limbs:
  implemented via two limbs integer multiplication, with a final modular
  reduction;

- unreduced dot product fits in three limbs, moduli up to about `2^{62.5}`:
  implemented via two limbs integer multiplication, with intermediate
  accumulation of sub-products in two limbs, and overall accumulation in three
  limbs;

- unreduced dot product fits in three limbs, other moduli: implemented via two
  limbs integer multiplication, with accumulation in three limbs.


.. type:: dot_params_t

.. function:: dot_params_t _nmod_vec_dot_params(slong len, nmod_t mod)

    Returns a ``dot_params_t`` element. This element can be used as input for
    the dot product macros and functions that require it, for any dot product
    of vector with entries reduced modulo ``mod.n`` and whose length is less
    than or equal to ``len``.

    Internals, subject to change: its field ``method`` indicates the method that
    will be used to compute a dot product of this length ``len`` when working
    with the given ``mod``. Its field ``pow2_precomp`` is set to ``2**DOT_SPLIT_BITS
    % mod.n`` if ``method == _DOT2_SPLIT``, and to `0` otherwise.

.. function:: ulong _nmod_vec_dot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)

    Returns the dot product of (``vec1``, ``len``) and (``vec2``, ``len``). The
    input ``params`` has type ``dot_params_t`` and must have been computed via
    ``_nmod_vec_dot_params`` with the specified ``mod`` and with a length
    greater than or equal to ``len``.

.. function:: ulong _nmod_vec_dot_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)

    The same as ``_nmod_vec_dot``, but reverses ``vec2``.

.. function:: ulong _nmod_vec_dot_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod, dot_params_t params)

    Returns the dot product of (``vec1``, ``len``) and the values at
    ``vec2[i][offset]``. The input ``params`` has type ``dot_params_t`` and
    must have been computed via ``_nmod_vec_dot_params`` with the specified
    ``mod`` and with a length greater than or equal to ``len``.

.. macro:: NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, params)

    Effectively performs the computation::

        res = 0;
        for (i = 0; i < len; i++)
            res += (expr1) * (expr2);

    but with the arithmetic performed modulo ``mod``. The input ``params`` has
    type ``dot_params_t`` and must have been computed via
    ``_nmod_vec_dot_params`` with the specified ``mod`` and with a length
    greater than or equal to ``len``.

    ``nmod.h`` has to be included in order for this macro to work (order of
    inclusions does not matter).

.. function:: int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod)

    Returns the number of limbs (0, 1, 2 or 3) needed to represent the
    unreduced dot product of two vectors of length ``len`` having entries
    modulo ``mod.n``, assuming that ``len`` is nonnegative and that
    ``mod.n`` is nonzero. The computed bound is tight. In other words,
    this function returns the precise limb size of ``len`` times
    ``(mod.n - 1)**2``.

.. function:: int _nmod_vec_dot_bound_limbs_from_params(slong len, nmod_t mod, dot_params_t params)

    Same specification as ``_nmod_vec_dot_bound_limbs``, but uses the additional
    input ``params`` to reduce the amount of computations; for correctness
    ``params`` must have been computed for the specified ``len`` and ``mod``.

