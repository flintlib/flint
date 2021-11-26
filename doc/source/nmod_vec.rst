.. _nmod-vec:

**nmod_vec.h** -- vectors over integers mod n (word-size n)
===============================================================================

Memory management
--------------------------------------------------------------------------------


.. function:: mp_ptr _nmod_vec_init(slong len)

    Returns a vector of the given length. The entries are not necessarily
    zero.

.. function:: void _nmod_vec_clear(mp_ptr vec)

    Frees the memory used by the given vector.


Random functions
--------------------------------------------------------------------------------


.. function:: void _nmod_vec_randtest(mp_ptr vec, flint_rand_t state, slong len, nmod_t mod)

    Sets ``vec`` to a random vector of the given length with entries 
    reduced modulo ``mod.n``.


Basic manipulation and comparison
--------------------------------------------------------------------------------


.. function:: void _nmod_vec_set(mp_ptr res, mp_srcptr vec, slong len)

    Copies ``len`` entries from the vector ``vec`` to ``res``.

.. function:: void _nmod_vec_zero(mp_ptr vec, slong len)

    Zeros the given vector of the given length.

.. function:: void _nmod_vec_swap(mp_ptr a, mp_ptr b, slong length)

    Swaps the vectors ``a`` and ``b`` of length `n` by actually
    swapping the entries.

.. function:: void _nmod_vec_reduce(mp_ptr res, mp_srcptr vec, slong len, nmod_t mod)

    Reduces the entries of ``(vec, len)`` modulo ``mod.n`` and set 
    ``res`` to the result.

.. function:: flint_bitcnt_t _nmod_vec_max_bits(mp_srcptr vec, slong len)

    Returns the maximum number of bits of any entry in the vector.

.. function:: int _nmod_vec_equal(mp_srcptr vec, mp_srcptr vec2, slong len)

    Returns~`1` if ``(vec, len)`` is equal to ``(vec2, len)``, 
    otherwise returns~`0`.


Arithmetic operations
--------------------------------------------------------------------------------


.. function:: void _nmod_vec_add(mp_ptr res, mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod)

    Sets ``(res, len)`` to the sum of ``(vec1, len)`` 
    and ``(vec2, len)``.

.. function:: void _nmod_vec_sub(mp_ptr res, mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod)

    Sets ``(res, len)`` to the difference of ``(vec1, len)`` 
    and ``(vec2, len)``.

.. function:: void _nmod_vec_neg(mp_ptr res, mp_srcptr vec, slong len, nmod_t mod)

    Sets ``(res, len)`` to the negation of ``(vec, len)``.

.. function:: void _nmod_vec_scalar_mul_nmod(mp_ptr res, mp_srcptr vec, slong len, mp_limb_t c, nmod_t mod)

    Sets ``(res, len)`` to ``(vec, len)`` multiplied by `c`. The element
    `c` and all elements of `vec` are assumed to be less than `mod.n`.

.. function:: void _nmod_vec_scalar_mul_nmod_shoup(mp_ptr res, mp_srcptr vec, slong len, mp_limb_t c, nmod_t mod)

    Sets ``(res, len)`` to ``(vec, len)`` multiplied by `c` using
    :func:`n_mulmod_shoup`. `mod.n` should be less than `2^{\mathtt{FLINT\_BITS} - 1}`. `c` 
    and all elements of `vec` should be less than `mod.n`.

.. function:: void _nmod_vec_scalar_addmul_nmod(mp_ptr res, mp_srcptr vec, slong len, mp_limb_t c, nmod_t mod)

    Adds ``(vec, len)`` times `c` to the vector ``(res, len)``. The element
    `c` and all elements of `vec` are assumed to be less than `mod.n`.


Dot products
--------------------------------------------------------------------------------


.. function:: int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod)

    Returns the number of limbs (0, 1, 2 or 3) needed to represent the
    unreduced dot product of two vectors of length ``len`` having entries
    modulo ``mod.n``, assuming that ``len`` is nonnegative and that
    ``mod.n`` is nonzero. The computed bound is tight. In other words,
    this function returns the precise limb size of ``len`` times
    ``(mod.n - 1) ^ 2``.

.. function:: macro NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, nlimbs)

    Effectively performs the computation::

        res = 0;
        for (i = 0; i < len; i++)
            res += (expr1) * (expr2);

    but with the arithmetic performed modulo ``mod``.
    The ``nlimbs`` parameter should be 0, 1, 2 or 3, specifying the
    number of limbs needed to represent the unreduced result.

.. function:: mp_limb_t _nmod_vec_dot(mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod, int nlimbs)

    Returns the dot product of (``vec1``, ``len``) and
    (``vec2``, ``len``). The ``nlimbs`` parameter should be
    0, 1, 2 or 3, specifying the number of limbs needed to represent the
    unreduced result.

.. function:: mp_limb_t _nmod_vec_dot_rev(mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod, int nlimbs)

    The same as ``_nmod_vec_dot``, but reverses ``vec2``.

.. function:: mp_limb_t _nmod_vec_dot_ptr(mp_srcptr vec1, const mp_ptr * vec2, slong offset, slong len, nmod_t mod, int nlimbs)

    Returns the dot product of (``vec1``, ``len``) and the values at
    ``vec2[i][offset]``. The ``nlimbs`` parameter should be
    0, 1, 2 or 3, specifying the number of limbs needed to represent the
    unreduced result.
