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


Modular reduction and arithmetic
--------------------------------------------------------------------------------


.. function:: void nmod_init(nmod_t * mod, mp_limb_t n)

    Initialises the given ``nmod_t`` structure for reduction modulo `n`
    with a precomputed inverse.

.. macro:: NMOD_RED2(r, a_hi, a_lo, mod)

    Macro to set `r` to `a` reduced modulo ``mod.n``, where `a` 
    consists of two limbs ``(a_hi, a_lo)``. The ``mod`` parameter 
    must be a valid ``nmod_t`` structure. It is assumed that ``a_hi`` 
    is already reduced modulo ``mod.n``.

.. macro:: NMOD_RED(r, a, mod)

    Macro to set `r` to `a` reduced modulo ``mod.n``. The ``mod`` 
    parameter must be a valid ``nmod_t`` structure.

.. macro:: NMOD2_RED2(r, a_hi, a_lo, mod)

    Macro to set `r` to `a` reduced modulo ``mod.n``, where `a` 
    consists of two limbs ``(a_hi, a_lo)``. The ``mod`` parameter 
    must be a valid ``nmod_t`` structure. No assumptions are made 
    about ``a_hi``.

.. macro:: NMOD_RED3(r, a_hi, a_me, a_lo, mod)

    Macro to set `r` to `a` reduced modulo ``mod.n``, where `a` 
    consists of three limbs ``(a_hi, a_me, a_lo)``. The ``mod`` 
    parameter must be a valid ``nmod_t`` structure. It is assumed 
    that ``a_hi`` is already reduced modulo ``mod.n``.

.. macro:: NMOD_ADDMUL(r, a, b, mod)

    Macro to set `r` to `r + ab` reduced modulo ``mod.n``. The 
    ``mod`` parameter must be a valid ``nmod_t`` structure. It is 
    assumed that `r`, `a`, `b` are already reduced modulo ``mod.n``.

.. function:: mp_limb_t _nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)

    Returns `a + b` modulo ``mod.n``. It is assumed that ``mod`` is 
    no more than ``FLINT_BITS - 1`` bits. It is assumed that `a` and `b` 
    are already reduced modulo ``mod.n``.

.. function:: mp_limb_t nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)

    Returns `a + b` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` and `b` are already reduced 
    modulo ``mod.n``.

.. function:: mp_limb_t _nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)

    Returns `a - b` modulo ``mod.n``. It is assumed that ``mod`` 
    is no more than ``FLINT_BITS - 1`` bits. It is assumed that 
    `a` and `b` are already reduced modulo ``mod.n``.

.. function:: mp_limb_t nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)

    Returns `a - b` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` and `b` are already reduced 
    modulo ``mod.n``.

.. function:: mp_limb_t nmod_neg(mp_limb_t a, nmod_t mod)

    Returns `-a` modulo ``mod.n``. It is assumed that `a` is already 
    reduced modulo ``mod.n``, but no assumptions are made about the 
    latter.

.. function:: mp_limb_t nmod_mul(mp_limb_t a, mp_limb_t b, nmod_t mod)

    Returns `ab` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` and `b` are already reduced 
    modulo ``mod.n``.

.. function:: mp_limb_t nmod_inv(mp_limb_t a, nmod_t mod)

    Returns `a^{-1}` modulo ``mod.n``. The inverse is assumed to exist.

.. function:: mp_limb_t nmod_div(mp_limb_t a, mp_limb_t b, nmod_t mod)

    Returns `ab^{-1}` modulo ``mod.n``. The inverse of `b` is assumed to
    exist. It is assumed that `a` is already reduced modulo ``mod.n``.

.. function:: mp_limb_t nmod_pow_ui(mp_limb_t a, ulong e, nmod_t mod)

    Returns `a^e` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` is already reduced
    modulo ``mod.n``.

.. function:: mp_limb_t nmod_pow_fmpz(mp_limb_t a, const fmpz_t e, nmod_t mod)

    Returns `a^e` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` is already reduced
    modulo ``mod.n`` and that `e` is not negative.



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


Discrete Logarithms via Pohlig-Hellman
--------------------------------------------------------------------------------

.. function:: void nmod_discrete_log_pohlig_hellman_init(nmod_discrete_log_pohlig_hellman_t L)

    Initialize ``L``. Upon initialization ``L`` is not ready for computation.

.. function:: void nmod_discrete_log_pohlig_hellman_clear(nmod_discrete_log_pohlig_hellman_t L)

    Free any space used by ``L``.

.. function:: double nmod_discrete_log_pohlig_hellman_precompute_prime(nmod_discrete_log_pohlig_hellman_t L, mp_limb_t p)

    Configure ``L`` for discrete logarithms modulo ``p`` to an internally chosen base. It is assumed that ``p`` is prime.
    The return is an estimate on the number of multiplications needed for one run.

.. function:: mp_limb_t nmod_discrete_log_pohlig_hellman_primitive_root(const nmod_discrete_log_pohlig_hellman_t L)

    Return the internally stored base.

.. function:: ulong nmod_discrete_log_pohlig_hellman_run(const nmod_discrete_log_pohlig_hellman_t L, mp_limb_t y)

    Return the logarithm of ``y`` with respect to the internally stored base. ``y`` is expected to be reduced modulo the ``p``.
    The function is undefined if the logarithm does not exist.
