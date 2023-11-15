.. _nmod:

**nmod.h** -- integers mod n (word-size n)
===============================================================================

Modular reduction and arithmetic
--------------------------------------------------------------------------------

.. function:: void nmod_init(nmod_t * mod, mp_limb_t n)

    Initialises the given ``nmod_t`` structure for reduction modulo `n`
    with a precomputed inverse.

.. macro:: NMOD_BITS(mod)

    Macro giving the number of bits in ``mod.n``.

.. macro:: NMOD_CAN_USE_SHOUP(mod)

    Macro returning whether Shoup's algorithm can be used for
    preconditioned multiplication mod ``mod.n``.

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

.. macro:: NMOD_MUL_PRENORM(res, a, b, mod)

    Macro to set `r` to `ab` modulo ``mod.n``. The 
    ``mod`` parameter must be a valid ``nmod_t`` structure. It is 
    assumed that `a`, `b` are already reduced modulo ``mod.n``
    and that either `a` or `b` is prenormalised by left-shifting
    by ``mod.norm``.

.. macro:: NMOD_MUL_FULLWORD(res, a, b, mod)

    Macro to set `r` to `ab` modulo ``mod.n``. The 
    ``mod`` parameter must be a valid ``nmod_t`` structure. It is 
    assumed that `a`, `b` are already reduced modulo ``mod.n``
    and that ``mod.n`` is exactly ``FLINT_BITS`` bits large.

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

.. function:: mp_limb_t _nmod_mul_fullword(mp_limb_t a, mp_limb_t b, nmod_t mod)

    Returns `ab` modulo ``mod.n``. Requires that ``mod.n`` is exactly
    ``FLINT_BITS`` large. It is assumed that `a` and `b` are already
    reduced modulo ``mod.n``.

.. function:: mp_limb_t nmod_inv(mp_limb_t a, nmod_t mod)

    Returns `a^{-1}` modulo ``mod.n``. The inverse is assumed to exist.

.. function:: mp_limb_t nmod_div(mp_limb_t a, mp_limb_t b, nmod_t mod)

    Returns `ab^{-1}` modulo ``mod.n``. The inverse of `b` is assumed to
    exist. It is assumed that `a` is already reduced modulo ``mod.n``.

.. function:: int nmod_divides(mp_limb_t * a, mp_limb_t b, mp_limb_t c, nmod_t mod)

    If `a\cdot c = b \mod n` has a solution for `a` return `1` and set `a` to such a solution. Otherwise return `0` and leave `a` undefined.

.. function:: mp_limb_t nmod_pow_ui(mp_limb_t a, ulong e, nmod_t mod)

    Returns `a^e` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` is already reduced
    modulo ``mod.n``.

.. function:: mp_limb_t nmod_pow_fmpz(mp_limb_t a, const fmpz_t e, nmod_t mod)

    Returns `a^e` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` is already reduced
    modulo ``mod.n`` and that `e` is not negative.


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
