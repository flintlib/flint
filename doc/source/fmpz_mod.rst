.. _fmpz-mod:

**fmpz_mod.h** -- arithmetic modulo integers
===============================================================================

Description.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: fmpz_mod_ctx_struct

.. type:: fmpz_mod_ctx_t

    The context object for arithmetic modulo integers.


Context object
--------------------------------------------------------------------------------


.. function:: void fmpz_mod_ctx_init(fmpz_mod_ctx_t ctx, const fmpz_t n)

    Initialise ``ctx`` for arithmetic modulo ``n``, which is expected to be positive.

.. function:: void fmpz_mod_ctx_clear(fmpz_mod_ctx_t ctx)

    Free any memory used by ``ctx``.

.. function:: void fmpz_mod_ctx_set_modulus(fmpz_mod_ctx_t ctx, const fmpz_t n)

    Reconfigure ``ctx`` for arithmetic modulo ``n``.


Conversions
-----------------------------------------------------------------------------------------------------------------------

.. function:: void fmpz_mod_set_fmpz(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)

    Set ``a`` to ``b`` after reduction modulo the modulus.


Arithmetic
--------------------------------------------------------------------------------

Unless specified otherwise all functions here expect their relevant arguments to be in the canonical range `[0,n)`.
Comparison of elements against each other or against zero can be accomplished with func::fmpz_equal or func::fmpz_is_zero without a context.

.. function:: int fmpz_mod_is_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx)

    Return ``1`` if `a` is in the canonical range `[0,n)` and ``0`` otherwise.

.. function:: int fmpz_mod_is_one(const fmpz_t a, const fmpz_mod_ctx_t ctx)

    Return ``1`` if `a` is `1` modulo `n` and return ``0`` otherwise.

.. function:: void fmpz_mod_add(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

    Set `a` to `b+c` modulo `n`.

.. function:: void fmpz_mod_add_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
              void fmpz_mod_add_ui(fmpz_t a, const fmpz_t b, ulong c, const fmpz_mod_ctx_t ctx);
              void fmpz_mod_add_si(fmpz_t a, const fmpz_t b, slong c, const fmpz_mod_ctx_t ctx);

    Set `a` to `b+c` modulo `n` where only `b` is assumed to be canonical.

.. function:: void fmpz_mod_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

    Set `a` to `b-c` modulo `n`.

.. function:: void fmpz_mod_sub_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
              void fmpz_mod_sub_ui(fmpz_t a, const fmpz_t b, ulong c, const fmpz_mod_ctx_t ctx);
              void fmpz_mod_sub_si(fmpz_t a, const fmpz_t b, slong c, const fmpz_mod_ctx_t ctx);

    Set `a` to `b-c` modulo `n` where only `b` is assumed to be canonical.

.. function:: void fmpz_mod_fmpz_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
              void fmpz_mod_ui_sub(fmpz_t a, ulong b, const fmpz_t c, const fmpz_mod_ctx_t ctx);
              void fmpz_mod_si_sub(fmpz_t a, slong b, const fmpz_t c, const fmpz_mod_ctx_t ctx);

    Set `a` to `b-c` modulo `n` where only `c` is assumed to be canonical.

.. function:: void fmpz_mod_neg(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)

    Set `a` to `-b` modulo `n`.

.. function:: void fmpz_mod_mul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

    Set `a` to `b*c` modulo `n`.

.. function:: void fmpz_mod_inv(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)

    Set `a` to `b^{-1}` modulo `n`.
    This function expects that `b` is invertible modulo `n` and throws if this not the case.
    Invertibility maybe tested with func:`fmpz_mod_pow_fmpz` or func:`fmpz_mod_divides`.

.. function:: int fmpz_mod_divides(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)

    If `a*c = b \mod n` has a solution for `a` return `1` and set `a` to such a solution. Otherwise return `0` and leave `a` undefined.

.. function:: void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong e, const fmpz_mod_ctx_t ctx)

    Set `a` to `b^e` modulo `n`.

.. function:: int fmpz_mod_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e, const fmpz_mod_ctx_t ctx)

    Try to set `a` to `b^e` modulo `n`.
    If `e < 0` and `b` is not invertible modulo `n`, the return is `0`. Otherwise, the return is `1`.


Discrete Logarithms via Pohlig-Hellman
--------------------------------------------------------------------------------

.. function:: void fmpz_mod_discrete_log_pohlig_hellman_init(fmpz_mod_discrete_log_pohlig_hellman_t L)

    Initialize ``L``. Upon initialization ``L`` is not ready for computation.

.. function:: void fmpz_mod_discrete_log_pohlig_hellman_clear(fmpz_mod_discrete_log_pohlig_hellman_t L)

    Free any space used by ``L``.

.. function:: double fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(fmpz_mod_discrete_log_pohlig_hellman_t L, const fmpz_t p)

    Configure ``L`` for discrete logarithms modulo ``p`` to an internally chosen base. It is assumed that ``p`` is prime.
    The return is an estimate on the number of multiplications needed for one run.

.. function:: const fmpz * fmpz_mod_discrete_log_pohlig_hellman_primitive_root(const fmpz_mod_discrete_log_pohlig_hellman_t L)

    Return the internally stored base.

.. function:: void fmpz_mod_discrete_log_pohlig_hellman_run(fmpz_t x, const fmpz_mod_discrete_log_pohlig_hellman_t L, const fmpz_t y)

    Set ``x`` to the logarithm of ``y`` with respect to the internally stored base. ``y`` is expected to be reduced modulo the ``p``.
    The function is undefined if the logarithm does not exist.


.. function:: int fmpz_next_smooth_prime(fmpz_t a, const fmpz_t b)

    Either return `1` and set `a` to a smooth prime strictly greater than `b`, or return `0` and set `a` to `0`.
    The smooth primes returned by this function currently have no prime factor of `a-1` greater than `23`, but this should not be relied upon.

