.. _nmod:

**nmod.h** -- integers mod n (word-size n)
===============================================================================

Generic rings
--------------------------------------------------------------------------------

.. function:: int gr_ctx_init_nmod(gr_ctx_t ctx, ulong n)
              int gr_ctx_init_nmod8(gr_ctx_t ctx, ulong n)
              int gr_ctx_init_nmod32(gr_ctx_t ctx, ulong n)
              int gr_ctx_init_nmod_redc(gr_ctx_t ctx, ulong n)
              int gr_ctx_init_nmod_redc_fast(gr_ctx_t ctx, ulong n)

    Initialize the context object ``ctx`` for generic arithmetic in the ring
    `R = \mathbb{Z} / n \mathbb{Z}` for word-size `n`
    using the respective representation.

    The element datatype is :type:`ulong` for all representations
    except for ``nmod8`` and ``nmod32`` where the type is
    :type:`uint8` or :type:`uint32` respectively.

    If the modulus is valid, returns ``GR_SUCCESS``. Otherwise, the context
    object is not initialized and an error code is returned:

    * ``GR_DOMAIN`` is returned if `n = 0`.
    * ``GR_UNABLE`` is returned if `n` does not satisfy the technical
      restrictions of a specific representation. For example, ``nmod8`` requires
      `n \le 255` and ``nmod_redc`` requires that `n` is odd.

.. note ::

    Some generic algorithms need to check whether
    `R` is a field, i.e. whether `n` is a prime number.
    To avoid expensive on-the-fly primality tests, it is recommended
    to call ``gr_ctx_set_is_field(ctx, T_TRUE)`` (if `n` is prime)
    ``gr_ctx_set_is_field(ctx, T_FALSE)`` (if `n` is composite)
    after constructing a context object.
    This is not done automatically as it would slow down creating a context
    object in the case where one is just interested in basic arithmetic.

.. note ::

    Presently, many operations for ``nmod8``, ``nmod32``, ``nmod_redc`` and
    ``nmod_redc_fast`` are not as optimized as those for the general-purpose
    ``nmod``. It is currently recommended to use ``nmod8`` and ``nmod32``
    only if one specifically wants to minimize memory usage.

Modular reduction and arithmetic
--------------------------------------------------------------------------------

.. function:: void nmod_init(nmod_t * mod, ulong n)

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

.. function:: ulong _nmod_add(ulong a, ulong b, nmod_t mod)

    Returns `a + b` modulo ``mod.n``. It is assumed that ``mod`` is 
    no more than ``FLINT_BITS - 1`` bits. It is assumed that `a` and `b` 
    are already reduced modulo ``mod.n``.

.. function:: ulong nmod_add(ulong a, ulong b, nmod_t mod)

    Returns `a + b` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` and `b` are already reduced 
    modulo ``mod.n``.

.. function:: ulong nmod_ui_add_ui(ulong a, ulong b, nmod_t mod)

    Returns `a + b` modulo ``mod.n``. Does not require that `a` and `b` are
    already reduced modulo ``mod.n``.

.. function:: ulong _nmod_sub(ulong a, ulong b, nmod_t mod)

    Returns `a - b` modulo ``mod.n``. It is assumed that ``mod`` 
    is no more than ``FLINT_BITS - 1`` bits. It is assumed that 
    `a` and `b` are already reduced modulo ``mod.n``.

.. function:: ulong nmod_sub(ulong a, ulong b, nmod_t mod)

    Returns `a - b` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` and `b` are already reduced 
    modulo ``mod.n``.

.. function:: ulong nmod_neg(ulong a, nmod_t mod)

    Returns `-a` modulo ``mod.n``. It is assumed that `a` is already 
    reduced modulo ``mod.n``, but no assumptions are made about the 
    latter.

.. function:: ulong nmod_mul(ulong a, ulong b, nmod_t mod)

    Returns `ab` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` and `b` are already reduced 
    modulo ``mod.n``.

.. function:: ulong _nmod_mul_fullword(ulong a, ulong b, nmod_t mod)

    Returns `ab` modulo ``mod.n``. Requires that ``mod.n`` is exactly
    ``FLINT_BITS`` large. It is assumed that `a` and `b` are already
    reduced modulo ``mod.n``.

.. function:: ulong nmod_ui_mul_ui(ulong a, ulong b, nmod_t mod)

    Returns `ab` modulo ``mod.n``. Does not require that `a` and `b` are
    already reduced modulo ``mod.n``.

.. function:: ulong nmod_inv(ulong a, nmod_t mod)

    Returns `a^{-1}` modulo ``mod.n``. The inverse is assumed to exist.

.. function:: ulong nmod_div(ulong a, ulong b, nmod_t mod)

    Returns `ab^{-1}` modulo ``mod.n``. The inverse of `b` is assumed to
    exist. It is assumed that `a` is already reduced modulo ``mod.n``.

.. function:: int nmod_divides(ulong * a, ulong b, ulong c, nmod_t mod)

    If `a\cdot c = b \mod n` has a solution for `a` return `1` and set `a` to such a solution. Otherwise return `0` and leave `a` undefined.

.. function:: ulong nmod_pow_ui(ulong a, ulong e, nmod_t mod)

    Returns `a^e` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` is already reduced
    modulo ``mod.n``.

.. function:: ulong nmod_pow_fmpz(ulong a, const fmpz_t e, nmod_t mod)

    Returns `a^e` modulo ``mod.n``. No assumptions are made about 
    ``mod.n``. It is assumed that `a` is already reduced
    modulo ``mod.n`` and that `e` is not negative.

.. function:: ulong nmod_ui_pow_ui(ulong a, ulong e, nmod_t mod)

    Returns `a^e` modulo ``mod.n``. Does not require that `a` is already reduced
    modulo ``mod.n``.

.. function:: ulong nmod_2_pow_ui(ulong e, nmod_t mod)

    Returns `2^e` modulo ``mod.n``.

Montgomery arithmetic
--------------------------------------------------------------------------------

Let `n` be an odd integer smaller than the machine word modulus
`R = 2^{32}` or `R = 2^{64}`.
The Montgomery representation of an integer `x` is the residue `x R \bmod n`.
We can use Montgomery representations
for doing arithmetic in `\mathbb{Z} / n \mathbb{Z}`,
following the rules

.. math ::

    (x + y) R \bmod n = x R + y R \bmod n,

.. math ::

    (x y) R \bmod n = (x R) (y R) / R \bmod n.

The advantage of using the Montgomery representation instead of
the standard representation `x \bmod n` is that it allows
for faster multiplication. Montgomery arithmetic is also known as
REDC arithmetic.

By default moduli up to the full word size are allowed and residues are
strictly canonicalised to `[0, n)`. We provide sets of alternative
methods for restricted moduli which, depending on the machine,
can speed up arithmetic:

* The ``half`` (and ``half_fast``) methods work with half-length
  `R`, avoiding the need for double-word operations in the internal
  arithmetic. They are not necessarily faster than the non-``half`` versions.
  Mixing ``half`` and non-``half`` methods and context objects
  is not allowed.

* The ``fast`` (and ``half_fast``) accept non-canonical residues in `[0, 2n)`
  and return non-canonical residues in `[0, 2n)`. The user must convert back
  to the canonical range `[0, n)` before using any non-``fast`` functions.
  The ``fast`` methods are generally faster than the non-``fast`` methods.

The restrictions are summarized in the following table, assuming a
64-bit word size:

+--------------------+-----------------+----------+----------------------+
| Prefix             |   Maximum `n`   |   `R`    |   Canonicalisation   |
+====================+=================+==========+======================+
| ``redc``           |  `2^{64} - 1`   | `2^{64}` |  `[0, n)`            |
+--------------------+-----------------+----------+----------------------+
| ``redc_fast``      |  `2^{62} - 1`   | `2^{64}` |  `[0, 2n)`           |
+--------------------+-----------------+----------+----------------------+
| ``redc_half``      |  `2^{31} - 1`   | `2^{32}` |  `[0, n)`            |
+--------------------+-----------------+----------+----------------------+
| ``redc_half_fast`` |  `2^{30} - 1`   | `2^{32}` |  `[0, 2n)`           |
+--------------------+-----------------+----------+----------------------+

For 32-bit machines, the maximum `n` become `2^{32} - 1`,
`2^{30} - 1`, `2^{15} - 1`, and `2^{14} - 1` respectively.

.. type:: nmod_redc_ctx_struct
          nmod_redc_ctx_t

    A context object for Montgomery arithmetic. This holds the same content
    as an ``nmod_t`` context, plus the precomputed constant `-1/n` modulo `R`.

.. function:: void nmod_redc_ctx_init_nmod(nmod_redc_ctx_t ctx, nmod_t mod)
              void nmod_redc_ctx_init_ui(nmod_redc_ctx_t ctx, ulong n)

    Initialize ``ctx`` for Montgomery arithmetic modulo the given modulus
    `n` which is required to be odd.

.. function:: ulong nmod_redc_set_nmod(ulong x, const nmod_redc_ctx_t ctx)

    Convert `x` from the standard representation to Montgomery representation.
    If `x` is viewed a residue in the standard representation,
    this returns `x R \bmod n`.

.. function:: ulong nmod_redc_set_ui(ulong x, const nmod_redc_ctx_t ctx)

    Convert `x` which is not necessarily reduced modulo `n`
    to Montgomery representation.

.. function:: ulong nmod_redc_get_nmod(ulong x, const nmod_redc_ctx_t ctx)

    Convert `x` back from Montgomery representation to standard representation.
    If `x` is viewed as a residue in the standard representation,
    this returns `x / R \bmod n`. This function allows `x \in [0, 2n)`.

.. function:: ulong nmod_redc_neg(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)

    Arithmetic in the Montgomery representation.

.. function:: int nmod_redc_can_use_fast(const nmod_redc_ctx_t ctx)

    Return whether the modulus is small enough to safely use fast operations.

.. function:: ulong nmod_redc_fast_normalise(ulong x, const nmod_redc_ctx_t ctx)

    Convert a non-canonical residue in `[0, 2n)` to a canonical
    residue in `[0, n)`.

.. function:: ulong nmod_redc_fast_neg(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_fast_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_fast_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_fast_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_fast_mul_two(ulong x, const nmod_redc_ctx_t ctx)

    Arithmetic in Montgomery representation, using non-canonical residues.

.. function:: ulong _nmod_redc_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx)
              ulong _nmod_redc_fast_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx)

    Return `a^{exp}`. The exponent is required to be positive.

.. function:: ulong _nmod_redc_2_pow_ui(ulong exp, const nmod_redc_ctx_t ctx)
              ulong _nmod_redc_fast_2_pow_ui(ulong exp, const nmod_redc_ctx_t ctx)

    Return `2^{exp}`. There are no restrictions on the exponent.

.. function:: void nmod_redc_half_ctx_init_nmod(nmod_redc_ctx_t ctx, nmod_t mod)
              void nmod_redc_half_ctx_init_ui(nmod_redc_ctx_t ctx, ulong n)

    Initialize context for use with ``half`` or ``half_fast`` methods.

.. function:: ulong nmod_redc_half_set_nmod(ulong x, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_half_set_ui(ulong x, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_half_get_nmod(ulong x, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_half_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_half_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_half_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              int nmod_redc_half_can_use_fast(const nmod_redc_ctx_t ctx)
              ulong nmod_redc_half_fast_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_half_fast_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
              ulong nmod_redc_half_fast_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)

    Methods analogous to their non-``half`` counterparts.

Discrete Logarithms via Pohlig-Hellman
--------------------------------------------------------------------------------

.. function:: void nmod_discrete_log_pohlig_hellman_init(nmod_discrete_log_pohlig_hellman_t L)

    Initialize ``L``. Upon initialization ``L`` is not ready for computation.

.. function:: void nmod_discrete_log_pohlig_hellman_clear(nmod_discrete_log_pohlig_hellman_t L)

    Free any space used by ``L``.

.. function:: double nmod_discrete_log_pohlig_hellman_precompute_prime(nmod_discrete_log_pohlig_hellman_t L, ulong p)

    Configure ``L`` for discrete logarithms modulo ``p`` to an internally chosen base. It is assumed that ``p`` is prime.
    The return is an estimate on the number of multiplications needed for one run.

.. function:: ulong nmod_discrete_log_pohlig_hellman_primitive_root(const nmod_discrete_log_pohlig_hellman_t L)

    Return the internally stored base.

.. function:: ulong nmod_discrete_log_pohlig_hellman_run(const nmod_discrete_log_pohlig_hellman_t L, ulong y)

    Return the logarithm of ``y`` with respect to the internally stored base. ``y`` is expected to be reduced modulo the ``p``.
    The function is undefined if the logarithm does not exist.
