.. _fft-small:

**fft_small.h** -- FFT modulo word-size primes
===============================================================================

This module currently requires building FLINT with support for
AVX2 or NEON instructions.

Integer multiplication
--------------------------------------------------------------------------------

.. type:: mpn_ctx_struct
          mpn_ctx_t

    Context object for multiplications allowing non-FFT moduli.
    The structure contains FFT context objects for multiple FFT primes
    (currently 8) together with tables for Chinese remaindering.

.. function:: void mpn_ctx_init(mpn_ctx_t R, ulong p)

    Initialize multiplication context object with initial prime ``p``.

.. function:: void mpn_ctx_clear(mpn_ctx_t R)

    Free memory allocated by the context object.

.. function:: mpn_ctx_struct * get_default_mpn_ctx(void)

    Return a pointer to a cached thread-local context object used by default
    for multiplications. Calling :func:`flint_cleanup` or :func:`flint_cleanup_master`
    frees the cache.

.. function:: void mpn_ctx_mpn_mul(mpn_ctx_t R, ulong * r1, const ulong * i1, ulong n1, const ulong * i2, ulong n2)
              void mpn_mul_default_mpn_ctx(mp_ptr r1, mp_srcptr i1, mp_size_t n1, mp_srcptr i2, mp_size_t n2)

    Writes to ``r1`` the product of the integers ``(i1, n1)`` and ``(i2, n2)``.
    Assumes that `n_1 \ge n_2 \ge 1`, respectively using a given context
    object ``R`` or the default thread-local object.

Polynomial arithmetic
---------------------------------------------------------------------------------

.. function:: void _nmod_poly_mul_mid_mpn_ctx(ulong * z, ulong zl, ulong zh, const ulong * a, ulong an, const ulong * b, ulong bn, nmod_t mod, mpn_ctx_t R)
              void _nmod_poly_mul_mid_default_mpn_ctx(mp_ptr res, slong zl, slong zh, mp_srcptr a, slong an, mp_srcptr b, slong bn, nmod_t mod)

    Writes to ``z`` the middle product containing coefficients in the
    range `[zl, zh)` of the product of the polynomials  ``(a, an)`` and ``(b, bn)``,
    respectively using a given context object ``R`` or the default thread-local object.
    Assumes that `an \ge bn \ge 1`.

.. function:: int _fmpz_poly_mul_mid_mpn_ctx(fmpz * z, ulong zl, ulong zh, const fmpz * a, ulong an, const fmpz * b, ulong bn, mpn_ctx_t R)
              int _fmpz_poly_mul_mid_default_mpn_ctx(fmpz * z, ulong zl, ulong zh, const fmpz * a, ulong an, const fmpz * b, ulong bn)

    Like the ``nmod`` functions. Performs the multiplication and returns 1
    if there are sufficiently many primes ``R`` to compute the result;
    otherwise returns 0 without touching the output.

.. function:: void _nmod_poly_divrem_mpn_ctx(ulong * q, ulong * r, const ulong * a, ulong an, const ulong * b, ulong bn, nmod_t mod, mpn_ctx_t R)

    Polynomial division with remainder.

Preconditioned polynomial arithmetic
---------------------------------------------------------------------------------

.. type:: mul_precomp_struct

.. function:: void _mul_precomp_init(mul_precomp_struct * M, const ulong * b, ulong bn, ulong btrunc, ulong depth, nmod_t mod, mpn_ctx_t R)
              void _mul_precomp_clear(mul_precomp_struct * M)

    Represents ``(b, bn)`` in transformed form for preconditioned multiplication.

.. function:: int _nmod_poly_mul_mid_precomp(ulong * z, ulong zl, ulong zh, const ulong * a, ulong an, mul_precomp_struct * M, nmod_t mod, mpn_ctx_t R)

    Polynomial multiplication given a precomputed transform ``M``.
    Returns 1 if successful, 0 if the precomputed transform is too short.

.. type:: nmod_poly_divrem_precomp_struct

.. function:: void _nmod_poly_divrem_precomp_init(nmod_poly_divrem_precomp_struct * M, const ulong * b, ulong bn, ulong Bn, nmod_t mod, mpn_ctx_t R)
              void _nmod_poly_divrem_precomp_clear(nmod_poly_divrem_precomp_struct * M)

    Represents ``(b, bn)`` and its inverse in transformed form for preconditioned multiplication.

.. function:: int _nmod_poly_divrem_precomp(ulong * q, ulong * r, const ulong * a, ulong an, nmod_poly_divrem_precomp_struct * M, nmod_t mod, mpn_ctx_t R)

    Polynomial multiplication given a precomputed transform ``M``.
    Returns 1 if successful, 0 if the precomputed transform is too short.
