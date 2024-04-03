.. _mpn-mod:

**mpn_mod.h** -- integers mod n (packed multi-word n)
===============================================================================

This module provides efficient arithmetic in rings
`R = \mathbb{Z} / n \mathbb{Z}` for medium-sized `n`.
Given an `\ell`-limb modulus `2^{\beta (\ell-1)} \le n < 2^{\beta \ell}`
where `\beta` is ``FLINT_BITS`` (32 or 64),
elements are represented as `\ell`-limb arrays (i.e. ``mp_limb_t[l]``),
zero-padded for values that happen to fit in less than `\ell` limbs,
which can be stack-allocated and packed consecutively
without indirection or memory allocation overhead.

This module is designed to use the :ref:`generics <gr>` interface.
As such, the ring is represented by a :type:`gr_ctx_t` context object,
methods return status flags, and one can use
generic structures such as :type:`gr_poly_t` for
polynomials and :type:`gr_mat_t` for matrices.

Types, macros and constants
-------------------------------------------------------------------------------

.. macro :: MPN_MOD_MIN_LIMBS
            MPN_MOD_MAX_LIMBS

    The number of limbs `\ell` permitted in a modulus. The current limits
    are `2 \le \ell \le 16`.
    We exclude single-limb moduli since these are covered by
    :ref:`nmod <nmod>` arithmetic, and this allows not bothering
    with various degenerate cases.
    The upper limit exists so that elements and temporary buffers
    are safe to allocate on the stack and so that simple operations
    like swapping or zeroing elements are not too expensive
    compared to a pointer-and-size representation.
    A second reason is that the algorithms in this module have
    been tuned only for moduli in a certain range.
    For larger moduli, one should use :ref:`fmpz_mod <fmpz-mod>` instead.
    The upper limit might be increased in the future.

.. function:: int gr_ctx_init_mpn_mod(gr_ctx_t ctx, const fmpz_t n)

    Initializes *ctx* to the ring `\mathbb{Z}/n\mathbb{Z}`
    of integers modulo *n* where elements are ``mp_limb_t`` arrays with
    the same number of limbs as *n*. This constructor does no
    initialization and returns
    ``GR_UNABLE`` if the modulus is not in bounds.

.. macro:: MPN_MOD_CTX_NLIMBS(ctx)

    Retrives the number of limbs `\ell` of the modulus.

.. macro:: MPN_MOD_CTX_MODULUS(ctx)

    Pointer to the limbs of the modulus.

.. macro:: MPN_MOD_CTX_NORM(ctx)

    An integer indicating the number of leading zero bits in the most
    significant limb of the modulus.

.. macro:: MPN_MOD_CTX_MODULUS_NORMED(ctx)

    Pointer to a copy of the modulus left-shifted so that the
    most significant bit is in a limb boundary.

.. macro:: MPN_MOD_CTX_MODULUS_PREINV(ctx)

    Pointer to a precomputed inverse of the (normed) modulus.

.. macro:: MPN_MOD_CTX_IS_PRIME(ctx)

    A :type:`truth_t` flag indicating whether `n` is prime.

