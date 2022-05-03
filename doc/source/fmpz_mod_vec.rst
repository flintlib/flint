.. _fmpz-mod-vec:

**fmpz_mod_vec.h** -- vectors over integers mod n
===============================================================================

Description.

Conversions
--------------------------------------------------------------------------------

.. function:: void _fmpz_mod_vec_set_fmpz_vec(fmpz * A, const fmpz * B, slong len, const fmpz_mod_ctx_t ctx)

    Set the `fmpz_mod_vec` `(A, len)` to the `fmpz_vec` `(B, len)` after
    reduction of each entry modulo the modulus..

Arithmetic
--------------------------------------------------------------------------------

.. function:: void _fmpz_mod_vec_neg(fmpz * A, const fmpz * B, slong len, const fmpz_mod_ctx_t ctx)

    Set `(A, len)` to `-(B, len)`.

.. function:: void _fmpz_mod_vec_add(fmpz * a, const fmpz * b, const fmpz * c, slong n, const fmpz_mod_ctx_t ctx)

    Set (A, len)` to `(B, len) + (C, len)`.

.. function:: void _fmpz_mod_vec_sub(fmpz * a, const fmpz * b, const fmpz * c, slong n, const fmpz_mod_ctx_t ctx)

    Set (A, len)` to `(B, len) - (C, len)`.


Scalar Multiplication
--------------------------------------------------------------------------------

.. function:: void _fmpz_mod_vec_scalar_mul_fmpz_mod(fmpz * A, const fmpz * B, slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx)

    Set `(A, len)` to `(B, len)*c`.

.. function:: void _fmpz_mod_vec_scalar_addmul_fmpz_mod(fmpz * A, const fmpz * B, slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx)

    Set `(A, len)` to `(A, len) + (B, len)*c`.

.. function:: void _fmpz_mod_vec_scalar_mul_fmpz_mod(fmpz * A, const fmpz * B, slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx)

    Set `(A, len)` to `(B, len)/c` assuming `c` is nonzero.


Dot Product
--------------------------------------------------------------------------------

.. function:: void _fmpz_vec_dot(fmpz_t d, const fmpz * A, const fmpz * B, slong len, const fmpz_mod_ctx_t ctx)

    Set `d` to the dot product of `(A, len)` with `(B, len)`.

.. function:: void _fmpz_vec_dot_rev(fmpz_t d, const fmpz * A, const fmpz * B, slong len, const fmpz_mod_ctx_t ctx)

    Set `d` to the dot product of `(A, len)` with the reverse of the vector `(B, len)`.


Multiplication
--------------------------------------------------------------------------------

.. function:: void _fmpz_mod_vec_mul(fmpz * A, const fmpz * B, const fmpz * C, slong len, const fmpz_mod_ctx_t ctx)

    Set `(A, len)` the pointwise multiplication of `(B, len)` and `(C, len)`.
