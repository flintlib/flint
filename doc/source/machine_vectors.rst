.. _machine-vectors:

**machine_vectors.h** -- SIMD-accelerated operations on fixed-length vectors
===============================================================================

This module currently requires building FLINT with support for
AVX2 or NEON instructions.

Some functions may require that vectors are aligned in memory.

Types
-------------------------------------------------------------------------------

.. type:: vec1n
          vec2n
          vec4n
          vec8n

    Vector with 1, 2, 4, or 8 :type:`ulong` entries.

.. type:: vec1d
          vec2d
          vec4d
          vec8d

    Vector with 1, 2, 4, or 8 ``double`` entries.

Printing
-------------------------------------------------------------------------------

.. function:: void vec4d_print(vec4d a)
              void vec4n_print(vec4n a)

Access and conversions
-------------------------------------------------------------------------------

.. function:: vec1d vec1d_load(const double * a)
              vec4d vec4d_load(const double * a)
              vec8d vec8d_load(const double * a)

.. function:: vec1d vec1d_load_aligned(const double * a)
              vec4d vec4d_load_aligned(const double * a)
              vec8d vec8d_load_aligned(const double * a)

.. function:: vec1d vec1d_load_unaligned(const double * a)
              vec4d vec4d_load_unaligned(const double * a)
              vec8d vec8d_load_unaligned(const double * a)
              vec4n vec4n_load_unaligned(const ulong * a)
              vec8n vec8n_load_unaligned(const ulong * a)

.. function:: void vec1d_store(double * z, vec1d a)
              void vec4d_store(double * z, vec4d a)
              void vec8d_store(double * z, vec8d a)

.. function:: void vec1d_store_aligned(double * z, vec1d a)
              void vec4d_store_aligned(double * z, vec4d a)
              void vec8d_store_aligned(double * z, vec8d a)

.. function:: void vec1d_store_unaligned(double * z, vec1d a)
              void vec4d_store_unaligned(double * z, vec4d a)
              void vec4n_store_unaligned(ulong * z, vec4n a)
              void vec8d_store_unaligned(double * z, vec8d a)

.. function:: double vec4d_get_index(vec4d a, const int i)
              double vec8d_get_index(vec8d a, int i)

    Extract the entry at index `i`.

.. function:: vec1d vec1d_set_d(double a)
              vec4d vec4d_set_d(double a)
              vec4n vec4n_set_n(ulong a)
              vec8d vec8d_set_d(double a)
              vec8n vec8n_set_n(ulong a)

    Set all entries to the same value.

.. function:: vec4d vec4d_set_d4(double a0, double a1, double a2, double a3)
              vec4n vec4n_set_n4(ulong a0, ulong a1, ulong a2, ulong a3)
              vec8d vec8d_set_d8(double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7)

    Create vector from distinct entries.

.. function:: vec4n vec4d_convert_limited_vec4n(vec4d a)
              vec8d vec8n_convert_limited_vec8d(vec8n a)

Permutations
-------------------------------------------------------------------------------

.. function:: vec4d vec4d_unpacklo(vec4d a, vec4d b)
              vec4d vec4d_unpackhi(vec4d a, vec4d b)
              vec4d vec4d_permute_0_2_1_3(vec4d a)
              vec4d vec4d_permute_3_1_2_0(vec4d a)
              vec4d vec4d_permute_3_2_1_0(vec4d a)
              vec4d vec4d_permute2_0_2(vec4d a, vec4d b)
              vec4d vec4d_permute2_1_3(vec4d a, vec4d b)
              vec4d vec4d_unpack_lo_permute_0_2_1_3(vec4d u, vec4d v)
              vec4d vec4d_unpack_hi_permute_0_2_1_3(vec4d u, vec4d v)
              vec4d vec4d_unpackhi_permute_3_1_2_0(vec4d u, vec4d v)
              vec4d vec4d_unpacklo_permute_3_1_2_0(vec4d u, vec4d v)

.. macro:: VEC4D_TRANSPOSE(z0, z1, z2, z3, a0, a1, a2, a3)

    Sets the rows ``z`` to the transpose of the 4x4 matrix
    given by rows ``a``.

Comparisons
-------------------------------------------------------------------------------

.. function:: int vec1d_same(double a, double b)
              int vec4d_same(vec4d a, vec4d b)
              int vec8d_same(vec8d a, vec8d b)

    Check whether the vectors are equal.

.. function:: vec4d vec4d_cmp_ge(vec4d a, vec4d b)
              vec4d vec4d_cmp_gt(vec4d a, vec4d b)

    Entrywise comparisons.

Arithmetic and basic operations
-------------------------------------------------------------------------------

.. function:: vec1d vec1d_round(vec1d a)
              vec4d vec4d_round(vec4d a)
              vec8d vec8d_round(vec8d a)

.. function:: vec1d vec1d_zero()
              vec4d vec4d_zero()
              vec8d vec8d_zero()

.. function:: vec1d vec1d_one()
              vec4d vec4d_one()
              vec8d vec8d_one()

.. function:: vec1d vec1d_add(vec1d a, vec1d b)
              vec1d vec1d_sub(vec1d a, vec1d b)
              vec4d vec4d_add(vec4d a, vec4d b)
              vec4d vec4d_sub(vec4d a, vec4d b)
              vec4n vec4n_add(vec4n a, vec4n b)
              vec4n vec4n_sub(vec4n a, vec4n b)
              vec8d vec8d_add(vec8d a, vec8d b)
              vec8d vec8d_sub(vec8d a, vec8d b)

.. function:: vec1d vec1d_addsub(vec1d a, vec1d b)
              vec4d vec4d_addsub(vec4d a, vec4d b)

.. function:: vec1d vec1d_neg(vec1d a)
              vec4d vec4d_neg(vec4d a)
              vec8d vec8d_neg(vec8d a)

.. function:: vec1d vec1d_abs(vec1d a)
              vec4d vec4d_abs(vec4d a)

.. function:: vec1d vec1d_max(vec1d a, vec1d b)
              vec1d vec1d_min(vec1d a, vec1d b)
              vec4d vec4d_max(vec4d a, vec4d b)
              vec4d vec4d_min(vec4d a, vec4d b)
              vec8d vec8d_max(vec8d a, vec8d b)
              vec8d vec8d_min(vec8d a, vec8d b)

.. function:: vec1d vec1d_mul(vec1d a, vec1d b)
              vec4d vec4d_mul(vec4d a, vec4d b)
              vec8d vec8d_mul(vec8d a, vec8d b)

.. function:: vec1d vec1d_half(vec1d a)
              vec4d vec4d_half(vec4d a)

.. function:: vec1d vec1d_div(vec1d a, vec1d b)
              vec4d vec4d_div(vec4d a, vec4d b)
              vec8d vec8d_div(vec8d a, vec8d b)

.. function:: vec1d vec1d_fmadd(vec1d a, vec1d b, vec1d c)
              vec4d vec4d_fmadd(vec4d a, vec4d b, vec4d c)
              vec8d vec8d_fmadd(vec8d a, vec8d b, vec8d c)

.. function:: vec1d vec1d_fmsub(vec1d a, vec1d b, vec1d c)
              vec4d vec4d_fmsub(vec4d a, vec4d b, vec4d c)
              vec8d vec8d_fmsub(vec8d a, vec8d b, vec8d c)

.. function:: vec1d vec1d_fnmadd(vec1d a, vec1d b, vec1d c)
              vec4d vec4d_fnmadd(vec4d a, vec4d b, vec4d c)
              vec8d vec8d_fnmadd(vec8d a, vec8d b, vec8d c)

.. function:: vec1d vec1d_fnmsub(vec1d a, vec1d b, vec1d c)
              vec4d vec4d_fnmsub(vec4d a, vec4d b, vec4d c)
              vec8d vec8d_fnmsub(vec8d a, vec8d b, vec8d c)

.. function:: vec1d vec1d_blendv(vec1d a, vec1d b, vec1d c)
              vec4d vec4d_blendv(vec4d a, vec4d b, vec4d c)
              vec8d vec8d_blendv(vec8d a, vec8d b, vec8d c)

.. function:: vec4n vec4n_bit_shift_right(vec4n a, ulong b)
              vec8n vec8n_bit_shift_right(vec8n a, ulong b)

.. function:: vec4n vec4n_bit_and(vec4n a, vec4n b)
              vec8n vec8n_bit_and(vec8n a, vec8n b)


Modular arithmetic
-------------------------------------------------------------------------------

These functions are used internally by the small-prime FFT.
Some ``double`` variants assume an odd modulus `n < 2^{50}`.
Other assumptions are not yet documented.

.. function:: int vec1d_same_mod(vec1d a, vec1d b, vec1d n, vec1d ninv)
              int vec4d_same_mod(vec4d a, vec4d b, vec4d n, vec4d ninv)

    Return whether `a` and `b` are the same mod `n`.

.. function:: vec1d vec1d_reduce_pm1no_to_0n(vec1d a, vec1d n)
              vec1d vec4d_reduce_pm1no_to_0n(vec4d a, vec4d n)
              vec8d vec8d_reduce_pm1no_to_0n(vec8d a, vec8d n)

    Return `a \bmod n` reduced to `[0,n)` assuming `a \in (-n,n)`.

.. function:: vec1d vec1d_reduce_to_pm1n(vec1d a, vec1d n, vec1d ninv)
              vec4d vec4d_reduce_to_pm1n(vec4d a, vec4d n, vec4d ninv)
              vec8d vec8d_reduce_to_pm1n(vec8d a, vec8d n, vec8d ninv)

    Return `a \bmod n` reduced to `[-n,n]`.

.. function:: vec1d vec1d_reduce_to_pm1no(vec1d a, vec1d n, vec1d ninv)
              vec4d vec4d_reduce_to_pm1no(vec4d a, vec4d n, vec4d ninv)
              vec8d vec8d_reduce_to_pm1no(vec8d a, vec8d n, vec8d ninv)

    Return `a \bmod n` reduced to `(-n,n)`.

.. function:: vec1d vec1d_reduce_0n_to_pmhn(vec1d a, vec1d n)
              vec4d vec4d_reduce_0n_to_pmhn(vec4d a, vec4d n)

    Return `a \bmod n` reduced to `[-n/2, n/2]` given `a \in [0,n]`.

.. function:: vec1d vec1d_reduce_pm1n_to_pmhn(vec1d a, vec1d n)
              vec4d vec4d_reduce_pm1n_to_pmhn(vec4d a, vec4d n)
              vec8d vec8d_reduce_pm1n_to_pmhn(vec8d a, vec8d n)

    Return `a \bmod n` reduced to `[-n/2, n/2]` given given `a \in [-n,n]`.

.. function:: vec1d vec1d_reduce_2n_to_n(vec1d a, vec1d n)
              vec4d vec4d_reduce_2n_to_n(vec4d a, vec4d n)
              vec8d vec8d_reduce_2n_to_n(vec8d a, vec8d n)

    Return `a \bmod n` reduced to `[0,n)` given given `a \in [0,2n)`.

.. function:: vec1d vec1d_reduce_to_0n(vec1d a, vec1d n, vec1d ninv)
              vec4d vec4d_reduce_to_0n(vec4d a, vec4d n, vec4d ninv)
              vec8d vec8d_reduce_to_0n(vec8d a, vec8d n, vec8d ninv)

    Return `a \bmod n` reduced to `[0,n)`.

.. function:: vec1d vec1d_mulmod(vec1d a, vec1d b, vec1d n, vec1d ninv)
              vec4d vec4d_mulmod(vec4d a, vec4d b, vec4d n, vec4d ninv)
              vec8d vec8d_mulmod(vec8d a, vec8d b, vec8d n, vec8d ninv)

    Return `ab \bmod n` in `[-n,n]` with assumptions.

.. function:: vec1d vec1d_nmulmod(vec1d a, vec1d b, vec1d n, vec1d ninv)
              vec4d vec4d_nmulmod(vec4d a, vec4d b, vec4d n, vec4d ninv)
              vec8d vec8d_nmulmod(vec8d a, vec8d b, vec8d n, vec8d ninv)

    Return `ab \bmod n` in `[-n,n]` with assumptions.

.. function:: vec4n vec4n_addmod(vec4n a, vec4n b, vec4n n)
              vec8n vec8n_addmod(vec8n a, vec8n b, vec8n n)

    Return `a + b \bmod n` in `[0,n)`

.. function:: vec4n vec4n_addmod_limited(vec4n a, vec4n b, vec4n n)
              vec8n vec8n_addmod_limited(vec8n a, vec8n b, vec8n n)

    Return `a + b \bmod n` in `[0,n)`, assuming that `n < 2^{63}`.
