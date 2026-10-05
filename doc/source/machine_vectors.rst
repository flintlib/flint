.. _machine-vectors:

**machine_vectors.h** -- SIMD-accelerated operations on fixed-length vectors
===============================================================================

Vector types and operations mapping onto the target's SIMD instructions,
together with the ``flint_sgemm`` and ``flint_dgemm`` matrix
multiplication kernels built on top of them.

Backends are selected automatically: AVX2 (and AVX512 for the ``vec8dz``
family) on x86, NEON on ARM, and otherwise a generic backend, using GNU
vector extensions where the compiler supports them and plain ISO C
structs elsewhere. The AVX2 and NEON backends require a 64-bit word
size, since their integer vectors have :type:`ulong` lanes; a 32-bit
build uses the generic backends, which provide only the floating-point
types. The generic backends implement only the subset of the
interface required by ``flint_sgemm``/``flint_dgemm`` and, for
:type:`vec4d`, the operations used by ``dfloat`` (see
:ref:`machine-vectors-vec4d-extra`); the full interface, as used by
``fft_small``, still requires AVX2 or NEON.

For the vector operations to use the target's instructions, FLINT must
be built with appropriate compiler flags. ``configure`` chooses these
from the detected CPU; ``--enable-avx2`` and ``--enable-avx512`` force
them on.

Defining ``FLINT_MACHINE_VECTORS_FORCE_GENERIC`` before including this
header selects the generic backend even on a target with AVX2 or NEON,
and ``FLINT_MACHINE_VECTORS_STRICT_C`` additionally selects the ISO C
tier over GNU vector extensions. These are intended for testing and
profiling the portable code paths.

The strict ISO C tier has no way to express a vector operation, so
whether its operations become SIMD instructions is entirely up to the
compiler's SLP vectorizer; its performance therefore varies between
compilers and compiler versions, unlike the other backends. It exists
so that the header works on compilers without GNU vector extensions,
and is not the fallback used on GCC or clang.

The generic backends express a fused multiply-add (``fmadd``, ``fnmadd``) as ``a * b + c``,
which a compiler fuses into an FMA instruction only when floating-point
contraction is enabled. GCC in a strict ISO mode, which is how FLINT is
built, does not contract by default; code using these operations in a
performance-critical loop should request contraction, as
``machine_vectors/gemm.c`` does with ``#pragma GCC optimize
("fp-contract=fast")``. The AVX2, AVX512 and NEON backends use fused
intrinsics and are unaffected. The ``_fused`` variants below round once
on every backend.

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

.. type:: vec1f
          vec4f
          vec8f
          vec16f

    Vector with 1, 4, 8, or 16 ``float`` entries.

.. type:: vec8dz
          vec16dz
          vec16fz
          vec32fz
          vec8nz

    Vectors backed by AVX512 registers, available only when building
    with AVX512F support: 8 or 16 ``double`` entries, 16 or 32 ``float``
    entries, and 8 :type:`ulong` entries respectively. The ``z`` suffix
    distinguishes these from the equally-named types built from pairs of
    narrower registers; ``vec8d``, for instance, remains a pair of AVX2
    registers, which gives more instruction level parallelism in
    existing code.

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

.. function:: vec1n vec1d_convert_limited_vec1n(vec1d a)
              vec2n vec2d_convert_limited_vec2n(vec2d a)
              vec4n vec4d_convert_limited_vec4n(vec4d a)

    Given that each entry in the input vector is an exact integer in
    ``[0, 2^{52})``, convert it to :type:`ulong`.
    Note that ``vec1d`` and ``vec2d`` functions are only available on NEON.

.. function:: vec2d vec2n_convert_limited_vec2d(vec2n a)
              vec4d vec4n_convert_limited_vec4d(vec4n a)
              vec8d vec8n_convert_limited_vec8d(vec8n a)

    The inverse conversion, from :type:`ulong` entries to exact
    ``double`` entries. Requires the same assumption.

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

    Entrywise comparisons (see below for the full set).

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

.. _machine-vectors-vec4d-extra:

Masks, bits, fused operations and gathers for vec4d
-------------------------------------------------------------------------------

These are provided by every backend (AVX2, NEON and both generic
tiers) with the same results, and are tested lane by lane against
scalar code in ``machine_vectors/test/t-vec4d.c``.

.. function:: vec4d vec4d_cmp_eq(vec4d a, vec4d b)
              vec4d vec4d_cmp_ne(vec4d a, vec4d b)
              vec4d vec4d_cmp_lt(vec4d a, vec4d b)
              vec4d vec4d_cmp_le(vec4d a, vec4d b)
              vec4d vec4d_cmp_nlt(vec4d a, vec4d b)
              vec4d vec4d_cmp_nle(vec4d a, vec4d b)
              vec4d vec4d_cmp_ngt(vec4d a, vec4d b)
              vec4d vec4d_cmp_nge(vec4d a, vec4d b)

    With ``vec4d_cmp_gt`` and ``vec4d_cmp_ge``: a mask with all bits
    set in the lanes where the predicate holds and all bits clear
    elsewhere. The predicates ``eq``, ``lt``, ``le``, ``gt``, ``ge``
    are false when a lane of *a* or *b* is a NaN, and their negations
    ``ne``, ``nlt``, ``nle``, ``ngt``, ``nge`` are true.

.. function:: int vec4d_movemask(vec4d a)

    The sign bits of the lanes, lane `i` giving bit `i`. (For
    ``vec4d_blendv``, which selects *b* where the sign bit of *c* is
    set, a mask from a comparison selects *b* where it holds.)

.. function:: vec4d vec4d_bit_and(vec4d a, vec4d b)
              vec4d vec4d_bit_or(vec4d a, vec4d b)
              vec4d vec4d_bit_xor(vec4d a, vec4d b)
              vec4d vec4d_bit_andnot(vec4d a, vec4d b)

    Bitwise operations; ``andnot`` computes ``(~a) & b``.

.. function:: vec4d vec4d_bits_set_u64(uint64_t a)
              vec4d vec4d_bits_add(vec4d a, vec4d b)
              vec4d vec4d_bits_sub(vec4d a, vec4d b)
              vec4d vec4d_bits_shl(vec4d a, int n)
              vec4d vec4d_bits_shr(vec4d a, int n)

    Operations on the lanes viewed as ``uint64_t`` bit patterns: a
    broadcast, addition and subtraction modulo `2^{64}`, and logical
    shifts by `0 \le n < 64`.

.. function:: vec4d vec4d_sqrt(vec4d a)

    The correctly rounded square root.

.. function:: vec4d vec4d_fmadd_fused(vec4d a, vec4d b, vec4d c)
              vec4d vec4d_fmsub_fused(vec4d a, vec4d b, vec4d c)
              vec4d vec4d_fnmadd_fused(vec4d a, vec4d b, vec4d c)

    `ab + c`, `ab - c` and `-ab + c` with a single rounding, like
    ``fma(a, b, c)``, ``fma(a, b, -c)`` and ``fma(-a, b, c)``
    (including the sign of an exact zero).

.. function:: vec4d vec4d_permute_0_1_1_0(vec4d a)
              vec4d vec4d_permute_0_2_0_2(vec4d a)
              vec4d vec4d_permute_1_3_1_3(vec4d a)

    Further permutations, like ``vec4d_permute_0_2_1_3``.

.. type:: vec4di

.. function:: vec4di vec4d_to_index(vec4d a)
              vec4d vec4d_gather(const double * t, vec4di i)

    Table indices converted from lanes holding integers in `[0, 2^{31})`,
    and the lanes ``t[i0], t[i1], t[i2], t[i3]`` at those indices.

In the generic backends, ``vec4d_min`` and ``vec4d_max`` return *b*
where a lane is unordered, as the AVX2 instructions do; NEON returns a
NaN.


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

    Return `a \bmod n` reduced to `[-n/2, n/2]` given `a \in [-n,n]`.

.. function:: vec1d vec1d_reduce_2n_to_n(vec1d a, vec1d n)
              vec4d vec4d_reduce_2n_to_n(vec4d a, vec4d n)
              vec8d vec8d_reduce_2n_to_n(vec8d a, vec8d n)

    Return `a \bmod n` reduced to `[0,n)` given `a \in [0,2n)`.

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

Matrix multiplication
-------------------------------------------------------------------------------

These functions compute a matrix product in single or double precision.
They are always available: when FLINT is built with BLAS, the default is
to call it, and otherwise FLINT's own kernels are used. The intended use
is as a building block for exact linear algebra over `\mathbb{Z}` and
`\mathbb{Z}/n\mathbb{Z}`, for example in :func:`nmod_mat_mul_blas` and
:func:`fmpz_mat_mul_blas`.

All of these functions compute `C = AB` for row-major matrices, where
*C* is *m* by *n*, *A* is *m* by *k* and *B* is *k* by *n*, with
*ldc*, *lda* and *ldb* the respective leading dimensions (the number of
entries between the start of consecutive rows, which must be at least
the number of columns). No transposition or accumulation is performed:
the previous contents of *C* are overwritten, and `k = 0` sets *C* to
zero. This is equivalent to ``cblas_sgemm`` or ``cblas_dgemm`` called
with ``CblasRowMajor``, ``CblasNoTrans``, ``CblasNoTrans``, ``alpha``
equal to 1 and ``beta`` equal to 0. Aliasing of *C* with *A* or *B* is
not allowed.

The FLINT kernels are multithreaded internally according to
:func:`flint_get_num_threads`, using FLINT's thread pool. They handle
arbitrary dimensions, including thin and unbalanced shapes, without
requiring any padding or alignment of the input.

.. function:: void flint_sgemm(slong m, slong k, slong n, const float * A, slong lda, const float * B, slong ldb, float * C, slong ldc)
              void flint_dgemm(slong m, slong k, slong n, const double * A, slong lda, const double * B, slong ldb, double * C, slong ldc)

    Sets `C = AB`, calling either the BLAS or the FLINT implementation
    according to :var:`flint_gemm_use_blas`.

.. function:: void flint_sgemm_blas(slong m, slong k, slong n, const float * A, slong lda, const float * B, slong ldb, float * C, slong ldc)
              void flint_dgemm_blas(slong m, slong k, slong n, const double * A, slong lda, const double * B, slong ldb, double * C, slong ldc)

    Sets `C = AB` using ``cblas_sgemm`` or ``cblas_dgemm``. These raise
    an exception if FLINT was built without BLAS support.

.. function:: void flint_sgemm_fallback(slong m, slong k, slong n, const float * A, slong lda, const float * B, slong ldb, float * C, slong ldc)
              void flint_dgemm_fallback(slong m, slong k, slong n, const double * A, slong lda, const double * B, slong ldb, double * C, slong ldc)

    Sets `C = AB` using FLINT's own kernels, which are always available.

.. var:: int flint_gemm_use_blas

    Selects the implementation used by :func:`flint_sgemm` and
    :func:`flint_dgemm`. It is initialized to 1 if FLINT was built with
    BLAS support and 0 otherwise, and may be set to either value at
    runtime, for example to compare the two implementations. Setting it
    to 1 in a build without BLAS support will result in an exception
    when a gemm is attempted.
