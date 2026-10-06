/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MAT_IMPL_H
#define NMOD_MAT_IMPL_H

#include "flint.h"
#include "flint-mparam.h"  /* FLINT_NMOD_MAT_MUL_FP50_U32_{MIN,MAX}_K */

/*
    Whether nmod_mat_mul_u52 has a kernel: it needs AVX512-IFMA (with F and
    DQ) at compile time and 64-bit words; otherwise it declines every
    multiplication. Both mul_u52.c and the dispatch in mul.c test this.
*/
#if FLINT_BITS == 64 && defined(__AVX512F__) && defined(__AVX512DQ__) \
        && defined(__AVX512IFMA__) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# define NMOD_MAT_HAVE_MUL_U52 1
#else
# define NMOD_MAT_HAVE_MUL_U52 0
#endif

/*
    Whether the double precision primitives of mul_fp_vec.h have a vector
    backend (AVX-512, AVX2 + FMA, AArch64 NEON), on which nmod_mat_mul_fp50
    and the delayed-reduction nmod_mat_nmod_vec_mul rely; mirrors the
    selection in that header.
*/
#if FLINT_BITS == 64 \
        && ((defined(__AVX512F__) && defined(__AVX512DQ__)) \
            || (defined(__AVX2__) && (defined(__FMA__) || defined(_MSC_VER))) \
            || (defined(__ARM_NEON) && defined(__aarch64__)) \
            || defined(_M_ARM64)) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C) \
        && !defined(NMOD_MAT_FP_FORCE_GENERIC)
# define NMOD_MAT_HAVE_FPV 1
#else
# define NMOD_MAT_HAVE_FPV 0
#endif

/*
    Whether the vector-matrix products of nmod_vec_mul.c have their
    integer path for moduli up to 2^32 on x86 without IFMA (vpmuludq
    products split at bit 56, as in the dot products of nmod_vec/dot.c).
*/
#if FLINT_BITS == 64 && !NMOD_MAT_HAVE_MUL_U52 \
        && ((defined(__AVX512F__) && defined(__AVX512DQ__)) \
            || defined(__AVX2__)) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# define NMOD_MAT_HAVE_VM_U32 1
#else
# define NMOD_MAT_HAVE_VM_U32 0
#endif

/* whether nmod_mat_nmod_vec_mul and _nmod_mat_mul_rows_simd delay their
   reductions for this modulus (SIMD accumulation), rather than one modular
   multiplication per entry */
#define NMOD_MAT_NMOD_VEC_MUL_IS_SIMD(n) \
    ((NMOD_MAT_HAVE_MUL_U52 && (n) <= (UWORD(1) << 52)) \
     || (NMOD_MAT_HAVE_VM_U32 && (n) <= (UWORD(1) << 32)) \
     || (NMOD_MAT_HAVE_FPV && (n) < (UWORD(1) << 50)))

/*
    Whether the few-columns products of mul_cols.c have a vectorized path
    for this modulus, and the most columns of B for which nmod_mat_mul
    uses them: n <= 2^52 on AVX512-IFMA (u52; u64 when the inner dimension
    is too large for u52), 2^32 < n < 2^50 with a vector backend of
    mul_fp_vec.h otherwise (other moduli go to mul_classical: transposed
    columns + nmod_vec_dot). Above 2^52, the u64 tier of mul_cols.c is not
    faster than mul_classical, whose dot products are the split-limbs / u64
    ones.
*/
#if NMOD_MAT_HAVE_MUL_U52
# define NMOD_MAT_MUL_COLS_IS_SIMD(n) ((n) <= (UWORD(1) << 52))
# define NMOD_MAT_MUL_COLS_MAX 8
#else
# define NMOD_MAT_MUL_COLS_IS_SIMD(n) \
    (NMOD_MAT_HAVE_FPV && (n) > (UWORD(1) << 32) && (n) < (UWORD(1) << 50))
# define NMOD_MAT_MUL_COLS_MAX 4
#endif

/* most rows of A for which nmod_mat_mul uses _nmod_mat_mul_rows_simd */
#if NMOD_MAT_HAVE_MUL_U52
# define NMOD_MAT_MUL_ROWS_MAX 8
#else
# define NMOD_MAT_MUL_ROWS_MAX 4
#endif

/* whether, for this modulus, nmod_mat_mul prefers the SIMD kernels to
   _nmod_mat_mul_rows_simd for the shapes both can take: moduli below 2^32
   when the rows engine only has its floating point tier for them (no
   IFMA, no x86 integer tier), where nmod_mat_mul_u32 is faster */
#if !NMOD_MAT_HAVE_MUL_U52 && !NMOD_MAT_HAVE_VM_U32
# define NMOD_MAT_MUL_ROWS_PREFER_KERNEL(n) ((n) <= (UWORD(1) << 32))
#else
# define NMOD_MAT_MUL_ROWS_PREFER_KERNEL(n) 0
#endif

/*
    Inner dimensions k for which nmod_mat_mul prefers nmod_mat_mul_fp50 to
    nmod_mat_mul_u32 for moduli up to 2^32 (flint-mparam.h parameters
    FLINT_NMOD_MAT_MUL_FP50_U32_MIN_K / _MAX_K, empty range if MIN > MAX):
    with so few products per entry the per-entry work of u32 dominates (it
    costs the same at 20 and 31 bits). 1000 x k x 1000, fp50 / u32 at
    k = 1, 2, 3: 0.5-0.9 on the AVX2 machines measured (Broadwell to Arrow
    Lake); with AVX-512 0.61-0.94 on Zen 4, 0.77-0.82 on Cascade Lake, 0.84,
    0.79, 0.88 on Ice Lake (Xeon), but 1.07, 0.73, 0.84 on Emerald Rapids
    and 0.77, 1.09, 1.39 on Tiger Lake (a single 512-bit FMA unit); 0.82,
    1.06, 1.27 on Apple M4. Typical values: k = 1..3 on AVX2, k = 1 with
    AVX-512 or NEON. Not used without vector code (no NMOD_MAT_HAVE_FPV).
*/

/* the parameters of nmod_mat_mul, which every flint-mparam.h must define */
#if !defined(FLINT_NMOD_MAT_MUL_SIMD_MIN_DIM) \
    || !defined(FLINT_NMOD_MAT_MUL_BLAS_1PASS_CUTOFF) \
    || !defined(FLINT_NMOD_MAT_MUL_BLAS_1PASS_CUTOFF_MT) \
    || !defined(FLINT_NMOD_MAT_MUL_SIMD_STRASSEN_CUTOFF) \
    || !defined(FLINT_NMOD_MAT_MUL_U52_MIN_BITS) \
    || !defined(FLINT_NMOD_MAT_MUL_U52_LO_MAX_BITS) \
    || !defined(FLINT_NMOD_MAT_MUL_U52_MIN_K) \
    || !defined(FLINT_NMOD_MAT_MUL_K52_MIN_BITS) \
    || !defined(FLINT_NMOD_MAT_MUL_K52_BLAS_CUTOFF) \
    || !defined(FLINT_NMOD_MAT_MUL_FP50_MAX_BITS) \
    || !defined(FLINT_NMOD_MAT_MUL_FP50_U32_MIN_K) \
    || !defined(FLINT_NMOD_MAT_MUL_FP50_U32_MAX_K)
# error "flint-mparam.h must define all FLINT_NMOD_MAT_MUL_* parameters (see src/mpn_extras/generic/flint-mparam.h)"
#endif

#include "nmod_types.h"

/*
    nmod_mat_mul_strassen, except that the products of the recursion whose
    dimensions are all at least cutoff (if cutoff > 0) use Strassen again
    instead of going back to nmod_mat_mul. This is for the tests, which use
    small cutoffs to exercise several levels and all parities of the
    dimensions.
*/
void _nmod_mat_mul_strassen_cutoff(nmod_mat_t C, const nmod_mat_t A,
                                   const nmod_mat_t B, slong cutoff);

/*
    c[s] = a[s] * B for 1 <= r <= NMOD_MAT_MUL_ROWS_MAX vectors a[s] of
    length len (1 <= len <= B->r), the rows of B being read once for all
    of them, with the delayed reductions of nmod_mat_nmod_vec_mul (which is
    the case r = 1). Returns 0 without doing anything when the modulus has
    no SIMD path (NMOD_MAT_NMOD_VEC_MUL_IS_SIMD). The c[s] must not overlap
    the a[s] nor B.
*/
int _nmod_mat_mul_rows_simd(ulong * const * c, const ulong * const * a,
                            slong r, slong len, const nmod_mat_t B);

/*
    C = A * B for B with 1 <= B->c <= NMOD_MAT_MUL_COLS_MAX columns and
    A->c >= 1, each row of A multiplied with all the columns of B at once
    (mul_cols.c). Returns 0 without doing anything when the modulus has no
    SIMD path (NMOD_MAT_MUL_COLS_IS_SIMD). C must not alias A or B.
*/
int _nmod_mat_mul_cols_simd(nmod_mat_t C, const nmod_mat_t A,
                            const nmod_mat_t B);

#endif
