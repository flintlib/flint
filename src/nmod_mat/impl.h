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

/* most rows of A for which nmod_mat_mul uses _nmod_mat_mul_rows_simd */
#if NMOD_MAT_HAVE_MUL_U52
# define NMOD_MAT_MUL_ROWS_MAX 8
#else
# define NMOD_MAT_MUL_ROWS_MAX 4
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

#endif
