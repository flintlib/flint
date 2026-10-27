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
# define NMOD_MAT_HAVE_VECMAT_U32 1
#else
# define NMOD_MAT_HAVE_VECMAT_U32 0
#endif

/*
    Moduli up to 2^32 on x86 without IFMA: whether the rows path takes the
    floating point products of the moduli below 2^50 (fp50, when there is a
    vector FMA) rather than the u32 path, for vectors of length len and B
    with ncols columns: u32 reduces each entry of the result with scalar
    two-word reductions, which dominate for short vectors, whereas fp50
    reduces them in vectors.
*/
#define NMOD_MAT_VECMAT_U32_TO_FP50(len, ncols) \
    (NMOD_MAT_HAVE_FPV \
     && (len) <= (((ncols) >= 16) ? FLINT_MIN((ncols), 32) : 6))

/*
    Whether the vector-matrix products of nmod_vec_mul.c have their u64
    path, for the moduli the paths above leave (on x86 with AVX2 or
    AVX-512: n >= 2^50, n > 2^52 with IFMA, and 2^32 < n < 2^50 without a
    vector FMA): the entries are split into two limbs of at most 32 bits
    whose vpmuludq products are accumulated in several lanes, as in the
    split limbs dot products of nmod_vec/dot_split_limbs.c. Then every
    modulus has a SIMD path.
*/
#if FLINT_BITS == 64 && defined(__AVX2__) \
        && !defined(FLINT_MACHINE_VECTORS_FORCE_GENERIC) \
        && !defined(FLINT_MACHINE_VECTORS_STRICT_C)
# define NMOD_MAT_HAVE_VECMAT_U64 1
#else
# define NMOD_MAT_HAVE_VECMAT_U64 0
#endif

/* whether _nmod_mat_mul_rows_simd has a SIMD path for this modulus (which
   nmod_mat_nmod_vec_mul uses when B is wide enough), rather than one modular
   multiplication per entry; _IS_U64: whether this is the u64 path */
#define _NMOD_MAT_NMOD_VEC_MUL_IS_SIMD_SMALL(n) \
    ((NMOD_MAT_HAVE_MUL_U52 && (n) <= (UWORD(1) << 52)) \
     || (NMOD_MAT_HAVE_VECMAT_U32 && (n) <= (UWORD(1) << 32)) \
     || (NMOD_MAT_HAVE_FPV && (n) < (UWORD(1) << 50)))
#define NMOD_MAT_NMOD_VEC_MUL_IS_SIMD(n) \
    (_NMOD_MAT_NMOD_VEC_MUL_IS_SIMD_SMALL(n) || NMOD_MAT_HAVE_VECMAT_U64)
#define NMOD_MAT_NMOD_VEC_MUL_IS_U64(n) \
    (NMOD_MAT_HAVE_VECMAT_U64 && !_NMOD_MAT_NMOD_VEC_MUL_IS_SIMD_SMALL(n))

/*
    The fewest columns of B for which nmod_mat_mul uses the u64 path of
    _nmod_mat_mul_rows_simd with 2 rows or more (one row goes through
    nmod_mat_nmod_vec_mul and its own rule; see also _nmod_mat_mul_rows_wanted
    in mul.c for 4 and 5 columns): from 6 columns on it was faster
    than nmod_mat_mul_classical for the inner dimensions k = 2..512 and 2 to
    4 rows (Zen 4, Ice Lake SP, Cascade Lake, 51 to 64 bits), and than what
    nmod_mat_mul did with 6 to 15 columns (also nmod_mat_mul_k52 at 51-52
    bits without IFMA) by 1.1-2.4x, except at 64 bits with 4 rows and
    k = 512 (up to 1.3x slower, Ice Lake SP); with 4 columns it won only
    up to k = 32, and it lost with 2.
*/
#define NMOD_MAT_MUL_ROWS_U64_MIN_COLS 6

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

/* whether nmod_mat_mul has SIMD kernels for this modulus (mul_u32 below
   2^32 on any 64-bit target, mul_u52 with IFMA, mul_fp50 or mul_k52 with
   a vector FMA, up to 2^52), which use several threads */
#if FLINT_BITS == 64
# define NMOD_MAT_MUL_HAVE_SIMD_KERNEL(n) \
    ((n) <= (UWORD(1) << 32) \
     || ((n) <= (UWORD(1) << 52) && (NMOD_MAT_HAVE_MUL_U52 || NMOD_MAT_HAVE_FPV)))
#else
# define NMOD_MAT_MUL_HAVE_SIMD_KERNEL(n) 0
#endif

/* whether, for this modulus, nmod_mat_mul prefers the SIMD kernels to
   _nmod_mat_mul_rows_simd for the shapes both can take: moduli below 2^32
   when the rows engine only has its floating point tier for them (no
   IFMA, no x86 integer tier), where nmod_mat_mul_u32 is faster */
#if !NMOD_MAT_HAVE_MUL_U52 && !NMOD_MAT_HAVE_VECMAT_U32
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
    of them, with the delayed reductions of nmod_mat_nmod_vec_mul (which
    uses the case r = 1 when B is wide enough for it, see
    vecmat_use_rows_simd in nmod_vec_mul.c). Returns 0 without doing
    anything when the modulus has
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

/*
    C = A * B for tiny matrices, all dimensions at most 8 (mul.c): each
    entry is a dot product accumulated in one, two or three words and
    reduced once. nmod_mat_mul uses it when all dimensions are at most
    NMOD_MAT_MUL_TINY_MAX. C may alias A or B.
*/
#define NMOD_MAT_MUL_TINY_MAX 3
void _nmod_mat_mul_tiny(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B);

/*
    PLE decomposition (ple*.c). The internal functions take
    Q = pivots_nonpivots (not NULL) and a bound max_nonpivots on the number
    of columns without pivot met while some rows are not yet pivot rows:
    if it is nonnegative and exceeded, they stop and return -1, leaving P,
    Q and A undefined. A negative bound means no limit. _nmod_mat_ple
    chooses the algorithm.
*/
slong _nmod_mat_ple(slong * P, slong * Q, nmod_mat_t A, slong max_nonpivots);
slong _nmod_mat_ple_classical(slong * P, slong * Q, nmod_mat_t A,
                              slong max_nonpivots);
slong _nmod_mat_ple_classical_delayed(slong * P, slong * Q, nmod_mat_t A,
                                      slong max_nonpivots);
slong _nmod_mat_ple_crout(slong * P, slong * Q, nmod_mat_t A,
                          slong max_nonpivots);
slong _nmod_mat_ple_recursive(slong * P, slong * Q, nmod_mat_t A,
                              slong max_nonpivots);

/* the recursion of _nmod_mat_ple_recursive, whose subproblems with
   min(m, n) >= cutoff recurse again and the others use the delayed
   algorithm; with cutoff <= 0, the subproblems use _nmod_mat_ple (this
   is _nmod_mat_ple_recursive). Small cutoffs are for the tests. */
slong _nmod_mat_ple_recursive_cutoff(slong * P, slong * Q, nmod_mat_t A,
                                     slong max_nonpivots, slong cutoff);

/* upward rotation: moves P[r] to position k and P[k..r-1] down to
   P[k+1..r] (k <= r) */
static inline void
_nmod_mat_ple_rotate(slong * P, slong k, slong r)
{
    slong i, t = P[r];

    for (i = r; i > k; i--)
        P[i] = P[i - 1];
    P[k] = t;
}

/*
    End of the base cases that swap rows (delayed, Crout). These store
    the pivot row k in row k of A, and the other rows in the rows rank, ...,
    m - 1 of A in any order, orig[i] being the index in the input of the
    row stored in row i of A. The pivot of a column is the row with a
    nonzero entry of smallest index in the input among the rows that are
    not pivot rows yet, which is the first one in the logical order of
    ple_classical.c (these rows are in increasing order there), and P[k] =
    orig[k] for k < rank. This function sets P[rank..m) to the indices of
    the other rows in increasing order, puts these rows in this order (orig
    is destroyed), and appends the columns without pivot to Q in
    increasing order.
*/
void _nmod_mat_ple_rows_finish(slong * P, slong * Q, nmod_mat_t A,
                               slong * orig, slong rank);

/* new row i of A = old row perm[i], in place (fixed rows not touched) */
void _nmod_mat_ple_permute_rows(nmod_mat_t A, const slong * perm);

/* end of the base cases: rows put in the order of P, columns without
   pivot appended to Q in increasing order */
void _nmod_mat_ple_finish(slong * P, slong * Q, nmod_mat_t A, slong rank);

/* the public interface of an algorithm: flags, NULL pivots_nonpivots */
slong _nmod_mat_ple_wrap(slong * P, slong * pivots_nonpivots, nmod_mat_t A,
        int flags, slong (* ple)(slong *, slong *, nmod_mat_t, slong));

/*
    Whether nmod_mat_addmul and nmod_mat_submul (single-threaded, below
    the Strassen cutoff) should compute A * B with nmod_mat_mul and add
    it, rather than with the fused classical _nmod_mat_mul_classical_op:
    when nmod_mat_mul has a faster kernel than the classical code for the
    modulus (mul_u8; with 64-bit words mul_u32, mul_u52,
    mul_fp50, mul_k52 and the row and column engines), except for tiny
    products and for B with fewer than 4 columns, where the temporary
    matrix and the extra pass cost as much as they save.
*/
static inline int
_nmod_mat_addmul_want_mul(slong m, slong k, slong n, nmod_t mod)
{
    if (n < 4 || m * k * FLINT_MIN(n, 128) < 128)
        return 0;

    /* nmod_mat_mul_u8, which does not pay for its conversions below
       dimension 16 */
    if (mod.n <= 255 && FLINT_MIN(m, k) >= 8 && n >= 8)
        return FLINT_MIN(m, k) >= 16 && n >= 16;

#if FLINT_BITS == 64
    if (mod.n <= (UWORD(1) << 32))
        return 1;

    if (mod.n <= (UWORD(1) << 52) && (NMOD_MAT_HAVE_MUL_U52 || NMOD_MAT_HAVE_FPV))
        return 1;
#endif

    return 0;
}

/*
    Triangular solves by substitution with delayed reductions
    (solve_tri_delayed.c), for an r x r triangular T (lower or upper) and B
    with k columns (T X = B, right = 0) or k rows (X T = B, right = 1).
    X may alias B; the diagonal of T is not read if unit is nonzero.

    _nmod_mat_solve_tri_rows_delayed computes X (or X^T, B being transposed
    by blocks) by rows, each as a combination of the previous ones: vectors
    of length k. _nmod_mat_solve_tri_right_axpy_delayed solves X T = B row
    by row, adding x[j] times the row j of T to the remaining entries:
    vectors of length at most r - 1.

    The choices between these kernels, the dot products and the recursive
    solves are below (_nmod_mat_solve_tri_param).
*/
void _nmod_mat_solve_tri_rows_delayed(nmod_mat_t X, const nmod_mat_t T,
        const nmod_mat_t B, int unit, int upper, int right);
void _nmod_mat_solve_tri_right_axpy_delayed(nmod_mat_t X, const nmod_mat_t T,
        const nmod_mat_t B, int unit, int upper);

/*
    Choices of the triangular solves, for a triangular matrix of dimension r
    and k columns (left) or rows (right) of B, by class of machine (AVX-512
    IFMA, AArch64 NEON, other) and range of the modulus n: n <= 3 * 2^29,
    n <= 2^31, n < 2^32, n <= 2^52, above (accumulations of the delayed
    kernels, dot products and matrix products differ between them), and
    n = 2 (last entry: recursion behaves differently there).
    - left classical: by rows of X if k >= max(lk, r / lc), otherwise one dot
      product per entry;
    - right classical: by rows of X^T if k >= max(rk, r / rc), otherwise axpy
      if r <= ar or k <= ak, otherwise dot products;
    - top level (nmod_mat_solve_tril etc.): recursive if r >= LAr and
      k >= LAk, or if r >= Lr, k >= Lk and r k >= Lp (left); r >= RAr and
      k >= RAk, or r >= Rr, k >= Rk and r k >= Rp (right); otherwise
      classical.
    Fitted on uniformly random inputs (Zen 4 and Ice Lake for IFMA, Apple
    M4 for NEON, Cascade Lake for the others); the recursion from timings
    of the classical solve and of one level of recursion (r = 16..1024,
    k = 4..1024), then checked against the previous choice with the whole
    recursion, keeping the previous values where they did better.
*/
typedef struct
{
    short lk, lc, rk, rc, ar, ak, LAr, LAk, Lr, Lk, Lp, RAr, RAk, Rr, Rk, Rp;
}
_nmod_mat_solve_tri_param_struct;

#define _TRI_NEVER 32767

#if NMOD_MAT_HAVE_MUL_U52
# define _NMOD_MAT_SOLVE_TRI_PARAMS \
    {{3, 6, 6, 4, 256, 32, 48, 48, 48, 12, 512, 64, 64, 64, 1, 2048}, \
     {3, 4, 6, 2, 128, 8, 64, 64, 24, 12, 0, 64, 64, 32, 1, 96}, \
     {3, 2, 6, 2, 64, 4, 64, 64, 8, 12, 128, 64, 64, 24, 1, 64}, \
     {3, 3, 4, 2, 64, 4, 64, 64, 8, 8, 128, 64, 64, 24, 1, 64}, \
     {3, 3, 4, 2, 64, 4, 64, 64, 768, 24, 0, 128, 128, 512, 1, 0}, \
     {3, 6, 6, 4, 256, 32, 64, 64, 48, 12, 0, 64, 64, 64, 24, 0}}
#elif (defined(__ARM_NEON) && defined(__aarch64__)) || defined(_M_ARM64)
# define _NMOD_MAT_SOLVE_TRI_PARAMS \
    {{2, 12, 4, 2, 16, 63, 64, 64, 32, 8, 512, 64, 64, 48, 16, 0}, \
     {2, 12, 4, 2, 16, 63, 64, 64, 16, 8, 0, 32, 32, 32, 8, 512}, \
     {2, 8, 4, 4, 16, 63, 64, 64, 24, 8, 0, 32, 32, 32, 8, 256}, \
     {3, 8, 4, 4, 256, 16, 64, 64, 384, 8, 0, 64, 64, 256, 12, 0}, \
     {3, 12, 4, 8, 256, 16, 256, 256, 512, 24, 0, _TRI_NEVER, _TRI_NEVER, 384, 12, 8192}, \
     {2, 12, 4, 2, 16, 63, 64, 64, 32, 8, 512, 64, 64, 48, 16, 0}}
#else
# define _NMOD_MAT_SOLVE_TRI_PARAMS \
    {{4, 4, 6, 4, 128, 16, 48, 48, 48, 4, 4096, 64, 24, 512, 1, 0}, \
     {6, 4, 6, 2, 64, 4, 64, 64, 24, 24, 1024, 32, 32, 32, 1, 512}, \
     {4, 2, 4, 1, 32, 4, 48, 32, 512, 1, 1024, 64, 32, 256, 2, 0}, \
     {3, 3, 4, 2, 256, 8, 32, 32, 32, 8, 256, 32, 32, 32, 1, 256}, \
     {3, 3, 4, 2, 64, 8, 64, 64, _TRI_NEVER, _TRI_NEVER, 0, 256, 256, 512, 8, 0}, \
     {4, 4, 6, 4, 128, 16, 64, 64, 384, 2, 1024, 64, 64, 64, 24, 0}}
#endif

FLINT_FORCE_INLINE const _nmod_mat_solve_tri_param_struct *
_nmod_mat_solve_tri_param(ulong n)
{
    static const _nmod_mat_solve_tri_param_struct tab[6] =
        _NMOD_MAT_SOLVE_TRI_PARAMS;
#if FLINT_BITS == 64
    return tab + ((n == 2) ? 5 : (n <= (UWORD(3) << 29)) ? 0 : (n <= (UWORD(1) << 31)) ? 1
                  : (n < (UWORD(1) << 32)) ? 2 : (n <= (UWORD(1) << 52)) ? 3 : 4);
#else
    return tab + ((n == 2) ? 5 : (n <= (UWORD(3) << 13)) ? 0 : (n <= (UWORD(1) << 15)) ? 1 : 2);
#endif
}

#define NMOD_MAT_SOLVE_TRI_LEFT_ROWS(r, k, n) \
    ((k) >= FLINT_MAX(_nmod_mat_solve_tri_param(n)->lk, \
                      (r) / _nmod_mat_solve_tri_param(n)->lc))
#define NMOD_MAT_SOLVE_TRI_RIGHT_ROWS(r, k, n) \
    ((k) >= FLINT_MAX(_nmod_mat_solve_tri_param(n)->rk, \
                      (r) / _nmod_mat_solve_tri_param(n)->rc))
#define NMOD_MAT_SOLVE_TRI_RIGHT_AXPY(r, k, n) \
    ((r) <= _nmod_mat_solve_tri_param(n)->ar \
     || (k) <= _nmod_mat_solve_tri_param(n)->ak)
#define NMOD_MAT_SOLVE_TRI_LEFT_REC(r, k, n) \
    (((r) >= _nmod_mat_solve_tri_param(n)->LAr \
      && (k) >= _nmod_mat_solve_tri_param(n)->LAk) \
     || ((r) >= _nmod_mat_solve_tri_param(n)->Lr \
         && (k) >= _nmod_mat_solve_tri_param(n)->Lk \
         && (r) * (k) >= _nmod_mat_solve_tri_param(n)->Lp))
#define NMOD_MAT_SOLVE_TRI_RIGHT_REC(r, k, n) \
    (((r) >= _nmod_mat_solve_tri_param(n)->RAr \
      && (k) >= _nmod_mat_solve_tri_param(n)->RAk) \
     || ((r) >= _nmod_mat_solve_tri_param(n)->Rr \
         && (k) >= _nmod_mat_solve_tri_param(n)->Rk \
         && (r) * (k) >= _nmod_mat_solve_tri_param(n)->Rp))

#endif
