/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Matrix multiplication over Z/nZ for n < 2^32 with integer SIMD.

    nmod_mat_mul_blas lifts the entries to floating point and relies on the
    exactness of a dot product of length k, which for a single double precision
    gemm requires k*(n/2)^2 < 2^53. Past roughly 20-bit moduli this forces a
    multimodular approach: several dgemm calls and a CRT reconstruction, which
    is what makes moduli between about 2^24 and 2^32 a weak spot of this
    strategy since there one may rather rely on integer SIMD instructions.

    This file instead uses the widening 32x32 -> 64 bit multiplication that
    every integer SIMD instruction set provides (vpmuldq on AVX2 and AVX-512,
    smlal on NEON). Entries are lifted to the symmetric range |a| <= n/2 as
    int32, products |ab| <= (n/2)^2 are accumulated exactly in signed int64
    lanes, and the accumulators are reduced only when they are about to
    overflow. This is one pass over the data whatever the modulus below 2^32,
    and no CRT.

    Delayed reduction. With P = floor(n/2)^2 the accumulator can absorb about
    2^63/P products before a reduction is needed, from 2 at 32-bit moduli
    through 8 at 31 bits, 32 at 30 bits, 128 at 29 bits, to more than a block's
    worth of k below 27 bits (where the microkernel then runs with no reduction
    at all). The in-kernel reduction is not towards the canonical remainder:
    writing acc = hi*2^32 + lo with hi the signed high word, acc is congruent
    to hi*c + lo where c = 2^32 mod n, and this replacement costs a shift, a
    mask, one more widening multiply and an add per accumulator. One round
    shrinks |acc| from 2^63 to about 2^31*c, a second (when c is large) to
    about n^2/2. The exact bounds are computed from n in u32_ctx_init and the
    cadence follows from them, so the code is correct for every n < 2^32 with
    no case analysis.

    A reduction to the canonical remainder is done at the end of each k-block,
    when a tile is written back to the output C: the block partial sums live in
    C itself as canonical residues, so there is no floating point or 64-bit
    temporary for C. That reduction is a Barrett step through double precision
    (measured at 1-5% of the total, since it runs once per output element per
    k-block rather than once per product).
    TODO : alternatives with an integer Barrett could be considered, see the
    discussion in
    https://github.com/flintlib/flint/issues/2597
    The straightforward version with a precomputed floor(2^64/n) is not much
    interesting as it would need a 64x64 -> high multiply, which none of AVX2,
    AVX-512 (without IFMA) or NEON provides. This is not critical since this
    is generally a negligible step of the whole process.

    Structure: BLIS-style blocking (NC/KC/MC) with packed A blocks and B panels
    of int32, a register tile microkernel of MR rows by NACC accumulator
    vectors, and for multithreading a plain split of C into independent row or
    column blocks over FLINT's thread pool. All of this is the template
    mul_blocked_templ.h (shared with the AVX512-IFMA kernel of mul_u52.c),
    instantiated here twice: for ulong entries, behind nmod_mat_mul_u32, and
    for uint32 entries, behind _nmod_mat_mul_u32 (used by the nmod32 ring of
    gr, and the natural input for multimodular algorithms that hold residues
    as uint32). The kernels never see the entry type: it only enters the
    packing and the loads and stores of C.

    Backends: the template is written against a small set of inline
    primitives (load, splat, widening multiply-accumulate, fold, canonical
    reduce) that each backend defines: AVX-512 (F+DQ), AVX2, AArch64 NEON,
    and a plain C fallback whose "vector" is a small struct the compiler is
    free to vectorize. All four are covered by the same test.

    TODO investigations for crafting these BLAS-like functions suggest possible
    changes to machine_vectors.h. In detail, the primitives are local to this
    file rather than taken from machine_vectors.h because that header has no
    signed widening multiply on any ISA (its vec4n_mul is _mm256_mul_epu32,
    AVX2 only), no arithmetic 64-bit shift and no 64-bit compare; and because
    the natural operand shape differs between the two families -- x86 reads the
    even 32-bit halves of 64-bit lanes, NEON takes narrow 2x32 operands -- so a
    common "vector of 4 ulong" abstraction would force a shuffle per operand on
    one side or the other. Promoting these to machine_vectors.h as a
    signed-integer tier is a worthwhile separate change. The header is still
    included so that this file follows FLINT's ISA detection and honours
    FLINT_MACHINE_VECTORS_FORCE_GENERIC.
*/

#include <string.h>
#include "nmod_mat.h"

#if FLINT_BITS == 64

#include "nmod.h"
#include "machine_vectors.h"
#include "thread_pool.h"
#include "thread_support.h"

/* backend selection *********************************************************/

/*
    NMOD_MAT_U32_FORCE_AVX2 / _FORCE_NEON / _FORCE_GENERIC select a
    backend other than the best one available, for testing and
    profiling; FLINT_MACHINE_VECTORS_FORCE_GENERIC (honoured by
    machine_vectors.h itself) also forces the plain C path here.
*/
#if defined(FLINT_MACHINE_VECTORS_GENERIC) \
        || defined(NMOD_MAT_U32_FORCE_GENERIC)
# define U32_GENERIC 1
#elif defined(__AVX512F__) && defined(__AVX512DQ__) \
        && !defined(NMOD_MAT_U32_FORCE_AVX2) \
        && !defined(NMOD_MAT_U32_FORCE_NEON)
# define U32_AVX512 1
#elif defined(__AVX2__) && !defined(NMOD_MAT_U32_FORCE_NEON)
# define U32_AVX2 1
#elif (defined(__ARM_NEON) && defined(__aarch64__)) || defined(_M_ARM64)
/* AArch64 only: the canonical reduction converts int64 <-> double in
   vector registers, which ARMv7 NEON cannot do */
# define U32_NEON 1
#else
# define U32_GENERIC 1
#endif

/*
    Tile geometry. An accumulator vector holds U32_VL int64 lanes; a
    tile is U32_MR rows by U32_NACC accumulator vectors, that is
    U32_NR = U32_VL * U32_NACC columns. The MR*NACC accumulators plus
    the NACC B operands and one A splat must fit the register file.
*/
/* the x86 backends interleave two column groups per packed row, so
   their U32_NACC is fixed at 2 (see u32_load_bstep) */
#if defined(U32_AVX512)
# define U32_VL 8
# define U32_NACC 2
# ifndef U32_MR
#  define U32_MR 12     /* 24 of 32 zmm */
# endif
#elif defined(U32_AVX2)
# define U32_VL 4
# define U32_NACC 2
# ifndef U32_MR
#  define U32_MR 6      /* 12 of 16 ymm */
# endif
#elif defined(U32_NEON)
# define U32_VL 2
# ifndef U32_NACC
#  define U32_NACC 4
# endif
# ifndef U32_MR
/*
    vmlal_n_s32 wants the A entry in a vector lane, so a k step needs
    MR*NACC accumulators plus NACC B operands plus MR A operands live at
    once; 4 x 4 is the widest that fits 32 q registers without spilling
    (checked on the generated code) and balances A against B reuse.
*/
#  define U32_MR 4
# endif
#else
# define U32_VL 4
# ifndef U32_NACC
#  define U32_NACC 2
# endif
# ifndef U32_MR
#  define U32_MR 4
# endif
#endif

#define U32_NR (U32_VL * U32_NACC)

/* blocking parameters, in elements */
#ifndef U32_KC
# define U32_KC 256
#endif
#ifndef U32_MC
# define U32_MC 96
#endif
#ifndef U32_NC
# define U32_NC 2048
#endif

/*
    products (m*k*n) per thread: the split uses at most work/U32_MT_MIN_WORK
    threads, so with the value below a second thread joins from about 100^3
    and a fourth from about 126^3, which is where nmod_mat_mul_blas and
    nmod_mat_mul_classical_threaded start to gain from their threads too
    (measured on four machines with 4 threads; below 100^3 the serial
    kernel still beats them all)
*/
#ifndef U32_MT_MIN_WORK
# define U32_MT_MIN_WORK 500000.0
#endif

/* modulus-derived parameters ************************************************/

/*
    c32      2^32 mod n, which is <= most 2^31 - 1, so it is
             nonnegative as a signed 32-bit multiplier.
    rounds   number of fold rounds (1 or 2) applied at each in-kernel
             reduction and at the end of a block.
    fold_bound  bound on |acc| after the fold rounds.
    cadence  number of products that may be accumulated between folds:
             fold_bound + cadence*P < 2^63 where P = floor(n/2)^2.
    ninv     1.0/n, for the final Barrett step.
*/
typedef struct
{
    ulong n;
    ulong c32;
    ulong fold_bound;
    slong cadence;
    int rounds;
    double ninv;
}
u32_ctx_struct;

static void
u32_ctx_init(u32_ctx_struct * ctx, ulong n)
{
    ulong P, c32, R1, R2, init, half;

    FLINT_ASSERT(n >= 1 && n < (UWORD(1) << 32));

    ctx->n = n;
    ctx->ninv = 1.0 / (double) n;

    c32 = (UWORD(1) << 32) % n;
    ctx->c32 = c32;
    FLINT_ASSERT(c32 < (UWORD(1) << 31));

    half = n / 2;
    P = half * half;

    /*
        Fold round: acc = hi*2^32 + lo, |hi| <= 2^31, 0 <= lo < 2^32,
        acc' = hi*c32 + lo, so |acc'| <= 2^31*c32 + 2^32 =: R1, which is
        below 2^62 + 2^32. A second round starts from |hi| <= R1/2^32 + 1
        and gives |acc''| <= (floor(R1/2^32) + 1)*c32 + 2^32 =: R2.
    */
    R1 = (UWORD(1) << 31) * c32 + (UWORD(1) << 32);
    R2 = ((R1 >> 32) + 1) * c32 + (UWORD(1) << 32);

    /* one round is enough when it already leaves a negligible residue;
       this is the case for n close to 2^32 or to 2^31 (small c32),
       e.g. 2^31 - 1 or 2^32 - 5 */
    if (R1 <= (UWORD(1) << 40))
    {
        ctx->rounds = 1;
        ctx->fold_bound = R1;
    }
    else
    {
        ctx->rounds = 2;
        ctx->fold_bound = R2;
    }

    /* the accumulator starts either from a fold (|acc| <= fold_bound)
       or from a canonical entry of C (< n) */
    init = FLINT_MAX(ctx->fold_bound, n);

    if (P == 0)
    {
        ctx->cadence = U32_KC;
    }
    else
    {
        ulong cad = ((UWORD(1) << 63) - 1 - init) / P;
        ctx->cadence = (slong) FLINT_MIN(cad, (ulong) U32_KC);
    }

    /* worst case is n = 2^32 - 1 with P about 2^62, giving 2 */
    FLINT_ASSERT(ctx->cadence >= 1);
}

/* symmetric lift of a canonical residue to int32: a - n if a > n/2 */
FLINT_FORCE_INLINE int32_t
u32_lift(ulong a, ulong n)
{
    return (int32_t) (a - (n & FLINT_SIGN_EXT(n / 2 - a)));
}

/*
    Backend primitives. Each backend provides

      u32_acc          accumulator vector, U32_VL signed 64-bit lanes
      u32_bv           B operand feeding one accumulator, U32_VL values
      u32_av           A operand (one lifted entry, splat or scalar)
      u32_consts       loop-invariant constants derived from the modulus

      u32_consts_init  build them
      u32_acc_zero     all-zero accumulator
      u32_load_c       U32_VL canonical residues of C -> accumulator
      u32_store_c      accumulator (canonical residues) -> C
      u32_load_c32, u32_store_c32   the same for uint32 entries
      u32_bslot        where column j of a packed B row is stored
      u32_load_bstep   the U32_NACC B operands of one k step
      u32_load_a       one int32 of a packed A panel
      u32_mul_add      acc + a*b, widening 32x32 -> 64 per lane
      u32_fold1        one non-canonical fold round
      u32_reduce       canonical residue, for |v| <= fold_bound

    Only these differ between AVX-512, AVX2, NEON and plain C; the
    microkernel, the packing and the blocked core below are shared.
*/

#if defined(U32_AVX512)

/* AVX-512 (F + DQ) *********************************************************/

typedef __m512i u32_acc;
typedef __m512i u32_bv;
typedef __m512i u32_av;

typedef struct
{
    __m512i mask32;
    __m512i c32v;
    __m512i nv;
    __m512d ninvv;
}
u32_consts;

FLINT_FORCE_INLINE u32_consts
u32_consts_init(const u32_ctx_struct * ctx)
{
    u32_consts C;

    C.mask32 = _mm512_set1_epi64(UWORD(0xFFFFFFFF));
    C.c32v = _mm512_set1_epi64(ctx->c32);
    C.nv = _mm512_set1_epi64(ctx->n);
    C.ninvv = _mm512_set1_pd(ctx->ninv);

    return C;
}

FLINT_FORCE_INLINE u32_acc u32_acc_zero(void) { return _mm512_setzero_si512(); }

FLINT_FORCE_INLINE u32_acc
u32_load_c(const ulong * p) { return _mm512_loadu_si512((const void *) p); }

FLINT_FORCE_INLINE void
u32_store_c(ulong * p, u32_acc a) { _mm512_storeu_si512((void *) p, a); }

/* uint32 entries: widen on load, narrow on store (vpmovzxdq, vpmovqd) */
FLINT_FORCE_INLINE u32_acc
u32_load_c32(const uint32_t * p)
{
    return _mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i *) p));
}

FLINT_FORCE_INLINE void
u32_store_c32(uint32_t * p, u32_acc a)
{
    _mm256_storeu_si256((__m256i *) p, _mm512_cvtepi64_epi32(a));
}

/*
    x86 reads the even 32-bit halves of the 64-bit lanes, so a packed B
    row holds the two column groups interleaved (slot 2j is column j,
    slot 2j+1 is column j + VL) and one aligned load plus one shift
    yields both operands. Sign extending two separate loads instead
    would cost the same instructions but put them on the shuffle port,
    which the A broadcasts already saturate.
*/
FLINT_FORCE_INLINE slong
u32_bslot(slong j)
{
    return (j < U32_VL) ? 2 * j : 2 * (j - U32_VL) + 1;
}

FLINT_FORCE_INLINE void
u32_load_bstep(u32_bv * bv, const int32_t * p)
{
    __m512i t = _mm512_load_si512((const void *) p);

    bv[0] = t;
    bv[1] = _mm512_srli_epi64(t, 32);
}

FLINT_FORCE_INLINE u32_av
u32_load_a(const int32_t * p) { return _mm512_set1_epi32(*p); }

FLINT_FORCE_INLINE u32_acc
u32_mul_add(u32_acc acc, u32_av a, u32_bv b)
{
    return _mm512_add_epi64(acc, _mm512_mul_epi32(a, b));
}

FLINT_FORCE_INLINE u32_acc
u32_fold1(u32_acc acc, const u32_consts * C)
{
    __m512i hi = _mm512_srli_epi64(acc, 32);
    __m512i lo = _mm512_and_si512(acc, C->mask32);

    return _mm512_add_epi64(_mm512_mul_epi32(hi, C->c32v), lo);
}

/* canonical residue of v, for |v| <= fold_bound < 2^62 */
FLINT_FORCE_INLINE u32_acc
u32_reduce(u32_acc v, const u32_consts * C)
{
    __m512d d = _mm512_cvtepi64_pd(v);
    __m512d qd = _mm512_mul_pd(d, C->ninvv);
    __m512i q, r;
    __mmask8 m;

    qd = _mm512_roundscale_pd(qd, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    q = _mm512_cvttpd_epi64(qd);
    r = _mm512_sub_epi64(v, _mm512_mullo_epi64(q, C->nv));

    m = _mm512_cmplt_epi64_mask(r, _mm512_setzero_si512());
    r = _mm512_mask_add_epi64(r, m, r, C->nv);
    m = _mm512_cmpge_epi64_mask(r, C->nv);
    r = _mm512_mask_sub_epi64(r, m, r, C->nv);

    return r;
}

#elif defined(U32_AVX2)

/* AVX2 *********************************************************************/

typedef __m256i u32_acc;
typedef __m256i u32_bv;
typedef __m256i u32_av;

typedef struct
{
    __m256i mask32;
    __m256i c32v;
    __m256i nv;
    __m256i nm1v;
    __m256d ninvv;
    int nbig;           /* n >= 2^31: vpmuldq reads n as n - 2^32 */
}
u32_consts;

FLINT_FORCE_INLINE u32_consts
u32_consts_init(const u32_ctx_struct * ctx)
{
    u32_consts C;

    C.mask32 = _mm256_set1_epi64x(UWORD(0xFFFFFFFF));
    C.c32v = _mm256_set1_epi64x(ctx->c32);
    C.nv = _mm256_set1_epi64x(ctx->n);
    C.nm1v = _mm256_set1_epi64x(ctx->n - 1);
    C.ninvv = _mm256_set1_pd(ctx->ninv);
    C.nbig = (ctx->n >= (UWORD(1) << 31));

    return C;
}

FLINT_FORCE_INLINE u32_acc u32_acc_zero(void) { return _mm256_setzero_si256(); }

FLINT_FORCE_INLINE u32_acc
u32_load_c(const ulong * p) { return _mm256_loadu_si256((const __m256i *) p); }

FLINT_FORCE_INLINE void
u32_store_c(ulong * p, u32_acc a) { _mm256_storeu_si256((__m256i *) p, a); }

/* uint32 entries: widen on load; AVX2 has no narrowing store, so the
   low halves (the residues are below 2^32) are gathered by a permute */
FLINT_FORCE_INLINE u32_acc
u32_load_c32(const uint32_t * p)
{
    return _mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i *) p));
}

FLINT_FORCE_INLINE void
u32_store_c32(uint32_t * p, u32_acc a)
{
    const __m256i idx = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);

    _mm_storeu_si128((__m128i *) p,
            _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(a, idx)));
}

/* interleaved column groups; see the AVX-512 comment */
FLINT_FORCE_INLINE slong
u32_bslot(slong j)
{
    return (j < U32_VL) ? 2 * j : 2 * (j - U32_VL) + 1;
}

FLINT_FORCE_INLINE void
u32_load_bstep(u32_bv * bv, const int32_t * p)
{
    __m256i t = _mm256_load_si256((const __m256i *) p);

    bv[0] = t;
    bv[1] = _mm256_srli_epi64(t, 32);
}

FLINT_FORCE_INLINE u32_av
u32_load_a(const int32_t * p) { return _mm256_set1_epi32(*p); }

FLINT_FORCE_INLINE u32_acc
u32_mul_add(u32_acc acc, u32_av a, u32_bv b)
{
    return _mm256_add_epi64(acc, _mm256_mul_epi32(a, b));
}

FLINT_FORCE_INLINE u32_acc
u32_fold1(u32_acc acc, const u32_consts * C)
{
    __m256i hi = _mm256_srli_epi64(acc, 32);
    __m256i lo = _mm256_and_si256(acc, C->mask32);

    return _mm256_add_epi64(_mm256_mul_epi32(hi, C->c32v), lo);
}

/*
    Canonical residue of v, for |v| <= fold_bound < 2^62.

    AVX2 has no int64 <-> double conversion, so v = hi*2^32 + lo is
    converted in halves with the 2^52 exponent trick (hi is signed and
    is offset by 2^31 first). The quotient q = round(v/n) satisfies
    |q| < 2^31 (fold_bound/n is below 2^30 + 2^32/n), so it is exact in
    32 bits and q*n is computed by vpmuldq; that reads n as a signed
    32-bit value, which for n >= 2^31 is n - 2^32, corrected by adding
    q*2^32. The error in d is at most 2^9 absolute and the double bound
    is exact below 2^53, so |q - v/n| < 1/2 + 2^9/n, which keeps
    r = v - q*n in (-n, n); a second correction is included anyway, it
    is cheap here.
*/
FLINT_FORCE_INLINE u32_acc
u32_reduce(u32_acc v, const u32_consts * C)
{
    const __m256i magic = _mm256_set1_epi64x(UWORD(0x4330000000000000)); /* 2^52 */
    const __m256i sign32 = _mm256_set1_epi64x(UWORD(0x80000000));
    const __m256d magic_d = _mm256_castsi256_pd(magic);
    const __m256d magic_hi_d = _mm256_set1_pd(0x1.0p52 + 0x1.0p31);
    __m256i hi, lo, q64, r, msk;
    __m256d lo_d, hi_d, d, qd;
    __m128i q32;

    hi = _mm256_srli_epi64(v, 32);
    lo = _mm256_and_si256(v, C->mask32);
    lo_d = _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(lo, magic)), magic_d);
    hi_d = _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(
                            _mm256_xor_si256(hi, sign32), magic)), magic_hi_d);
    d = _mm256_add_pd(_mm256_mul_pd(hi_d, _mm256_set1_pd(0x1.0p32)), lo_d);

    qd = _mm256_mul_pd(d, C->ninvv);
    qd = _mm256_round_pd(qd, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    q32 = _mm256_cvttpd_epi32(qd);
    q64 = _mm256_cvtepi32_epi64(q32);

    r = _mm256_sub_epi64(v, _mm256_mul_epi32(q64, C->nv));
    if (C->nbig)
        r = _mm256_sub_epi64(r, _mm256_slli_epi64(q64, 32));

    msk = _mm256_cmpgt_epi64(_mm256_setzero_si256(), r);
    r = _mm256_add_epi64(r, _mm256_and_si256(msk, C->nv));
    msk = _mm256_cmpgt_epi64(r, C->nm1v);
    r = _mm256_sub_epi64(r, _mm256_and_si256(msk, C->nv));

    return r;
}

#elif defined(U32_NEON)

/* AArch64 NEON *************************************************************/

/*
    NEON widens from narrow operands rather than from the even halves of
    64-bit lanes: smlal (vmlal_n_s32) takes an int32x2_t and a scalar
    and accumulates into int64x2_t, so the multiply-accumulate is a
    single instruction and no splat of A into a vector is needed.
*/
typedef int64x2_t u32_acc;
typedef int32x2_t u32_bv;
typedef int32_t u32_av;

typedef struct
{
    int64x2_t mask32;
    int64x2_t nv;
    int64x2_t nm1v;
    int64x2_t zero;
    float64x2_t ninvv;
    int32_t c32;
    int32_t n32;
    int nbig;           /* n >= 2^31: smull reads n as n - 2^32 */
}
u32_consts;

FLINT_FORCE_INLINE u32_consts
u32_consts_init(const u32_ctx_struct * ctx)
{
    u32_consts C;

    C.mask32 = vdupq_n_s64((int64_t) UWORD(0xFFFFFFFF));
    C.nv = vdupq_n_s64((int64_t) ctx->n);
    C.nm1v = vdupq_n_s64((int64_t) (ctx->n - 1));
    C.zero = vdupq_n_s64(0);
    C.ninvv = vdupq_n_f64(ctx->ninv);
    C.c32 = (int32_t) ctx->c32;
    C.n32 = (int32_t) ctx->n;
    C.nbig = (ctx->n >= (UWORD(1) << 31));

    return C;
}

FLINT_FORCE_INLINE u32_acc u32_acc_zero(void) { return vdupq_n_s64(0); }

FLINT_FORCE_INLINE u32_acc
u32_load_c(const ulong * p)
{
    return vreinterpretq_s64_u64(vld1q_u64((const uint64_t *) p));
}

FLINT_FORCE_INLINE void
u32_store_c(ulong * p, u32_acc a)
{
    vst1q_u64((uint64_t *) p, vreinterpretq_u64_s64(a));
}

/* uint32 entries: widen on load, narrow on store */
FLINT_FORCE_INLINE u32_acc
u32_load_c32(const uint32_t * p)
{
    return vreinterpretq_s64_u64(vmovl_u32(vld1_u32(p)));
}

FLINT_FORCE_INLINE void
u32_store_c32(uint32_t * p, u32_acc a)
{
    vst1_u32(p, vmovn_u64(vreinterpretq_u64_s64(a)));
}

/* NEON takes narrow operands, so the columns stay in order */
FLINT_FORCE_INLINE slong
u32_bslot(slong j) { return j; }

FLINT_FORCE_INLINE void
u32_load_bstep(u32_bv * bv, const int32_t * p)
{
    slong v;

    for (v = 0; v < U32_NACC; v++)
        bv[v] = vld1_s32(p + v * U32_VL);
}

FLINT_FORCE_INLINE u32_av
u32_load_a(const int32_t * p) { return *p; }

FLINT_FORCE_INLINE u32_acc
u32_mul_add(u32_acc acc, u32_av a, u32_bv b) { return vmlal_n_s32(acc, b, a); }

FLINT_FORCE_INLINE u32_acc
u32_fold1(u32_acc acc, const u32_consts * C)
{
    int64x2_t hi = vshrq_n_s64(acc, 32);                /* arithmetic */
    int64x2_t lo = vandq_s64(acc, C->mask32);

    /* |hi| <= 2^31 fits in int32, so narrow and widen back through the
       32x32 -> 64 multiply */
    return vaddq_s64(vmull_n_s32(vmovn_s64(hi), C->c32), lo);
}

/* canonical residue of v, for |v| <= fold_bound < 2^62; see the AVX2
   version for why |q| < 2^31 and why n >= 2^31 needs a correction */
FLINT_FORCE_INLINE u32_acc
u32_reduce(u32_acc v, const u32_consts * C)
{
    float64x2_t d = vcvtq_f64_s64(v);
    float64x2_t qd = vrndnq_f64(vmulq_f64(d, C->ninvv));
    int64x2_t q = vcvtq_s64_f64(qd);
    int64x2_t r;

    r = vsubq_s64(v, vmull_n_s32(vmovn_s64(q), C->n32));
    if (C->nbig)
        r = vsubq_s64(r, vshlq_n_s64(q, 32));

    r = vaddq_s64(r, vandq_s64(
            vreinterpretq_s64_u64(vcltq_s64(r, C->zero)), C->nv));
    r = vsubq_s64(r, vandq_s64(
            vreinterpretq_s64_u64(vcgtq_s64(r, C->nm1v)), C->nv));

    return r;
}

#else

/* plain C ******************************************************************/

/*
    The same algorithm on small structs, which the compiler is free to
    vectorize; U32_VL is a logical width, not a hardware one.
*/
typedef struct { slong v[U32_VL]; } u32_acc;
typedef struct { int32_t v[U32_VL]; } u32_bv;
typedef slong u32_av;

typedef struct
{
    ulong c32;
    slong n;
    double ninv;
}
u32_consts;

FLINT_FORCE_INLINE u32_consts
u32_consts_init(const u32_ctx_struct * ctx)
{
    u32_consts C;

    C.c32 = ctx->c32;
    C.n = (slong) ctx->n;
    C.ninv = ctx->ninv;

    return C;
}

FLINT_FORCE_INLINE u32_acc
u32_acc_zero(void)
{
    u32_acc a;
    slong i;

    for (i = 0; i < U32_VL; i++)
        a.v[i] = 0;

    return a;
}

FLINT_FORCE_INLINE u32_acc
u32_load_c(const ulong * p)
{
    u32_acc a;
    slong i;

    for (i = 0; i < U32_VL; i++)
        a.v[i] = (slong) p[i];

    return a;
}

FLINT_FORCE_INLINE void
u32_store_c(ulong * p, u32_acc a)
{
    slong i;

    for (i = 0; i < U32_VL; i++)
        p[i] = (ulong) a.v[i];
}

FLINT_FORCE_INLINE u32_acc
u32_load_c32(const uint32_t * p)
{
    u32_acc a;
    slong i;

    for (i = 0; i < U32_VL; i++)
        a.v[i] = (slong) p[i];

    return a;
}

FLINT_FORCE_INLINE void
u32_store_c32(uint32_t * p, u32_acc a)
{
    slong i;

    for (i = 0; i < U32_VL; i++)
        p[i] = (uint32_t) a.v[i];
}

FLINT_FORCE_INLINE slong
u32_bslot(slong j) { return j; }

FLINT_FORCE_INLINE void
u32_load_bstep(u32_bv * bv, const int32_t * p)
{
    slong v, i;

    for (v = 0; v < U32_NACC; v++)
        for (i = 0; i < U32_VL; i++)
            bv[v].v[i] = p[v * U32_VL + i];
}

FLINT_FORCE_INLINE u32_av
u32_load_a(const int32_t * p) { return (slong) *p; }

FLINT_FORCE_INLINE u32_acc
u32_mul_add(u32_acc acc, u32_av a, u32_bv b)
{
    slong i;

    for (i = 0; i < U32_VL; i++)
        acc.v[i] += a * (slong) b.v[i];

    return acc;
}

FLINT_FORCE_INLINE u32_acc
u32_fold1(u32_acc acc, const u32_consts * C)
{
    slong i;

    for (i = 0; i < U32_VL; i++)
    {
        slong hi = acc.v[i] >> 32;                    /* arithmetic */
        ulong lo = (ulong) acc.v[i] & UWORD(0xFFFFFFFF);

        acc.v[i] = hi * (slong) C->c32 + (slong) lo;
    }

    return acc;
}

FLINT_FORCE_INLINE u32_acc
u32_reduce(u32_acc acc, const u32_consts * C)
{
    slong i;

    for (i = 0; i < U32_VL; i++)
    {
        slong v = acc.v[i];
        slong q = (slong) ((double) v * C->ninv);
        slong r = v - q * C->n;

        /* q is v/n truncated towards zero after an error below
           2^9/n + eps, so it is within 2 of floor(v/n) */
        if (r < 0) r += C->n;
        if (r < 0) r += C->n;
        if (r >= C->n) r -= C->n;
        if (r >= C->n) r -= C->n;

        acc.v[i] = r;
    }

    return acc;
}

#endif

/* shared across backends ****************************************************/

FLINT_FORCE_INLINE u32_acc
u32_fold(u32_acc acc, const u32_consts * C, int rounds)
{
    acc = u32_fold1(acc, C);
    if (rounds > 1)
        acc = u32_fold1(acc, C);

    return acc;
}

/*
    The packing, the microkernel, the blocked core and the thread split
    come from mul_blocked_templ.h, instantiated once for ulong entries
    (u32_64_*, behind nmod_mat_mul_u32) and once for uint32 entries
    (u32_32_*, behind _nmod_mat_mul_u32). Only the C accessors differ.
*/
#define BT_PACKED int32_t
#define BT_CTX u32_ctx_struct
#define BT_MR U32_MR
#define BT_VL U32_VL
#define BT_NACC U32_NACC
#define BT_KC U32_KC
#define BT_MC U32_MC
#define BT_NC U32_NC
#define BT_MT_MIN_WORK U32_MT_MIN_WORK
#define BT_LIFT(x, ctx) u32_lift((ulong) (x), (ctx)->n)
#define BT_BSLOT(j) u32_bslot(j)
#define BT_ACC u32_acc
#define BT_BV u32_bv
#define BT_AV u32_av
#define BT_CONSTS u32_consts
#define BT_CONSTS_INIT(ctx) u32_consts_init(ctx)
#define BT_ACC_ZERO() u32_acc_zero()
#define BT_LOAD_BSTEP(bv, p) u32_load_bstep(bv, p)
#define BT_LOAD_A(p) u32_load_a(p)
#define BT_MUL_ADD(acc, a, b, C) u32_mul_add(acc, a, b)
#define BT_CADENCE(ctx) ((ctx)->cadence)
#define BT_FOLD(acc, C, ctx) u32_fold(acc, C, (ctx)->rounds)
#define BT_FINISH(acc, C, ctx) u32_reduce(u32_fold(acc, C, (ctx)->rounds), C)

/* ulong entries */
#define BT_NAME(x) u32_64_##x
#define BT_ENTRY ulong
#define BT_LOAD_C(p) u32_load_c(p)
#define BT_STORE_C(p, acc) u32_store_c(p, acc)
#include "mul_blocked_templ.h"
#undef BT_NAME
#undef BT_ENTRY
#undef BT_LOAD_C
#undef BT_STORE_C

/* uint32 entries */
#define BT_NAME(x) u32_32_##x
#define BT_ENTRY uint32_t
#define BT_LOAD_C(p) u32_load_c32(p)
#define BT_STORE_C(p, acc) u32_store_c32(p, acc)
#include "mul_blocked_templ.h"
#undef BT_NAME
#undef BT_ENTRY
#undef BT_LOAD_C
#undef BT_STORE_C

/* public entries ************************************************************/

int
nmod_mat_mul_u32(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    ulong modn = C->mod.n;
    u32_ctx_struct ctx;

    FLINT_ASSERT(C->r == A->r);
    FLINT_ASSERT(C->c == B->c);
    FLINT_ASSERT(A->c == B->r);

    if (modn >= (UWORD(1) << 32))
        return 0;

    if (m <= 0 || n <= 0)
        return 1;

    if (k <= 0 || modn == 1)
    {
        nmod_mat_zero(C);
        return 1;
    }

    if (C == A || C == B)
    {
        nmod_mat_t T;
        nmod_mat_init(T, m, n, modn);
        nmod_mat_mul_u32(T, A, B);
        nmod_mat_swap_entrywise(C, T);
        nmod_mat_clear(T);
        return 1;
    }

    u32_ctx_init(&ctx, modn);

    u32_64_core_mt(C->entries, C->stride, A->entries, A->stride,
                   B->entries, B->stride, m, k, n, &ctx,
                   flint_get_num_threads());

    return 1;
}

int
_nmod_mat_mul_u32(uint32_t * C, slong Cstride,
                  const uint32_t * A, slong Astride,
                  const uint32_t * B, slong Bstride,
                  slong m, slong k, slong n, nmod_t mod)
{
    u32_ctx_struct ctx;
    slong i;

    if (mod.n >= (UWORD(1) << 32))
        return 0;

    if (m <= 0 || n <= 0)
        return 1;

    if (k <= 0 || mod.n == 1)
    {
        for (i = 0; i < m; i++)
            memset(C + i * Cstride, 0, n * sizeof(uint32_t));
        return 1;
    }

    if (C == A || C == B)
    {
        uint32_t * T = flint_malloc(m * n * sizeof(uint32_t));

        _nmod_mat_mul_u32(T, n, A, Astride, B, Bstride, m, k, n, mod);
        for (i = 0; i < m; i++)
            memcpy(C + i * Cstride, T + i * n, n * sizeof(uint32_t));

        flint_free(T);
        return 1;
    }

    u32_ctx_init(&ctx, mod.n);

    u32_32_core_mt(C, Cstride, A, Astride, B, Bstride, m, k, n, &ctx,
                   flint_get_num_threads());

    return 1;
}

#else

/* the packing and folds here assume 64-bit words */
int
nmod_mat_mul_u32(nmod_mat_t FLINT_UNUSED(C),
                 const nmod_mat_t FLINT_UNUSED(A),
                 const nmod_mat_t FLINT_UNUSED(B))
{
    return 0;
}

int
_nmod_mat_mul_u32(uint32_t * FLINT_UNUSED(C), slong FLINT_UNUSED(Cstride),
                  const uint32_t * FLINT_UNUSED(A), slong FLINT_UNUSED(Astride),
                  const uint32_t * FLINT_UNUSED(B), slong FLINT_UNUSED(Bstride),
                  slong FLINT_UNUSED(m), slong FLINT_UNUSED(k),
                  slong FLINT_UNUSED(n), nmod_t FLINT_UNUSED(mod))
{
    return 0;
}

#endif
