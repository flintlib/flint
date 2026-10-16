/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    flint_sgemm and flint_dgemm: always-available single/double precision
    matrix multiplication, so that FLINT modules can use a GEMM building
    block whether or not an external BLAS is configured.

    C = A * B, row-major, no transposes, no accumulation (alpha = 1,
    beta = 0), arbitrary dimensions and leading dimensions. This is the
    subset of xGEMM that FLINT actually uses in nmod_mat and fmpz_mat.

    Structure: BLIS-style three-level blocking (NC/KC/MC), packed A
    blocks and B panels, a register-tile microkernel, plus streaming
    no-pack paths for thin and small shapes. All SIMD goes through
    machine_vectors.h. Multithreading uses FLINT's thread pool: one
    thread_pool_wake/wait per worker per call, with the phases inside
    separated by spin barriers.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
/*
    The generic vector backends express a fused multiply-add as the
    expression a * b + c, which the compiler contracts into an FMA only
    if it is allowed to. GCC in a strict ISO mode (FLINT builds with
    -std=c11) defaults to -ffp-contract=on and does not contract, which
    costs about 1.4x here; clang contracts under its default. Requesting
    fast contraction for this file restores it. Contraction is safe for
    everything in this file: no operation depends on the intermediate
    product being rounded, and an FMA is at least as accurate.

    The AVX2/AVX512/NEON backends are unaffected either way, since they
    use fused intrinsics directly.
*/
#if defined(__GNUC__) && !defined(__clang__)
# pragma GCC optimize("fp-contract=fast")
#endif

#include "flint.h"

/*
    For profiling and testing the portable code paths without rebuilding
    FLINT: uncomment one of these to override the vector backend that
    machine_vectors.h would select for this file. FLINT_MACHINE_VECTORS_-
    FORCE_GENERIC gives the GNU vector extension backend (the fallback on
    targets with no AVX2/NEON support in machine_vectors.h) and
    FLINT_MACHINE_VECTORS_STRICT_C gives the strict ISO C backend. Only
    this translation unit is affected, so the rest of FLINT keeps using
    the native vectors.

#define FLINT_MACHINE_VECTORS_FORCE_GENERIC
#define FLINT_MACHINE_VECTORS_STRICT_C
*/

#include "machine_vectors.h"
#include "thread_pool.h"
#include "thread_support.h"

/*
    The single point in FLINT where a BLAS header is included. Everything
    else calls flint_sgemm/flint_dgemm (or the _blas/_fallback variants)
    unconditionally.
*/
#if FLINT_USES_BLAS
# include <cblas.h>
#endif

/* Blocking parameters, in elements. */
#define SG_KC 256
#define SG_MC 96
#define SG_NC 2048

/* Widest vector per element type, and microkernel rows; MR x 2 vectors
   of accumulators plus the two B vectors must fit the register file. */
#if defined(FLINT_MACHINE_VECTORS_GENERIC)
/*
    Generic backends. The logical vector width follows the target: the
    GNU vector extension types lower directly onto whatever registers
    exist, and the ISO C struct loops are what the auto-vectorizer is
    given to work with, so a width narrower than the hardware wastes the
    datapath. 32 bytes is the largest used here, for two reasons: only
    the flat types (vec2d/vec4d, vec4f/vec8f) are single objects, the
    wider ones being structs of two, which doubles register pressure and
    defeats auto-vectorization of the ISO C tier; and MR x 2 accumulators
    of 32 bytes plus two B vectors is exactly the 16 architectural
    registers of pre-AVX512 x86.
*/
# if defined(FLINT_MACHINE_VECTORS_STRICT_C)
#  define GEMM_VD vec4d
#  define GEMM_VF GEMM_SC_VF
#  define GEMM_VLD 4
#  define GEMM_VLF GEMM_SC_VLF
# elif defined(__AVX__)
#  define GEMM_VD vec4d
#  define GEMM_VF vec8f
#  define GEMM_VLD 4
#  define GEMM_VLF 8
# else
#  define GEMM_VD vec2d
#  define GEMM_VF vec4f
#  define GEMM_VLD 2
#  define GEMM_VLF 4
# endif
# ifndef GEMM_MRD
#  define GEMM_MRD GEMM_GEN_MR
# endif
# ifndef GEMM_MRF
#  define GEMM_MRF GEMM_GEN_MR
# endif
#elif defined(__AVX512F__)
# define GEMM_VD vec8dz
# define GEMM_VF vec16fz
# define GEMM_VLD 8
# define GEMM_VLF 16
# define GEMM_MRD 8
# define GEMM_MRF 8
#elif defined(__AVX2__)
# define GEMM_VD vec4d
# define GEMM_VF vec8f
# define GEMM_VLD 4
# define GEMM_VLF 8
# define GEMM_MRD 6
# define GEMM_MRF 6
#elif defined(__ARM_NEON) || defined(_M_ARM64)
# define GEMM_VD vec2d
# define GEMM_VF vec4f
# define GEMM_VLD 2
# define GEMM_VLF 4
# define GEMM_MRD 8
# define GEMM_MRF 8
#else
/* Generic backends: 128-bit logical vectors. Every target reaching this
   branch (neither AVX2 nor NEON) has at most 128-bit SIMD and 16 vector
   registers in practice, and wider logical vectors then spill. */
# define GEMM_VD vec2d
# define GEMM_VF vec4f
# define GEMM_VLD 2
# define GEMM_VLF 4
# define GEMM_MRD 6
# define GEMM_MRF 6
#endif

/* The direct no-pack kernel for small problems is a win when the vector
   layer is real SIMD, and a loss for the strict ISO C structs, where the
   compiler vectorizes the k-loop over a strided B instead of the lanes. */
#ifndef GEMM_SMALL_DIM
# if defined(FLINT_MACHINE_VECTORS_STRICT_C)
#  define GEMM_SMALL_DIM 0
# else
#  define GEMM_SMALL_DIM 96
# endif
#endif

/* work below which threading cannot pay for itself */
#define SG_MT_FLOP_PER_THREAD 1.0e6

/*
    Two parallel drivers are available.

    GEMM_MT_ATOMIC: workers cooperate on one packed B panel at a time,
    coordinated by atomic work queues and sense-reversing spin barriers.
    This is the faster of the two, but it exercises C11 atomics and
    relaxed/acquire-release ordering that we have so far only stress
    tested on x86-64.

    Otherwise: the row range of C is split into one independent serial
    multiplication per worker. No atomics, no barriers, nothing shared
    but read-only inputs and disjoint output blocks. Slightly slower at
    high thread counts (B is packed once per worker rather than once per
    team) but portable with no assumptions beyond the thread pool.

    Define GEMM_MT_FORCE_ATOMIC or GEMM_MT_FORCE_SPLIT to override.
*/
#if defined(GEMM_MT_FORCE_ATOMIC)
# define GEMM_MT_ATOMIC 1
#elif defined(GEMM_MT_FORCE_SPLIT)
# undef GEMM_MT_ATOMIC
#elif defined(__x86_64__) || defined(_M_X64)
# define GEMM_MT_ATOMIC 1
#endif

#if defined(GEMM_MT_ATOMIC)
# include <stdatomic.h>
# if !defined(_WIN32) && !defined(_MSC_VER)
#  include <sched.h>
# endif
#endif

#ifndef GEMM_GEN_MR
# define GEMM_GEN_MR 4
#endif
#ifndef GEMM_SC_VF
# define GEMM_SC_VF vec4f
# define GEMM_SC_VLF 4
#endif

#define GEMM_CAT_(a, b) a##b
#define GEMM_CAT(a, b) GEMM_CAT_(a, b)

/*
    Packing scratch. With thread-local storage this is a grow-only buffer
    reused across calls, which matters for small matrices where an
    allocation would otherwise dominate; without TLS the serial core must
    allocate per call, since it may be entered concurrently.
*/
#if FLINT_USES_TLS

static FLINT_TLS_PREFIX char * gemm_scratch_buf = NULL;
static FLINT_TLS_PREFIX slong gemm_scratch_sz = 0;

static void
gemm_scratch_cleanup(void)
{
    flint_aligned_free(gemm_scratch_buf);
    gemm_scratch_buf = NULL;
    gemm_scratch_sz = 0;
}

#if !defined(GEMM_MT_ATOMIC)
/* release the calling thread's scratch: thread pool threads outlive the
   call and do not run registered cleanup functions, so the split driver
   makes each worker drop its cache before returning */
static void
gemm_scratch_free_local(void)
{
    flint_aligned_free(gemm_scratch_buf);
    gemm_scratch_buf = NULL;
    gemm_scratch_sz = 0;
}
#endif

static char *
gemm_get_scratch(slong bytes)
{
    if (bytes > gemm_scratch_sz)
    {
        slong sz = (bytes + 4095) & ~(slong) 4095;

        if (gemm_scratch_buf == NULL)
            flint_register_cleanup_function(gemm_scratch_cleanup);
        else
            flint_aligned_free(gemm_scratch_buf);

        gemm_scratch_buf = flint_aligned_alloc(64, sz);
        gemm_scratch_sz = sz;
    }

    return gemm_scratch_buf;
}

# define gemm_put_scratch(ptr) do { } while (0)

#else

static char *
gemm_get_scratch(slong bytes)
{
    return flint_aligned_alloc(64, (bytes + 63) & ~(slong) 63);
}

# define gemm_put_scratch(ptr) flint_aligned_free(ptr)
# define gemm_scratch_free_local() do { } while (0)

#endif

#if defined(GEMM_MT_ATOMIC)

/* pause/yield hints for the spin barrier */
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64)
# if defined(__GNUC__)
#  define GEMM_CPU_RELAX() __builtin_ia32_pause()
# else
#  include <intrin.h>
#  define GEMM_CPU_RELAX() _mm_pause()
# endif
#elif defined(__aarch64__) || defined(__arm__)
# define GEMM_CPU_RELAX() __asm__ __volatile__("yield")
#else
# define GEMM_CPU_RELAX() ((void) 0)
#endif

static void
gemm_thread_yield(void)
{
#if defined(_WIN32) || defined(_MSC_VER)
    /* no portable yield; spinning is acceptable here since the barrier
       waits are short and oversubscription is the only case that hits
       this path */
#else
    sched_yield();
#endif
}

/* sense-reversing barrier shared by both stamps */
static void
gemm_mt_barrier(_Atomic int * count, _Atomic int * bsense,
                _Atomic slong * q1, _Atomic slong * q2,
                int nt, int * sense)
{
    *sense ^= 1;
    if (atomic_fetch_add_explicit(count, 1, memory_order_acq_rel) == nt - 1)
    {
        atomic_store_explicit(count, 0, memory_order_relaxed);
        atomic_store_explicit(q1, 0, memory_order_relaxed);
        atomic_store_explicit(q2, 0, memory_order_relaxed);
        atomic_store_explicit(bsense, *sense, memory_order_release);
    }
    else
    {
        unsigned spins = 0;

        while (atomic_load_explicit(bsense, memory_order_acquire) != *sense)
        {
            if (++spins < 2048)
                GEMM_CPU_RELAX();
            else
            {
                spins = 0;
                gemm_thread_yield();
            }
        }
    }
}

#endif /* GEMM_MT_ATOMIC */

/* single precision **********************************************************/

#define GEMM_PART_SERIAL
#define GEMM_PART_MT
#define GT_PRE flint_gemm_f32_
#define GT_T float
#define GT_VEC GEMM_VF
#define GT_VL GEMM_VLF
#define GT_MR GEMM_MRF
#define GT_LOADU(p) GEMM_CAT(GEMM_VF, _load_unaligned)(p)
#define GT_LOADA(p) GEMM_CAT(GEMM_VF, _load_aligned)(p)
#define GT_STOREU(p, v) GEMM_CAT(GEMM_VF, _store_unaligned)(p, v)
#define GT_SETZERO() GEMM_CAT(GEMM_VF, _zero)()
#define GT_BCAST(p) GEMM_CAT(GEMM_VF, _set_f)(*(p))
#define GT_FMADD(a, b, c) GEMM_CAT(GEMM_VF, _fmadd)(a, b, c)
#include "gemm_templ.h"
#undef GT_PRE
#undef GT_T
#undef GT_VEC
#undef GT_VL
#undef GT_MR
#undef GT_LOADU
#undef GT_LOADA
#undef GT_STOREU
#undef GT_SETZERO
#undef GT_BCAST
#undef GT_FMADD

/* double precision **********************************************************/

#define GT_PRE flint_gemm_f64_
#define GT_T double
#define GT_VEC GEMM_VD
#define GT_VL GEMM_VLD
#define GT_MR GEMM_MRD
#define GT_LOADU(p) GEMM_CAT(GEMM_VD, _load_unaligned)(p)
#define GT_LOADA(p) GEMM_CAT(GEMM_VD, _load_aligned)(p)
#define GT_STOREU(p, v) GEMM_CAT(GEMM_VD, _store_unaligned)(p, v)
#define GT_SETZERO() GEMM_CAT(GEMM_VD, _zero)()
#define GT_BCAST(p) GEMM_CAT(GEMM_VD, _set_d)(*(p))
#define GT_FMADD(a, b, c) GEMM_CAT(GEMM_VD, _fmadd)(a, b, c)
#include "gemm_templ.h"
#undef GT_PRE
#undef GT_T
#undef GT_VEC
#undef GT_VL
#undef GT_MR
#undef GT_LOADU
#undef GT_LOADA
#undef GT_STOREU
#undef GT_SETZERO
#undef GT_BCAST
#undef GT_FMADD
#undef GEMM_PART_SERIAL
#undef GEMM_PART_MT

/* fallback implementations *************************************************/

void
flint_sgemm_fallback(slong m, slong k, slong n,
                     const float * A, slong lda,
                     const float * B, slong ldb,
                     float * C, slong ldc)
{
    slong i;

    if (m <= 0 || n <= 0)
        return;

    if (k <= 0)
    {
        for (i = 0; i < m; i++)
            memset(C + i * ldc, 0, n * sizeof(float));
        return;
    }

#if defined(GEMM_MT_ATOMIC)
    flint_gemm_f32_core_mt(C, ldc, A, lda, B, ldb, m, k, n,
                           flint_get_num_threads());
#else
    flint_gemm_f32_core_split_mt(C, ldc, A, lda, B, ldb, m, k, n,
                                 flint_get_num_threads());
#endif
}

void
flint_dgemm_fallback(slong m, slong k, slong n,
                     const double * A, slong lda,
                     const double * B, slong ldb,
                     double * C, slong ldc)
{
    slong i;

    if (m <= 0 || n <= 0)
        return;

    if (k <= 0)
    {
        for (i = 0; i < m; i++)
            memset(C + i * ldc, 0, n * sizeof(double));
        return;
    }

#if defined(GEMM_MT_ATOMIC)
    flint_gemm_f64_core_mt(C, ldc, A, lda, B, ldb, m, k, n,
                           flint_get_num_threads());
#else
    flint_gemm_f64_core_split_mt(C, ldc, A, lda, B, ldb, m, k, n,
                                 flint_get_num_threads());
#endif
}

/* BLAS wrappers ************************************************************/

void
flint_sgemm_blas(slong m, slong k, slong n,
                 const float * A, slong lda,
                 const float * B, slong ldb,
                 float * C, slong ldc)
{
#if FLINT_USES_BLAS
    if (m <= 0 || n <= 0)
        return;

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0f, A, lda, B, ldb, 0.0f, C, ldc);
#else
    (void) m; (void) k; (void) n;
    (void) A; (void) lda; (void) B; (void) ldb; (void) C; (void) ldc;

    flint_throw(FLINT_ERROR, "flint_sgemm_blas: FLINT was built without BLAS\n");
#endif
}

void
flint_dgemm_blas(slong m, slong k, slong n,
                 const double * A, slong lda,
                 const double * B, slong ldb,
                 double * C, slong ldc)
{
#if FLINT_USES_BLAS
    if (m <= 0 || n <= 0)
        return;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0, A, lda, B, ldb, 0.0, C, ldc);
#else
    (void) m; (void) k; (void) n;
    (void) A; (void) lda; (void) B; (void) ldb; (void) C; (void) ldc;

    flint_throw(FLINT_ERROR, "flint_dgemm_blas: FLINT was built without BLAS\n");
#endif
}

/* dispatch *****************************************************************/

/* FLINT_USES_BLAS is left undefined rather than 0 when BLAS is absent,
   so it cannot be used directly as an initializer */
#if FLINT_USES_BLAS
int flint_gemm_use_blas = 1;
#else
int flint_gemm_use_blas = 0;
#endif

void
flint_sgemm(slong m, slong k, slong n,
            const float * A, slong lda,
            const float * B, slong ldb,
            float * C, slong ldc)
{
    if (flint_gemm_use_blas)
        flint_sgemm_blas(m, k, n, A, lda, B, ldb, C, ldc);
    else
        flint_sgemm_fallback(m, k, n, A, lda, B, ldb, C, ldc);
}

void
flint_dgemm(slong m, slong k, slong n,
            const double * A, slong lda,
            const double * B, slong ldb,
            double * C, slong ldc)
{
    if (flint_gemm_use_blas)
        flint_dgemm_blas(m, k, n, A, lda, B, ldb, C, ldc);
    else
        flint_dgemm_fallback(m, k, n, A, lda, B, ldb, C, ldc);
}
