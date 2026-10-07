/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Template for the dN / dNb arithmetic, instantiated by d1.c ... d4.c
   with N, DF(name) -> dN_name and DB(name) -> dNb_name defined.

   Every ball operation follows one pattern: the midpoint operation is
   carried out with error-free transformations (TwoSum, TwoProd) that
   produce a list of exact terms, the list is renormalized to N
   components, and the magnitude of everything dropped is added to a
   running error bound.  The only a-priori bounds are for partial
   products of order >= N in multiplication, which are bounded by their
   own (rounded) magnitude.  This makes the radius zero exactly when
   the inputs are exact and no rounding occurred, and it makes rigor
   independent of any assumption on the inputs apart from the
   no-underflow precondition of TwoProd (checked; violations take the
   scaled slow path in slow.c). */

#define TINY_THRESHOLD 0x1p-480   /* components below this (and nonzero)
                                     could make a product underflow */

/* the kernels are written as loops with constant bounds; they must be
   fully unrolled to keep the terms in registers */
#if defined(__clang__)
/* (full unrolling is an optimization; clang reports the loops where it
   gives up, which is harmless) */
#pragma clang diagnostic ignored "-Wpass-failed"
#define UNROLL _Pragma("clang loop unroll(full)")
#elif defined(__GNUC__)
#define UNROLL _Pragma("GCC unroll 64")
#else
#define UNROLL
#endif
#if defined(__GNUC__) && !defined(__clang__)
#define DFLOAT_NO_UNROLL _Pragma("GCC unroll 1")
#else
#define DFLOAT_NO_UNROLL
#endif

/* a ulong below 2^52 converts to double exactly (every ulong does on
   a 32-bit target, whose doubles hold 53 bits) */
#define DFLOAT_ULONG_IS_EXACT(u) (FLINT_BITS == 32 || (u) < (UWORD(1) << (FLINT_BITS - 12)))

#if defined(__GNUC__)
#define DFLOAT_PREFETCH(p) __builtin_prefetch(p)
#else
#define DFLOAT_PREFETCH(p) ((void) 0)
#endif

/* Code-size control: functions that are large and either expensive
   (so that a call costs nothing in relative terms) or cold, and that
   would otherwise be inlined at several call sites, are kept out of
   line.  The ball operations take their error mode (automatic or
   fast) as a run-time argument: the automatic mode's body contains
   both the tracking and the bounding kernels, so the fast mode reuses
   it at the cost of a predictable branch. */
#if defined(__GNUC__) && !defined(__clang__)
#define DFLOAT_NOINLINE __attribute__((noinline, noclone))
#elif defined(__GNUC__)
#define DFLOAT_NOINLINE __attribute__((noinline))
#else
#define DFLOAT_NOINLINE
#endif

/* the dot products use the SIMD loop from this length on (for the
   plain d1 the scalar loop is cheap enough that the SIMD setup and the
   final combination of the lanes only pay off later) */
#define DFLOAT_DOT_SIMD_MIN 8
#ifndef DFLOAT_DOT1_SIMD_MIN
#define DFLOAT_DOT1_SIMD_MIN 32
#endif
#ifndef DFLOAT_DOT23_SIMD_MIN
#define DFLOAT_DOT23_SIMD_MIN 10
#endif
#define DFLOAT_DOT_PLAIN_SIMD_MIN ((N == 1) ? DFLOAT_DOT1_SIMD_MIN : ((N <= 3) ? DFLOAT_DOT23_SIMD_MIN : DFLOAT_DOT_SIMD_MIN))

/* the public elementary functions: never inlined into their callers
   (the gr wrappers and the compositions), where a call costs nothing
   in relative terms but an inlined copy costs 0.5-1.5 KB each */
#if defined(__GNUC__)
#define DFLOAT_FUNC_NOINLINE __attribute__((noinline))
#else
#define DFLOAT_FUNC_NOINLINE
#endif

/* The small operations are written as kernels, with the public
   functions as wrappers, and the code in this file (the gr wrappers,
   the vector loops, the dot products, the complex arithmetic) calls the
   kernels: in a shared library, a call to a public function from the
   same file cannot be inlined, or even made directly, since the symbol
   might be interposed; it goes through the PLT.  The trivial kernels
   DF(_name_inline) are always inlined (FLINT_FORCE_INLINE), as are the
   division and square roots for N = 1.  The arithmetic DF(_name_impl)
   (DFLOAT_ARITH_IMPL) is inlined for N = 1, where it is a few
   instructions; for N >= 2, where its bodies are larger (0.3-4 KB) and
   their cost dominates a call, it is a single static copy that is
   called directly. */
#if N == 1
#define DFLOAT_ARITH_IMPL FLINT_FORCE_INLINE
#else
#define DFLOAT_ARITH_IMPL static DFLOAT_NOINLINE
#endif

/* ------------------------------------------------------------------ */
/* Core kernels                                                        */
/* ------------------------------------------------------------------ */

/* Slow path of the renormalization: VecSumErrBranch (with TwoSum in
   place of Fast2Sum) over the residuals t[1..M-1] of a VecSum: terms
   are absorbed into the running component until a nonzero residual
   appears, which starts the next component.  Emitted components are
   nonzero, so the head is zero only when the whole sum is zero. */
static void
DF(_renorm_branchy)(double * res, double * err, double * t, int M, int track)
{
    int i, j;
    double f, eps, eps_new, e;

    eps = t[0];
    j = 0;

    for (i = 1; i < M; i++)
    {
        DFLOAT_TWO_SUM(f, eps_new, eps, t[i]);
        if (eps_new != 0.0)
        {
            res[j] = f;
            j++;
            if (j == N)
            {
                /* fold everything that remains into the last
                   component, dropping only the residuals (the terms
                   need not be sorted, so a late term may be large) */
                e = 0.0;
                f = res[N - 1];
                DFLOAT_TWO_SUM(f, eps, f, eps_new);
                e += fabs(eps);
                for (i = i + 1; i < M; i++)
                {
                    DFLOAT_TWO_SUM(f, eps, f, t[i]);
                    e += fabs(eps);
                }
                res[N - 1] = f;
                if (track)
                    *err += e;
                return;
            }
            eps = eps_new;
        }
        else
            eps = f;
    }

    res[j] = eps;
    for (j = j + 1; j < N; j++)
        res[j] = 0.0;
}

/* Renormalize the M terms t[0..M-1] (roughly sorted by decreasing
   magnitude) to N nonoverlapping components res[0..N-1]; the sum is
   preserved exactly except for the dropped terms, whose magnitudes are
   added to *err when track is set.  t is destroyed.

   Bottom-up VecSum makes t[0] the (nearly) correctly rounded sum of
   everything, with t[1..M-1] holding residuals.  The top-down pass
   condenses the residuals with a straight chain of N TwoSums (the
   first could be a Fast2Sum, but TwoSum keeps exactness independent
   of any ordering).  When a residual in that chain is zero the chain
   wastes a slot; if terms were dropped in that case, the branchy
   VecSumErrBranch pass is run instead, which packs the terms tightly
   (needed for exact results that use all 53 N bits). */
FLINT_FORCE_INLINE void
DF(_renorm_terms)(double * res, double * err, double * t, int M, int track)
{
    int i, bad;
    double g, e, r[N], sfx[N + 1];

    UNROLL
    for (i = M - 2; i >= 0; i--)
        DFLOAT_TWO_SUM(t[i], t[i + 1], t[i], t[i + 1]);

    if (M == 1)
    {
        res[0] = t[0];
        for (i = 1; i < N; i++)
            res[i] = 0.0;
        return;
    }

    /* suffix sums of the magnitudes: sfx[i] = |t[i+2]| + ... + |t[M-1]| */
    sfx[N] = 0.0;
    UNROLL
    for (i = M - 1; i >= N + 2; i--)
        sfx[N] += fabs(t[i]);
    UNROLL
    for (i = N - 1; i >= 0; i--)
        sfx[i] = sfx[i + 1] + ((i + 2 < M) ? fabs(t[i + 2]) : 0.0);

    g = t[1];
    bad = 0;
    UNROLL
    for (i = 0; i < N; i++)
    {
        if (i == 0)
        {
            DFLOAT_TWO_SUM(r[0], g, t[0], t[1]);
        }
        else if (i + 1 < M)
        {
            DFLOAT_TWO_SUM(r[i], g, g, t[i + 1]);
        }
        else
        {
            r[i] = g;
            g = 0.0;
        }
        /* a zero residual followed by nonzero terms leaves a gap */
        bad |= (g == 0.0) & (sfx[i] != 0.0);
    }

    if (FLINT_UNLIKELY(bad || r[0] == 0.0))
    {
        /* t[1..M-1] are still the VecSum residuals: repack */
        DF(_renorm_branchy)(res, err, t, M, track);
        return;
    }

    /* the remaining terms are negligible when the list was sorted by
       magnitude; otherwise fold them into the last component, dropping
       only the residuals */
    e = fabs(g);
    if (FLINT_UNLIKELY(sfx[N - 1] > 0x1p-40 * fabs(r[N - 1])))
    {
        UNROLL
        for (i = N + 1; i < M; i++)
        {
            DFLOAT_TWO_SUM(r[N - 1], g, r[N - 1], t[i]);
            e += fabs(g);
        }
    }
    else
        e += sfx[N - 1];

    UNROLL
    for (i = 0; i < N; i++)
        res[i] = r[i];

    if (track)
        *err += e;
}

/* (the error modes DFLOAT_MODE_AUTO and DFLOAT_MODE_FAST of the ball
   products are defined in dfloat.h; see DB(_mul_m)) */
#define TRACK_MODE(mode, exact) ((mode) == DFLOAT_MODE_FAST ? 2 : ((exact) ? 1 : 2))

/* Experiment knob: build the kernels with the strong renormalization
   pass (see _renorm_orders in kernels.inc) instead of applying it
   afterwards in the strong contexts. */
#ifndef DFLOAT_KERNEL_STRONG
#define DFLOAT_KERNEL_STRONG 0
#endif

/* The arithmetic kernels (kernels.inc), instantiated for scalars at
   N components (prefix _s) and at every smaller number of components
   M (prefix _sM), which the tapered algorithms use inline. */
#define K_ADD(a, b) ((a) + (b))
#define K_SUB(a, b) ((a) - (b))
#define K_MUL(a, b) ((a) * (b))
#define K_FMS(a, b, c) fma((a), (b), -(c))
#define K_NEG(a) (-(a))
#define K_ABS(a) fabs(a)
#define K_ZERO 0.0
#define K_TWICE(a) ((a) + (a))
#define K_SET_D(d) (d)
#define K_BOUND_BUMP_IF(a) (((a) != 0.0) ? K_BOUND_BUMP : 0.0)
#define K_ANY_ZERO(a) ((a) == 0.0)
#define K_ANY_GAP(a, b) ((a) == 0.0 && (b) != 0.0)
#define KV double
#define KM N
#define KN(name) DF(_s##name)
#include "kernels.inc"
#undef KM
#undef KN
#if N > 1
#define KM 1
#define KN(name) DF(_s1##name)
#include "kernels.inc"
#undef KM
#undef KN
#endif
#if N > 2
#define KM 2
#define KN(name) DF(_s2##name)
#include "kernels.inc"
#undef KM
#undef KN
#endif
#if N > 3
#define KM 3
#define KN(name) DF(_s3##name)
#include "kernels.inc"
#undef KM
#undef KN
#endif
#undef KV
#undef K_ADD
#undef K_SUB
#undef K_MUL
#undef K_FMS
#undef K_NEG
#undef K_ABS
#undef K_ZERO
#undef K_TWICE
#undef K_SET_D
#undef K_BOUND_BUMP_IF
#undef K_ANY_ZERO
#undef K_ANY_GAP

/* the scalar kernels at level M <= N (M constant after unrolling),
   with the zero-head repair; the higher components of the output are
   zeroed */
/* the names of the level kernels (those at or above N alias the
   N-level ones, so that the switches below compile; they are never
   reached) */
#if N > 1
#define DFS1(fn) DF(_s1##fn)
#define DFV1(fn) DF(_v1##fn)
#else
#define DFS1(fn) DF(_s##fn)
#define DFV1(fn) DF(_v##fn)
#endif
#if N > 2
#define DFS2(fn) DF(_s2##fn)
#define DFV2(fn) DF(_v2##fn)
#else
#define DFS2(fn) DF(_s##fn)
#define DFV2(fn) DF(_v##fn)
#endif
#if N > 3
#define DFS3(fn) DF(_s3##fn)
#define DFV3(fn) DF(_v3##fn)
#else
#define DFS3(fn) DF(_s##fn)
#define DFV3(fn) DF(_v##fn)
#endif

#define SCALAR_LEVEL_OP(flag, fn, M, ...) \
    do { \
        switch (M) \
        { \
            case 1: (flag) = DFS1(fn)(__VA_ARGS__); break; \
            case 2: (flag) = DFS2(fn)(__VA_ARGS__); break; \
            case 3: (flag) = DFS3(fn)(__VA_ARGS__); break; \
            default: (flag) = DF(_s##fn)(__VA_ARGS__); break; \
        } \
    } while (0)

FLINT_FORCE_INLINE void
DF(_shift_zero_head_m)(double * res, int M)
{
    int i, j = 0;
    for (i = 0; i < M; i++)
        if (res[i] != 0.0)
            res[j++] = res[i];
    for (; j < N; j++)
        res[j] = 0.0;
}

FLINT_FORCE_INLINE void
DF(_level_finish)(double * res, int flag, int M)
{
    int i;
    if (FLINT_UNLIKELY(flag))
        DF(_shift_zero_head_m)(res, M);
    for (i = 0; i < N; i++)
        if (i >= M)
            res[i] = 0.0;
}

FLINT_FORCE_INLINE void
DF(_level_mul_d)(double * res, const double * x, double y, int M)
{
    int flag;
    SCALAR_LEVEL_OP(flag, _mul_d_kernel, M, res, NULL, x, y, 0, 0);
    DF(_level_finish)(res, flag, M);
}

FLINT_FORCE_INLINE void
DF(_level_addsub)(double * res, const double * x, const double * y, int negate, int M)
{
    int flag;
    SCALAR_LEVEL_OP(flag, _addsub_kernel, M, res, NULL, x, y, negate, 0, 0);
    DF(_level_finish)(res, flag, M);
}

FLINT_FORCE_INLINE void
DF(_level_mul)(double * res, const double * x, const double * y, int M)
{
    int flag;
    SCALAR_LEVEL_OP(flag, _mul_kernel, M, res, NULL, x, y, 0, 0);
    DF(_level_finish)(res, flag, M);
}

FLINT_FORCE_INLINE void
DF(_level_sqr)(double * res, const double * x, int M)
{
    int flag;
    SCALAR_LEVEL_OP(flag, _sqr_kernel, M, res, NULL, x, 0, 0);
    DF(_level_finish)(res, flag, M);
}

/* zero components before nonzero ones (after exact cancellation of an
   order sum): squeeze them out, so that the components decrease
   geometrically (each is a rounding error of the sum of the previous
   ones) and the tapered kernels can rely on their positions */
FLINT_FORCE_INLINE void
DF(_shift_zero_head)(double * res)
{
    int i, j = 0;
    for (i = 0; i < N; i++)
        if (res[i] != 0.0)
            res[j++] = res[i];
    for (; j < N; j++)
        res[j] = 0.0;
}

/* nonzero iff |t[i+1]| <= ulp(t[i])/2 for all i (with a tie possibly
   counted as a failure), which implies the canonical property
   |t[i+1]| <= ulp(t[i]) with zeros only at the end (t[i] = 0 forces
   t[i+1] = 0). */
FLINT_FORCE_INLINE int
DF(_is_canonical)(const double * t)
{
    int i, ok = 1;
    UNROLL
    for (i = 0; i + 1 < N; i++)
        ok &= (t[i] + t[i + 1] == t[i]);
    return ok;
}

/* Strong renormalization: res = x exactly, with the components in
   canonical form, |res[i+1]| <= ulp(res[i]) and zeros only at the end.
   The input is any expansion overlapping by at most a few bits (in
   particular any output of the kernels).

   A top-down VecSumErr chain (TwoSums, so that the value is preserved
   exactly whatever the input) produces the canonical form except in
   rare cases (about 1e-5 of the results of the kernels for N = 3 and
   1e-4 for N = 4) where a tiny component is followed by a larger one;
   those are detected by _is_canonical and repaired by a
   VecSumErrBranch pass (which merges components until a nonzero
   residual appears), repeated while needed (one pass has always
   sufficed in testing).  A preceding bottom-up VecSum pass was found
   not to reduce the failure rate and is omitted. */
static void
DF(_canonicalise_slow)(double * t)
{
    int pass;

    for (pass = 0; pass < 8 && !DF(_is_canonical)(t); pass++)
    {
        double eps, f, eps_new, r[N];
        int i, j = 0;

        eps = t[0];
        for (i = 1; i < N; i++)
        {
            DFLOAT_TWO_SUM(f, eps_new, eps, t[i]);
            if (eps_new != 0.0)
            {
                r[j++] = f;
                eps = eps_new;
            }
            else
                eps = f;
        }
        r[j++] = eps;
        for (; j < N; j++)
            r[j] = 0.0;
        for (i = 0; i < N; i++)
            t[i] = r[i];
    }
}

FLINT_FORCE_INLINE void
DF(_canonicalise)(double * t)
{
    int i;

    if (N > 2)
    {
        UNROLL
        for (i = 0; i < N - 1; i++)
            DFLOAT_TWO_SUM(t[i], t[i + 1], t[i], t[i + 1]);
    }
    if (N > 1 && FLINT_UNLIKELY(!DF(_is_canonical)(t)))
        DF(_canonicalise_slow)(t);
}

FLINT_FORCE_INLINE void
DF(_add_kernel)(double * res, double * err, const double * x, const double * y, int track)
{
    if (FLINT_UNLIKELY(DF(_s_addsub_kernel)(res, err, x, y, 0, track, DFLOAT_KERNEL_STRONG)))
        DF(_shift_zero_head)(res);
}

FLINT_FORCE_INLINE void
DF(_sub_kernel)(double * res, double * err, const double * x, const double * y, int track)
{
    if (FLINT_UNLIKELY(DF(_s_addsub_kernel)(res, err, x, y, 1, track, DFLOAT_KERNEL_STRONG)))
        DF(_shift_zero_head)(res);
}

FLINT_FORCE_INLINE void
DF(_mul_kernel)(double * res, double * err, const double * x, const double * y, int track)
{
    if (FLINT_UNLIKELY(DF(_s_mul_kernel)(res, err, x, y, track, DFLOAT_KERNEL_STRONG)))
        DF(_shift_zero_head)(res);
}

FLINT_FORCE_INLINE void
DF(_sqr_kernel)(double * res, double * err, const double * x, int track)
{
    if (FLINT_UNLIKELY(DF(_s_sqr_kernel)(res, err, x, track, DFLOAT_KERNEL_STRONG)))
        DF(_shift_zero_head)(res);
}

FLINT_FORCE_INLINE void
DF(_mul_d_kernel)(double * res, double * err, const double * x, double y, int track)
{
    if (FLINT_UNLIKELY(DF(_s_mul_d_kernel)(res, err, x, y, track, DFLOAT_KERNEL_STRONG)))
        DF(_shift_zero_head)(res);
}

/* ---- SIMD instantiation (4 elements per vector) ---- */

/*
    The vector code is written against machine_vectors.h and needs vec4d
    with fused multiply-adds. It is enabled with AVX2 and FMA and with
    NEON on AArch64. The generic backends can run it too, but whether
    that beats the scalar code depends on how well the compiler maps the
    lane-wise code (and fma) onto the target, so there it is opt-in, as
    it is with MSVC (whose /arch:AVX2 implies FMA without defining
    __FMA__, and which does not define __ARM_NEON on ARM64; the vector
    code has not been tested with MSVC). Define
    DFLOAT_FORCE_SIMD to use the vector code with whatever backend
    machine_vectors.h selects (together with
    FLINT_MACHINE_VECTORS_FORCE_GENERIC this tests the generic code on
    any machine), or DFLOAT_NO_SIMD to use only the scalar code.
*/
#if !defined(DFLOAT_NO_SIMD) && (defined(DFLOAT_FORCE_SIMD) \
    || (!defined(FLINT_MACHINE_VECTORS_GENERIC) \
        && ((defined(__AVX2__) && defined(__FMA__)) \
            || defined(__ARM_NEON))))
#define DFLOAT_HAVE_SIMD 1

#define K_ADD(a, b) vec4d_add((a), (b))
#define K_SUB(a, b) vec4d_sub((a), (b))
#define K_MUL(a, b) vec4d_mul((a), (b))
#define K_FMS(a, b, c) vec4d_fmsub_fused((a), (b), (c))
#define K_NEG(a) vec4d_neg(a)
#define K_ABS(a) vec4d_abs(a)
#define K_ZERO vec4d_zero()
#define K_TWICE(a) vec4d_add((a), (a))
#define K_SET_D(d) vec4d_set_d(d)
#define K_BOUND_BUMP_IF(a) vec4d_bit_and(vec4d_set_d(K_BOUND_BUMP), vec4d_cmp_ne((a), vec4d_zero()))
#define K_ANY_ZERO(a) (vec4d_movemask(vec4d_cmp_eq((a), vec4d_zero())) != 0)
#define K_ANY_GAP(a, b) (vec4d_movemask(vec4d_bit_and(vec4d_cmp_eq((a), vec4d_zero()), vec4d_cmp_ne((b), vec4d_zero()))) != 0)
#define KV vec4d
#define KM N
#define KN(name) DF(_v##name)
#include "kernels.inc"
#undef KM
#undef KN
#if N > 1
#define KM 1
#define KN(name) DF(_v1##name)
#include "kernels.inc"
#undef KM
#undef KN
#endif
#if N > 2
#define KM 2
#define KN(name) DF(_v2##name)
#include "kernels.inc"
#undef KM
#undef KN
#endif
#if N > 3
#define KM 3
#define KN(name) DF(_v3##name)
#include "kernels.inc"
#undef KM
#undef KN
#endif
#undef KV
#undef K_ADD
#undef K_SUB
#undef K_MUL
#undef K_FMS
#undef K_NEG
#undef K_ABS
#undef K_ZERO
#undef K_TWICE
#undef K_SET_D
#undef K_BOUND_BUMP_IF
#undef K_ANY_ZERO
#undef K_ANY_GAP


#define VANY(mask) (vec4d_movemask(mask) != 0)

#define V_TWO_SUM(s, e, a, b) \
    do { \
        vec4d _a = (a), _b = (b), _bb; \
        (s) = vec4d_add(_a, _b); \
        _bb = vec4d_sub((s), _a); \
        (e) = vec4d_add(vec4d_sub(_a, vec4d_sub((s), _bb)), vec4d_sub(_b, _bb)); \
    } while (0)

#define V_NONZERO(a) vec4d_cmp_ne((a), vec4d_zero())

/* 4 consecutive balls (or plain values, stride N) to and from SoA */
FLINT_FORCE_INLINE void
DF(_vload4)(vec4d * m, const double * x, int stride)
{
    int i;
    UNROLL
    for (i = 0; i < N; i++)
        m[i] = vec4d_set_d4(x[i], x[stride + i], x[2 * stride + i], x[3 * stride + i]);
}

/* the same for 4 consecutive elements in reverse order: lane k gets
   element 3 - k (for the reversed dot products) */
FLINT_FORCE_INLINE void
DF(_vload4_rev)(vec4d * m, const double * x, int stride)
{
    int i;
    UNROLL
    for (i = 0; i < N; i++)
        m[i] = vec4d_set_d4(x[3 * stride + i], x[2 * stride + i], x[stride + i], x[i]);
}

FLINT_FORCE_INLINE void
DF(_vstore4)(double * res, const vec4d * m, int stride)
{
    int i;
    UNROLL
    for (i = 0; i < N; i++)
    {
        double t[4];
        vec4d_store_unaligned(t, m[i]);
        res[i] = t[0];
        res[stride + i] = t[1];
        res[2 * stride + i] = t[2];
        res[3 * stride + i] = t[3];
    }
}

FLINT_FORCE_INLINE vec4d
DF(_vabs_sum)(const vec4d * m)
{
    vec4d s = vec4d_abs(m[0]);
    int i;
    UNROLL
    for (i = 1; i < N; i++)
        s = vec4d_add(s, vec4d_abs(m[i]));
    return s;
}

/* mask of lanes with a nonzero component below the underflow guard */
FLINT_FORCE_INLINE vec4d
DF(_vtiny_mask)(const vec4d * m)
{
    vec4d bad = vec4d_zero(), zero = vec4d_zero(), thr = vec4d_set_d(TINY_THRESHOLD);
    int i;
    UNROLL
    for (i = 0; i < N; i++)
    {
        vec4d a = vec4d_abs(m[i]);
        bad = vec4d_bit_or(bad, vec4d_bit_and(vec4d_cmp_gt(a, zero), vec4d_cmp_gt(thr, a)));
    }
    return bad;
}

/* radius product a b with the semantics of _dfloat_rad_mul; sets
   *bad if some lane would need the underflow bump */
FLINT_FORCE_INLINE vec4d
DF(_vrad_mul)(vec4d a, vec4d b, vec4d * bad)
{
    vec4d zero = vec4d_zero(), thr = vec4d_set_d(DFLOAT_RAD_UNDERFLOW_THRESHOLD);
    vec4d p = vec4d_mul(a, b);
    vec4d zin = vec4d_bit_or(vec4d_cmp_eq(a, zero), vec4d_cmp_eq(b, zero));
    vec4d small = vec4d_cmp_nge(p, thr);   /* p < thr or nan */
    *bad = vec4d_bit_or(*bad, vec4d_bit_andnot(zin, small));
    return vec4d_bit_andnot(zin, p);   /* zero where an input is zero */
}

/* 4 ball sums / products in SoA form: returns nonzero if some lane
   hit a special case (zero head after cancellation, radius underflow,
   tiny components, nonfinite result), in which case the results are
   unspecified and the caller redoes the block in scalar */
FLINT_FORCE_INLINE int
DB(_vaddsub4_soa)(vec4d * mr, vec4d * rrad, const vec4d * mx, vec4d rx,
    const vec4d * my, vec4d ry, int negate)
{
    vec4d err = vec4d_zero(), rad, chk;
    int flag;

    flag = DF(_v_addsub_kernel)(mr, &err, mx, my, negate, 1, DFLOAT_KERNEL_STRONG);
    rad = vec4d_mul(vec4d_add(vec4d_add(rx, ry), err), vec4d_set_d(DFLOAT_RAD_SLACK));
    chk = vec4d_add(DF(_vabs_sum)(mr), rad);
    *rrad = rad;
    return flag || VANY(vec4d_cmp_nle(chk, vec4d_set_d(DBL_MAX)));
}

/* mode: see DB(_mul_m); in the automatic mode a block with both exact
   and inexact lanes is redone in scalar */
FLINT_FORCE_INLINE int
DB(_vmul4_soa)(vec4d * mr, vec4d * rrad, const vec4d * mx, vec4d rx,
    const vec4d * my, vec4d ry, int mode)
{
    vec4d err = vec4d_zero(), rad, chk, bad, ax, ay;
    int flag, exact;

    bad = vec4d_bit_or(DF(_vtiny_mask)(mx), DF(_vtiny_mask)(my));
    ax = DF(_vabs_sum)(mx);
    ay = DF(_vabs_sum)(my);
    rad = vec4d_add(DF(_vrad_mul)(ax, ry, &bad), DF(_vrad_mul)(ay, rx, &bad));
    rad = vec4d_add(rad, DF(_vrad_mul)(rx, ry, &bad));
    if (FLINT_UNLIKELY(VANY(bad)))
        return 1;

    if (mode == DFLOAT_MODE_FAST)
        exact = 0;
    else
    {
        exact = vec4d_movemask(vec4d_cmp_eq(rad, vec4d_zero()));
        if (FLINT_UNLIKELY(exact != 0 && exact != 15))
            return 1;
    }
    if (exact)
        flag = DF(_v_mul_kernel)(mr, &err, mx, my, 1, DFLOAT_KERNEL_STRONG);
    else
        flag = DF(_v_mul_kernel)(mr, &err, mx, my, 2, DFLOAT_KERNEL_STRONG);
    rad = vec4d_mul(vec4d_add(rad, err), vec4d_set_d(DFLOAT_RAD_SLACK));
    chk = vec4d_add(DF(_vabs_sum)(mr), rad);
    *rrad = rad;
    return flag || VANY(vec4d_cmp_nle(chk, vec4d_set_d(DBL_MAX)));
}

FLINT_FORCE_INLINE void
DB(_vload_rad)(vec4d * rad, DB(srcptr) x)
{
    *rad = vec4d_set_d4(x[0].rad, x[1].rad, x[2].rad, x[3].rad);
}

/* lanes whose computed radius is below the underflow threshold while it
   should be nonzero (some ingredient nonzero), or nonzero: the vector
   radius arithmetic (plain products and quotients) may have lost or
   underflowed a tiny radius, which the scalar functions bound with
   bumps */
FLINT_FORCE_INLINE vec4d
DB(_v_rad_underflow)(vec4d rad, vec4d inexact)
{
    return vec4d_bit_and(vec4d_cmp_lt(rad, vec4d_set_d(DFLOAT_RAD_UNDERFLOW_THRESHOLD)),
        vec4d_bit_or(inexact, vec4d_cmp_ne(rad, vec4d_zero())));
}

FLINT_FORCE_INLINE void
DB(_vload_rad_rev)(vec4d * rad, DB(srcptr) x)
{
    *rad = vec4d_set_d4(x[3].rad, x[2].rad, x[1].rad, x[0].rad);
}

FLINT_FORCE_INLINE void
DB(_vstore4)(DB(ptr) res, const vec4d * mr, vec4d rad)
{
    double t[4];
    DF(_vstore4)(res->d, mr, N + 1);
    vec4d_store_unaligned(t, rad);
    res[0].rad = t[0]; res[1].rad = t[1]; res[2].rad = t[2]; res[3].rad = t[3];
}

/* AoS block operations: returns 0 if the SIMD result was stored, 1 if
   the block must be redone in scalar */
FLINT_FORCE_INLINE int
DB(_vaddsub4)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, int negate)
{
    vec4d mx[N], my[N], mr[N], rx, ry, rad;
    DF(_vload4)(mx, x->d, N + 1);
    DF(_vload4)(my, y->d, N + 1);
    DB(_vload_rad)(&rx, x);
    DB(_vload_rad)(&ry, y);
    if (FLINT_UNLIKELY(DB(_vaddsub4_soa)(mr, &rad, mx, rx, my, ry, negate)))
        return 1;
    DB(_vstore4)(res, mr, rad);
    return 0;
}

FLINT_FORCE_INLINE int
DB(_vmul4)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, int mode)
{
    vec4d mx[N], my[N], mr[N], rx, ry, rad;
    DF(_vload4)(mx, x->d, N + 1);
    DF(_vload4)(my, y->d, N + 1);
    DB(_vload_rad)(&rx, x);
    DB(_vload_rad)(&ry, y);
    if (FLINT_UNLIKELY(DB(_vmul4_soa)(mr, &rad, mx, rx, my, ry, mode)))
        return 1;
    DB(_vstore4)(res, mr, rad);
    return 0;
}

/* plain: 4 at a time, nonfinite results folded into the head by the
   scalar fallback */
FLINT_FORCE_INLINE int
DF(_vaddsub4)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, int negate)
{
    vec4d mx[N], my[N], mr[N], chk;
    int flag;
    DF(_vload4)(mx, x->d, N);
    DF(_vload4)(my, y->d, N);
    flag = DF(_v_addsub_kernel)(mr, NULL, mx, my, negate, 0, DFLOAT_KERNEL_STRONG);
    chk = DF(_vabs_sum)(mr);
    if (FLINT_UNLIKELY(flag || VANY(vec4d_cmp_nle(chk, vec4d_set_d(DBL_MAX)))))
        return 1;
    DF(_vstore4)(res->d, mr, N);
    return 0;
}

FLINT_FORCE_INLINE int
DF(_vmul4)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y)
{
    vec4d mx[N], my[N], mr[N], chk;
    int flag;
    DF(_vload4)(mx, x->d, N);
    DF(_vload4)(my, y->d, N);
    flag = DF(_v_mul_kernel)(mr, NULL, mx, my, 0, DFLOAT_KERNEL_STRONG);
    chk = DF(_vabs_sum)(mr);
    if (FLINT_UNLIKELY(flag || VANY(vec4d_cmp_nle(chk, vec4d_set_d(DBL_MAX)))))
        return 1;
    DF(_vstore4)(res->d, mr, N);
    return 0;
}

#endif

/* nonzero iff some component is nonzero and below the threshold;
   integer arithmetic on the bit patterns: |x| in (0, T) iff
   (bits(|x|) - 1) < bits(T) - 1 as unsigned integers */
FLINT_FORCE_INLINE int
DF(_has_tiny)(const double * x)
{
    int bad = 0;
    int i;
    UNROLL
    for (i = 0; i < N; i++)
    {
        double_uint64_u u;
        u.f = x[i];
        bad |= ((u.i & UINT64_C(0x7fffffffffffffff)) - 1) < (UINT64_C(0x21f0000000000000) - 1);
    }
    return bad;
}

/* nonzero iff all components are finite */
FLINT_FORCE_INLINE int
DF(_is_finite)(const double * x)
{
    return _dfloat_abs_sum(x, N) <= DBL_MAX;
}

/* upper bound (up to slack) on |x[1]| + ... + |x[N-1]| */
FLINT_FORCE_INLINE double
DF(_tail_abs_sum)(const double * x)
{
    double s = 0.0;
    int i;
    UNROLL
    for (i = 1; i < N; i++)
        s += fabs(x[i]);
    return s;
}

/* for the plain type: fold a nonfinite result into the head */
FLINT_FORCE_INLINE void
DF(_fix_nonfinite)(double * r)
{
    if (FLINT_UNLIKELY(!(_dfloat_abs_sum(r, N) <= DBL_MAX)))
    {
        double s = r[0];
        int i;
        UNROLL
        for (i = 1; i < N; i++)
            s += r[i];
        r[0] = s;
        UNROLL
        for (i = 1; i < N; i++)
            r[i] = 0.0;
    }
}

/* the same with the value to use when the result is nonfinite (the
   double operation on the heads, which distinguishes inf + 1 = inf
   from the nan that the error-free transformations produce) */
FLINT_FORCE_INLINE void
DF(_fix_nonfinite_v)(double * r, double v)
{
    if (FLINT_UNLIKELY(!(_dfloat_abs_sum(r, N) <= DBL_MAX)))
    {
        int i;
        r[0] = v;
        UNROLL
        for (i = 1; i < N; i++)
            r[i] = 0.0;
    }
}

/* for balls: a result with a nonfinite midpoint component or radius
   (or one whose midpoint plus radius overflows) becomes the whole real
   line in its normal form [0 +/- inf] */
FLINT_FORCE_INLINE void
DB(_fix_nonfinite)(double * r, double * rad)
{
    if (FLINT_UNLIKELY(!(_dfloat_abs_sum(r, N) + *rad <= DBL_MAX)))
    {
        int i;
        UNROLL
        for (i = 0; i < N; i++)
            r[i] = 0.0;
        *rad = D_INF;
    }
}

/* whether x is the whole real line (any nonfinite component) */
FLINT_FORCE_INLINE int
DB(_is_W)(DB(srcptr) x)
{
    return _dfloat_ball_is_whole(x->d, N, x->rad);
}

/* Non-inline entry points to the kernels for the generic (runtime N)
   code in scaled.c, which dispatches on N (the general renormalization
   _dfloat_renorm is implemented there directly, for all N). */
void DFV(add_tracked)(double * res, double * err, const double * x, const double * y)
{
    DF(_add_kernel)(res, err, x, y, 1);
}

void DFV(sub_tracked)(double * res, double * err, const double * x, const double * y)
{
    DF(_sub_kernel)(res, err, x, y, 1);
}

void DFV(mul_tracked)(double * res, double * err, const double * x, const double * y)
{
    DF(_mul_kernel)(res, err, x, y, 1);
}

void DFV(sqr_tracked)(double * res, double * err, const double * x)
{
    DF(_sqr_kernel)(res, err, x, 1);
}

/* ------------------------------------------------------------------ */
/* Plain arithmetic                                                    */
/* ------------------------------------------------------------------ */

FLINT_FORCE_INLINE void
DF(_set_d_inline)(DF(t) res, double x)
{
    int i;
    res->d[0] = x;
    for (i = 1; i < N; i++)
        res->d[i] = 0.0;
}
void DF(set_d)(DF(t) res, double x) { DF(_set_d_inline)(res, x); }

void DF(set_dd)(DF(t) res, double x0, double x1)
{
    double t[2] = { x0, x1 };
    DF(renorm)(res, t, 2);
}

FLINT_FORCE_INLINE void DF(_zero_inline)(DF(t) res) { DF(_set_d_inline)(res, 0.0); }
void DF(zero)(DF(t) res) { DF(_zero_inline)(res); }
FLINT_FORCE_INLINE void DF(_one_inline)(DF(t) res) { DF(_set_d_inline)(res, 1.0); }
void DF(one)(DF(t) res) { DF(_one_inline)(res); }

FLINT_FORCE_INLINE void
DF(_set_inline)(DF(t) res, DF(srcptr) x)
{
    *res = *x;
}
void DF(set)(DF(t) res, DF(srcptr) x) { DF(_set_inline)(res, x); }

FLINT_FORCE_INLINE void
DF(_neg_inline)(DF(t) res, DF(srcptr) x)
{
    int i;
    for (i = 0; i < N; i++)
        res->d[i] = -x->d[i];
}
void DF(neg)(DF(t) res, DF(srcptr) x) { DF(_neg_inline)(res, x); }

FLINT_FORCE_INLINE void
DF(_abs_inline)(DF(t) res, DF(srcptr) x)
{
    if (x->d[0] < 0.0 || (x->d[0] == 0.0 && signbit(x->d[0])))
        DF(_neg_inline)(res, x);
    else
        DF(_set_inline)(res, x);
}
void DF(abs)(DF(t) res, DF(srcptr) x) { DF(_abs_inline)(res, x); }

FLINT_FORCE_INLINE void
DF(_canonicalise_inplace)(double * x)
{
    if (FLINT_LIKELY(fabs(x[0]) <= DBL_MAX))
        DF(_canonicalise)(x);
}

FLINT_FORCE_INLINE void
DF(_canonicalise_inline)(DF(t) res, DF(srcptr) x)
{
    DF(_set_inline)(res, x);
    DF(_canonicalise_inplace)(res->d);
}
void DF(canonicalise)(DF(t) res, DF(srcptr) x) { DF(_canonicalise_inline)(res, x); }

DFLOAT_ARITH_IMPL void
DF(_add_impl)(DF(t) res, DF(srcptr) x, DF(srcptr) y)
{
    double v = x->d[0] + y->d[0];
    DF(_add_kernel)(res->d, NULL, x->d, y->d, 0);
    DF(_fix_nonfinite_v)(res->d, v);
}
void DF(add)(DF(t) res, DF(srcptr) x, DF(srcptr) y) { DF(_add_impl)(res, x, y); }

DFLOAT_ARITH_IMPL void
DF(_sub_impl)(DF(t) res, DF(srcptr) x, DF(srcptr) y)
{
    double v = x->d[0] - y->d[0];
    DF(_sub_kernel)(res->d, NULL, x->d, y->d, 0);
    DF(_fix_nonfinite_v)(res->d, v);
}
void DF(sub)(DF(t) res, DF(srcptr) x, DF(srcptr) y) { DF(_sub_impl)(res, x, y); }

DFLOAT_ARITH_IMPL void
DF(_mul_impl)(DF(t) res, DF(srcptr) x, DF(srcptr) y)
{
    double v = x->d[0] * y->d[0];
    DF(_mul_kernel)(res->d, NULL, x->d, y->d, 0);
    DF(_fix_nonfinite_v)(res->d, v);
}
void DF(mul)(DF(t) res, DF(srcptr) x, DF(srcptr) y) { DF(_mul_impl)(res, x, y); }

DFLOAT_ARITH_IMPL void
DF(_sqr_impl)(DF(t) res, DF(srcptr) x)
{
    double v = x->d[0] * x->d[0];
    DF(_sqr_kernel)(res->d, NULL, x->d, 0);
    DF(_fix_nonfinite_v)(res->d, v);
}
void DF(sqr)(DF(t) res, DF(srcptr) x) { DF(_sqr_impl)(res, x); }

DFLOAT_ARITH_IMPL void
DF(_mul_d_impl)(DF(t) res, DF(srcptr) x, double y)
{
    double v = x->d[0] * y;
    DF(_mul_d_kernel)(res->d, NULL, x->d, y, 0);
    DF(_fix_nonfinite_v)(res->d, v);
}
void DF(mul_d)(DF(t) res, DF(srcptr) x, double y) { DF(_mul_d_impl)(res, x, y); }

void DF(renorm)(DF(t) res, const double * t, slong len)
{
    _dfloat_renorm(res->d, N, NULL, (double *) t, len);
}

/* a compiled lower-precision function (dM_name) at level M < N on
   the first M components; used for the division inside the square
   root (M is a constant after unrolling) */
#define LEVEL_CALL(name, M, ...) \
    do { \
        switch (M) \
        { \
            case 1: d1_##name((d1_ptr) __VA_ARGS__); break; \
            case 2: d2_##name((d2_ptr) __VA_ARGS__); break; \
            default: d3_##name((d3_ptr) __VA_ARGS__); break; \
        } \
    } while (0)

FLINT_FORCE_INLINE void
DF(_level_add3)(double * res, const double * x, const double * y, const double * z, int M)
{
    int flag;
    SCALAR_LEVEL_OP(flag, _add3_kernel, M, res, NULL, x, y, z, 0, 0);
    DF(_level_finish)(res, flag, M);
}

/* Long division at level L with D <= L digits: the quotient digits
   q_k = r_0 (1/y_0) (within 2u + u^2 of r_0 / y_0, which the next
   remainder absorbs) and the remainders r <- r - q_k y.  The
   remainder after step k is about u^(k+1) |x| and is only needed to
   absolute precision u^L |x|, so it is kept with L - k - 1
   components (tapered precision).  The cancellation in the head is
   done exactly: q_k y_0 = p + e (TwoProd), d_0 = r_0 - p (exact by
   Sterbenz's lemma), and r - q_k y = (d_0 - e) + (r_1 + ... ) -
   q_k (y_1 + ...) is a three-term sum of expansions of the size of
   the new remainder, so that its components are properly tapered (a
   subtraction with a cancelling head would leave the information
   spread over one more component than the level keeps).  The D
   digits are renormalized into q[0..D-1] (the rest is zeroed). */
FLINT_FORCE_INLINE void
DF(_div_level)(double * q, const double * x, const double * y, int L, int D)
{
    /* (for N = 1, X[1] appears in code that is dead) */
    double r[N], X[N == 1 ? 2 : N], Z[N], qk[N], y1[N], yinv, p, e, d0;
    int i, k;

    for (i = 0; i < N; i++)
        r[i] = (i < L) ? x[i] : 0.0;
    /* (y_1, ..., y_{N-1}, 0): the kernels at level M - 1 <= L - 1 read
       only y_1, ..., y_{M-1}, but padding to N components makes this
       visible to the compiler (no out-of-bounds warnings) */
    for (i = 0; i < N; i++)
        y1[i] = (i + 1 < N) ? y[i + 1] : 0.0;
    yinv = 1.0 / y[0];

    UNROLL
    for (k = 0; k < D; k++)
    {
        const int M = L - k;
        qk[k] = (M == 1) ? r[0] / y[0] : r[0] * yinv;
        if (k + 1 < D)
        {
            DFLOAT_TWO_PROD(p, e, qk[k], y[0]);
            d0 = r[0] - p;
            DFLOAT_TWO_SUM(X[0], X[1], d0, -e);
            for (i = 2; i < N; i++)
                X[i] = 0.0;
            /* Y = (r_1, ..., r_{M-1}): shift r down in place */
            for (i = 0; i + 1 < M; i++)
                r[i] = r[i + 1];
            for (i = 0; i < N; i++)
                if (i >= M - 1)
                    r[i] = 0.0;
            DF(_level_mul_d)(Z, y1, -qk[k], M - 1);
            DF(_level_add3)(r, X, r, Z, M - 1);
        }
    }

    for (i = 0; i < N; i++)
        if (i >= D)
            qk[i] = 0.0;

    if (D == N)
        DF(_renorm_terms)(q, NULL, qk, N, 0);
    else
    {
        int flag;
        SCALAR_LEVEL_OP(flag, _renorm_orders, D, q, qk, 0);
        DF(_level_finish)(q, flag, D);
    }
}

/* The kernel of the division is inlined for N <= 2; for N >= 3, where
   it is large (1-2 KB) and expensive, it is a single static copy shared
   by the plain division, the approximate entry point and the elementary
   functions that use it (the ball division inlines its own). */
#if N <= 2
#define DFLOAT_DIV_KERNEL FLINT_FORCE_INLINE
#else
#define DFLOAT_DIV_KERNEL static DFLOAT_NOINLINE
#endif

DFLOAT_DIV_KERNEL void
DF(_div_kernel)(double * q, const double * x, const double * y)
{
    DF(_div_level)(q, x, y, N, N);
}

static DFLOAT_NOINLINE void
DF(_div_ool)(DF(t) res, DF(srcptr) x, DF(srcptr) y)
{
    if (y->d[0] == 0.0 || !DF(_is_finite)(x->d) || !DF(_is_finite)(y->d))
    {
        DF(_set_d_inline)(res, x->d[0] / y->d[0]);
    }
    else
    {
        DF(_div_kernel)(res->d, x->d, y->d);
        DF(_fix_nonfinite)(res->d);
    }
}

FLINT_FORCE_INLINE void
DF(_div_inline)(DF(t) res, DF(srcptr) x, DF(srcptr) y)
{
    if (N == 1)
        res->d[0] = x->d[0] / y->d[0];
    else
        DF(_div_ool)(res, x, y);
}
void DF(div)(DF(t) res, DF(srcptr) x, DF(srcptr) y) { DF(_div_inline)(res, x, y); }

FLINT_FORCE_INLINE void
DF(_inv_inline)(DF(t) res, DF(srcptr) x)
{
    DF(t) one;
    DF(_one_inline)(one);
    DF(_div_inline)(res, one, x);
}
void DF(inv)(DF(t) res, DF(srcptr) x) { DF(_inv_inline)(res, x); }

/* Newton iteration for the square root: s <- s + (x - s^2) / (2 s),
   starting from the double square root and doubling the number of
   correct components each step (M -> Mn = min(2M, N)).  With s
   correct to M components, t = x - s^2 is about u^M |x|; it is
   computed as (x_0 - s_0^2) - (s - s_0)(s + s_0) + (x_1 + ...) with
   the head exact (TwoProd and Sterbenz) and the rest as a three-term
   sum at level Mn - 1, which is accurate to u^Mn |x| but, because of
   the cancellation, spreads the information over all Mn - 1
   components; the correction t / (2 s) is therefore a division at
   level Mn - 1 producing only the M digits it needs.  x > 0 is
   assumed. */
FLINT_FORCE_INLINE void
DF(_sqrt_kernel)(double * s, const double * x)
{
    /* (for N = 1, X[1] and t[1] appear in code that is dead) */
    double t[N == 1 ? 2 : N], u[N], X[N == 1 ? 2 : N], A[N], B[N], P[N], two_s[N], p, e, d0;
    int i, M, Mn;

    s[0] = sqrt(x[0]);
    for (i = 1; i < N; i++)
        s[i] = 0.0;

    UNROLL
    for (M = 1; M < N; M = Mn)
    {
        Mn = FLINT_MIN(2 * M, N);
        if (M == 1)
        {
            /* t = x - s_0^2 exactly in the head, then t / (2 s_0) */
            DFLOAT_TWO_PROD(p, e, s[0], s[0]);
            d0 = x[0] - p;
            DFLOAT_TWO_SUM(X[0], X[1], d0, -e);
            for (i = 2; i < N; i++)
                X[i] = 0.0;
            for (i = 0; i < N; i++)
                A[i] = (i + 1 < N && i + 1 < Mn) ? x[i + 1] : 0.0;
            if (Mn - 1 == 1)
                t[0] = X[0] + A[0], t[1] = 0.0;
            else
                DF(_level_addsub)(t, X, A, 0, Mn - 1);
            u[0] = t[0] / (2.0 * s[0]);
            for (i = 1; i < N; i++)
                u[i] = 0.0;
        }
        else
        {
            DFLOAT_TWO_PROD(p, e, s[0], s[0]);
            d0 = x[0] - p;
            DFLOAT_TWO_SUM(X[0], X[1], d0, -e);
            for (i = 2; i < N; i++)
                X[i] = 0.0;
            /* (s - s_0)(s + s_0) = (s_1 + ...)(2 s_0 + s_1 + ...) */
            for (i = 0; i < N; i++)
            {
                A[i] = (i + 1 < N && i + 1 < M) ? s[i + 1] : 0.0;
                B[i] = (i == 0) ? 2.0 * s[0] : ((i < M) ? s[i] : 0.0);
            }
            DF(_level_mul)(P, A, B, Mn - 1);
            for (i = 0; i < N; i++)
            {
                P[i] = -P[i];
                A[i] = (i + 1 < N && i + 1 < Mn) ? x[i + 1] : 0.0;
            }
            DF(_level_add3)(t, X, A, P, Mn - 1);
            for (i = 0; i < N; i++)
                two_s[i] = (i < M) ? 2.0 * s[i] : 0.0;
            DF(_div_level)(u, t, two_s, Mn - 1, M);
        }
        DF(_level_addsub)(s, s, u, 0, Mn);
    }
}

static DFLOAT_NOINLINE void
DF(_sqrt_ool)(DF(t) res, DF(srcptr) x)
{
    if (!(x->d[0] > 0.0) || !DF(_is_finite)(x->d))
    {
        DF(_set_d_inline)(res, sqrt(x->d[0]));
    }
    else
    {
        /* the kernel writes the output before it has read all of the
           input: a copy for res == x */
        double xc[N];
        memcpy(xc, x->d, sizeof(xc));
        DF(_sqrt_kernel)(res->d, xc);
        DF(_fix_nonfinite)(res->d);
    }
}

FLINT_FORCE_INLINE void
DF(_sqrt_inline)(DF(t) res, DF(srcptr) x)
{
    if (N == 1)
        res->d[0] = sqrt(x->d[0]);
    else
        DF(_sqrt_ool)(res, x);
}
void DF(sqrt)(DF(t) res, DF(srcptr) x) { DF(_sqrt_inline)(res, x); }

void DFV(div_approx)(double * q, const double * x, const double * y)
{
    DF(_div_kernel)(q, x, y);
}

/* (through the full function, which is the kernel for x > 0, so that
   the kernel is not compiled a second time) */
void DFV(sqrt_approx)(double * s, const double * x)
{
    DF(_sqrt_inline)((DF(ptr)) s, (DF(srcptr)) x);
}

/* Newton iteration for the reciprocal square root: r <- r + r (1 - x
   r^2) / 2, starting from the double reciprocal root and doubling
   the number of correct components each step (M -> Mn = min(2M, N)).
   With r correct to M components, t = 1 - x r^2 is about u^M and is
   needed to absolute accuracy u^Mn.  In the first step (M = 1) it is
   assembled in doubles from the exact products r_0^2 = p + e and
   x_0 p = q + f: t = (1 - q) - f - x_0 e - x_1 (p + e), with 1 - q
   exact; the terms left out are O(u^2).  In the later steps r^2 and
   x r^2 are formed at level Mn (the head of x r^2 is within a few
   ulps of 1) and 1 - x r^2 as a two-term sum at level Mn - 1 whose
   head cancels exactly.  The correction r t / 2 only needs Mn - M
   correct components (it is formed at level Mn - 1, which reads all
   the components of t).  x > 0 is assumed. */
FLINT_FORCE_INLINE void
DF(_rsqrt_kernel)(double * r, const double * x)
{
    double P[N], Q[N], X[N], Y[N], t[N], U[N];
    int i, M, Mn;

    r[0] = 1.0 / sqrt(x[0]);
    for (i = 1; i < N; i++)
        r[i] = 0.0;
    if (N == 1)
        return;

    {
        double p, e, q, f, t0, u0;
        DFLOAT_TWO_PROD(p, e, r[0], r[0]);
        DFLOAT_TWO_PROD(q, f, x[0], p);
        t0 = (1.0 - q) - f;
        t0 = fma(-x[0], e, t0);
        t0 = fma(-x[1], p + e, t0);
        u0 = 0.5 * r[0] * t0;
        DFLOAT_TWO_SUM(r[0], r[1], r[0], u0);
    }

    /* (at most one pass for N <= 4, which the compilers peel without
       an unrolling pragma; clang fails to apply one here) */
    for (M = 2; M < N; M = Mn)
    {
        Mn = FLINT_MIN(2 * M, N);
        DF(_level_sqr)(P, r, Mn);
        DF(_level_mul)(Q, x, P, Mn);
        X[0] = 1.0 - Q[0];
        for (i = 1; i < N; i++)
            X[i] = 0.0;
        for (i = 0; i < N; i++)
            Y[i] = (i + 1 < N && i + 1 < Mn) ? Q[i + 1] : 0.0;
        DF(_level_addsub)(t, X, Y, 1, Mn - 1);
        DF(_level_mul)(U, r, t, Mn - 1);
        for (i = 0; i < N; i++)
            U[i] *= 0.5;
        DF(_level_addsub)(r, r, U, 0, Mn);
    }
}

static DFLOAT_NOINLINE void
DF(_rsqrt_ool)(DF(t) res, DF(srcptr) x)
{
    if (!(x->d[0] > 0.0) || !DF(_is_finite)(x->d))
    {
        DF(_set_d_inline)(res, 1.0 / sqrt(x->d[0]));
    }
    else
    {
        /* the kernel writes the output before it has read all of the
           input: a copy for res == x */
        double xc[N];
        memcpy(xc, x->d, sizeof(xc));
        DF(_rsqrt_kernel)(res->d, xc);
        DF(_fix_nonfinite)(res->d);
    }
}

FLINT_FORCE_INLINE void
DF(_rsqrt_inline)(DF(t) res, DF(srcptr) x)
{
    if (N == 1)
        res->d[0] = 1.0 / sqrt(x->d[0]);
    else
        DF(_rsqrt_ool)(res, x);
}
void DF(rsqrt)(DF(t) res, DF(srcptr) x) { DF(_rsqrt_inline)(res, x); }

void DFV(rsqrt_approx)(double * r, const double * x)
{
    DF(_rsqrt_inline)((DF(ptr)) r, (DF(srcptr)) x);
}

DFLOAT_ARITH_IMPL void
DF(_mul_2exp_si_impl)(DF(t) res, DF(srcptr) x, slong e)
{
    int i;
    if (e >= -1000 && e <= 1000)
    {
        double s = d_mul_2exp(1.0, e);
        for (i = 0; i < N; i++)
            res->d[i] = x->d[i] * s;
    }
    else
    {
        for (i = 0; i < N; i++)
            res->d[i] = ldexp(x->d[i], FLINT_MAX(FLINT_MIN(e, 5000), -5000));
    }
    DF(_fix_nonfinite)(res->d);
}
void DF(mul_2exp_si)(DF(t) res, DF(srcptr) x, slong e) { DF(_mul_2exp_si_impl)(res, x, e); }

FLINT_FORCE_INLINE int DF(_is_zero_inline)(DF(srcptr) x) { return x->d[0] == 0.0; }
int DF(is_zero)(DF(srcptr) x) { return DF(_is_zero_inline)(x); }

FLINT_FORCE_INLINE int
DF(_is_one_inline)(DF(srcptr) x)
{
    int i;
    if (x->d[0] != 1.0)
        return 0;
    for (i = 1; i < N; i++)
        if (x->d[i] != 0.0)
            return 0;
    return 1;
}
int DF(is_one)(DF(srcptr) x) { return DF(_is_one_inline)(x); }

FLINT_FORCE_INLINE int DF(_is_finite_inline)(DF(srcptr) x) { return DF(_is_finite)(x->d); }
int DF(is_finite)(DF(srcptr) x) { return DF(_is_finite_inline)(x); }
FLINT_FORCE_INLINE int DF(_is_nan_inline)(DF(srcptr) x) { return x->d[0] != x->d[0]; }
int DF(is_nan)(DF(srcptr) x) { return DF(_is_nan_inline)(x); }

/* The comparisons and predicates: bounds on the difference x - y (or
   x - c), computed once by a tracked subtraction d + [-e, e]:
   *lo <= |x - y| <= finish(*hi) (lo may be negative; hi is returned
   before the slack factor, for the callers to add to), *d0 the head of d
   (whose sign is that of the difference when lo > 0), and the return
   value is nonzero iff the difference is exactly zero.  Shared out of
   line by all the predicates (it is the only place they inline the
   subtraction kernel). */
static DFLOAT_NOINLINE int
DF(_diff_bounds)(double * lo, double * hi, double * d0, const double * x, const double * y)
{
    double d[N], e = 0.0;
    DF(_sub_kernel)(d, &e, x, y, 1);
    *d0 = d[0];
    *lo = fabs(d[0]) - _dfloat_rad_finish(DF(_tail_abs_sum)(d) + e);
    *hi = _dfloat_abs_sum(d, N) + e;
    return d[0] == 0.0 && e == 0.0;
}

static int
DF(_diff_bounds_d)(double * lo, double * hi, double * d0, const double * x, double c)
{
    double cc[N];
    int i;
    cc[0] = c;
    for (i = 1; i < N; i++)
        cc[i] = 0.0;
    return DF(_diff_bounds)(lo, hi, d0, x, cc);
}

/* exact comparison through a tracked subtraction: the difference is
   d + [-e, e] */
FLINT_FORCE_INLINE int
DF(_cmp_inline)(DF(srcptr) x, DF(srcptr) y)
{
    double lo, hi, d0;
    if (N == 1)
        return (x->d[0] > y->d[0]) - (x->d[0] < y->d[0]);
    DF(_diff_bounds)(&lo, &hi, &d0, x->d, y->d);
    if (d0 == 0.0)
        return 0;   /* renorm keeps the head nonzero when the sum is nonzero */
    if (lo > 0.0)
        return (d0 > 0.0) ? 1 : -1;
    /* not decidable from the truncated difference: fall back to arf */
    {
        arf_t a, b;
        int c;
        arf_init(a); arf_init(b);
        DF(get_arf)(a, x); DF(get_arf)(b, y);
        c = arf_cmp(a, b);
        arf_clear(a); arf_clear(b);
        return c;
    }
}
int DF(cmp)(DF(srcptr) x, DF(srcptr) y) { return DF(_cmp_inline)(x, y); }

FLINT_FORCE_INLINE int
DF(_equal_inline)(DF(srcptr) x, DF(srcptr) y)
{
    if (!DF(_is_finite)(x->d) || !DF(_is_finite)(y->d))
        return x->d[0] == y->d[0];
    return DF(_cmp_inline)(x, y) == 0;
}
int DF(equal)(DF(srcptr) x, DF(srcptr) y) { return DF(_equal_inline)(x, y); }

int DF(sgn)(DF(srcptr) x)
{
    return (x->d[0] > 0.0) - (x->d[0] < 0.0);
}

double DF(get_d)(DF(srcptr) x) { return x->d[0]; }

void DF(get_arf)(arf_t res, DF(srcptr) x) { _dfloat_get_arf(res, x->d, N); }
void DF(set_arf)(DF(t) res, const arf_t x) { _dfloat_set_arf(res->d, N, NULL, x); }
void DF(get_arb)(arb_t res, DF(srcptr) x) { _dfloat_get_arb(res, x->d, N, 0.0); }

void DF(set_arb)(DF(t) res, const arb_t x)
{
    _dfloat_set_arf(res->d, N, NULL, arb_midref(x));
}

void DF(set_fmpz)(DF(t) res, const fmpz_t x)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, x);
    DF(set_arf)(res, t);
    arf_clear(t);
}

void DF(set_fmpq)(DF(t) res, const fmpq_t x)
{
    arb_t t;
    arb_init(t);
    arb_set_fmpq(t, x, 53 * N + 64);
    DF(set_arb)(res, t);
    arb_clear(t);
}

int DF(set_str)(DF(t) res, const char * s)
{
    return _dfloat_set_str(res->d, NULL, N, s);
}

char * DF(get_str)(DF(srcptr) x, slong digits)
{
    return _dfloat_get_str(x->d, N, -1.0, digits);
}

void DF(print)(DF(srcptr) x)
{
    char * s = DF(get_str)(x, 0);
    flint_printf("%s", s);
    flint_free(s);
}

void DF(randtest)(DF(t) res, flint_rand_t state) { _dfloat_randtest(res->d, N, state); }
void DF(randtest_special)(DF(t) res, flint_rand_t state) { _dfloat_randtest_special(res->d, N, state); }

#if N > 1
static int DF(_exp_core)(double * res, const double * x);
#else
FLINT_FORCE_INLINE int DF(_exp_core)(double * res, const double * x);
#endif
FLINT_FORCE_INLINE int DF(_exp_core_m)(double * res, const double * x, int m1);

DFLOAT_FUNC_NOINLINE void DF(exp)(DF(t) res, DF(srcptr) x)
{
    double xc[N], x0;
    int i;

#if N == 1
    x0 = x->d[0];
    if (FLINT_LIKELY(fabs(x0) <= 708.0))
    {
        /* the fast path: a normal result, scaled by an exponent
           addition */
        double_uint64_u u;
        if (FLINT_UNLIKELY(DF(_exp_core)(&u.f, x->d)))
            _dfloat_exp_slow(res->d, N, x->d);
        else
            res->d[0] = u.f;
        return;
    }
    xc[0] = x0;
#else
    for (i = 0; i < N; i++)
        xc[i] = x->d[i];
    DF(_canonicalise_inplace)(xc);
    x0 = xc[0];
#endif

    if (!(fabs(x0) <= DBL_MAX) || x0 > 709.9 || x0 < -746.0)
    {
        /* inf, -inf, nan, overflow and underflow as for doubles */
        res->d[0] = exp(x0);
        for (i = 1; i < N; i++)
            res->d[i] = 0.0;
        return;
    }
    if (FLINT_UNLIKELY(!DF(_is_canonical)(xc) || DF(_exp_core)(res->d, xc)))
        _dfloat_exp_slow(res->d, N, xc);
}

/* expm1: the kernel's variant with the tables of exp(i/256) - 1 for
   |x| <= 0.34 (relative accuracy near 0), exp(x) - 1 beyond (where
   |exp(x) - 1| >= 0.29 and the subtraction loses nothing); the
   guarantee is |res - expm1(x)| <= EPS |expm1(x)| with EPS =
   DFLOAT_EXPM1_EPS_N */
#define EXPM1_XMAX 0.34
/* shared by the plain and ball functions; inlined for N = 1, where the
   call would cost 15% */
#if N >= 2
static DFLOAT_NOINLINE int
#else
FLINT_FORCE_INLINE int
#endif
DF(_expm1_core)(double * res, const double * xc)
{
    if (fabs(xc[0]) <= EXPM1_XMAX)
        return DF(_exp_core_m)(res, xc, 1);
    else
    {
        double y[N], one[N];
        int idx;
        if (DF(_exp_core)(y, xc))
            return 1;
        one[0] = 1.0;
        for (idx = 1; idx < N; idx++)
            one[idx] = 0.0;
        DF(_sub_kernel)(res, NULL, y, one, 0);
        return 0;
    }
}

DFLOAT_FUNC_NOINLINE void DF(expm1)(DF(t) res, DF(srcptr) x)
{
    double xc[N], x0;
    int i;

    for (i = 0; i < N; i++)
        xc[i] = x->d[i];
    DF(_canonicalise_inplace)(xc);
    x0 = xc[0];
    if (!(fabs(x0) <= DBL_MAX) || x0 > 709.9 || x0 < -746.0)
    {
        /* inf, nan, overflow; -1 exactly below -746 (exp(x) < 2^-1075
           is far below the relative accuracy) */
        res->d[0] = (x0 < -746.0) ? -1.0 : expm1(x0);
        for (i = 1; i < N; i++)
            res->d[i] = 0.0;
        return;
    }
    if (FLINT_UNLIKELY(!DF(_is_canonical)(xc) || DF(_expm1_core)(res->d, xc)))
        _dfloat_expm1_slow(res->d, N, xc);
}

/* vectors */

void DFV(vec_add)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
        if (FLINT_UNLIKELY(DF(_vaddsub4)(res + i, x + i, y + i, 0)))
        {
            DF(_add_impl)(res + i, x + i, y + i);
            DF(_add_impl)(res + i + 1, x + i + 1, y + i + 1);
            DF(_add_impl)(res + i + 2, x + i + 2, y + i + 2);
            DF(_add_impl)(res + i + 3, x + i + 3, y + i + 3);
        }
#endif
    for (; i < len; i++)
        DF(_add_impl)(res + i, x + i, y + i);
}

void DFV(vec_sub)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
        if (FLINT_UNLIKELY(DF(_vaddsub4)(res + i, x + i, y + i, 1)))
        {
            DF(_sub_impl)(res + i, x + i, y + i);
            DF(_sub_impl)(res + i + 1, x + i + 1, y + i + 1);
            DF(_sub_impl)(res + i + 2, x + i + 2, y + i + 2);
            DF(_sub_impl)(res + i + 3, x + i + 3, y + i + 3);
        }
#endif
    for (; i < len; i++)
        DF(_sub_impl)(res + i, x + i, y + i);
}

void DFV(vec_mul)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
        if (FLINT_UNLIKELY(DF(_vmul4)(res + i, x + i, y + i)))
        {
            DF(_mul_impl)(res + i, x + i, y + i);
            DF(_mul_impl)(res + i + 1, x + i + 1, y + i + 1);
            DF(_mul_impl)(res + i + 2, x + i + 2, y + i + 2);
            DF(_mul_impl)(res + i + 3, x + i + 3, y + i + 3);
        }
#endif
    for (; i < len; i++)
        DF(_mul_impl)(res + i, x + i, y + i);
}

void DFV(vec_mul_scalar)(DF(ptr) res, DF(srcptr) x, slong len, DF(srcptr) c)
{
    slong i;
    for (i = 0; i < len; i++)
        DF(_mul_impl)(res + i, x + i, c);
}

void DFV(vec_addmul_scalar)(DF(ptr) res, DF(srcptr) x, slong len, DF(srcptr) c)
{
    slong i;
    DF(t) t;
    for (i = 0; i < len; i++)
    {
        DF(_mul_impl)(t, x + i, c);
        DF(_add_impl)(res + i, res + i, t);
    }
}

void DFV(vec_submul_scalar)(DF(ptr) res, DF(srcptr) x, slong len, DF(srcptr) c)
{
    slong i;
    DF(t) t;
    for (i = 0; i < len; i++)
    {
        DF(_mul_impl)(t, x + i, c);
        DF(_sub_impl)(res + i, res + i, t);
    }
}

/* the dot products, forward (rev = 0) or with y reversed (rev = 1, the
   kernel of the classical polynomial multiplication and of the series
   basecases): four lane accumulators in SIMD registers, combined at
   the end; the reversed one loads the blocks of y in reverse lane
   order, so that both run the same loop */
/* one block of 4 terms of a plain dot product into the lane
   accumulators ACC (a block with a special case is redone in scalar) */
#define DFV_DOT_BLOCK(ACC, I) \
    do { \
        vec4d _mx[N], _my[N], _prod[N], _sum[N], _chk; \
        int _k, _flag; \
        DF(_vload4)(_mx, x[I].d, N); \
        if (rev) \
            DF(_vload4_rev)(_my, y[len - 4 - (I)].d, N); \
        else \
            DF(_vload4)(_my, y[I].d, N); \
        _flag = DF(_v_mul_kernel)(_prod, NULL, _mx, _my, 0, DFLOAT_KERNEL_STRONG); \
        _flag |= DF(_v_addsub_kernel)(_sum, NULL, ACC, _prod, 0, 0, DFLOAT_KERNEL_STRONG); \
        _chk = DF(_vabs_sum)(_sum); \
        if (FLINT_LIKELY(!_flag && !VANY(vec4d_cmp_nle(_chk, vec4d_set_d(DBL_MAX))))) \
        { \
            for (_k = 0; _k < N; _k++) \
                ACC[_k] = _sum[_k]; \
        } \
        else \
        { \
            DF(_vstore4)(a4->d, ACC, N); \
            for (_k = 0; _k < 4; _k++) \
            { \
                DF(_mul_impl)(t, x + (I) + _k, rev ? y + len - 1 - (I) - _k : y + (I) + _k); \
                DF(_add_impl)(a4 + _k, a4 + _k, t); \
            } \
            DF(_vload4)(ACC, a4->d, N); \
        } \
    } while (0)

/* the dot products, forward (rev = 0) or with y reversed (rev = 1, the
   kernel of the classical polynomial multiplication and of the series
   basecases): four lane accumulators in SIMD registers, combined at
   the end (a second set of four for small N does not help: the loop is
   bound by its loads); the reversed one loads the blocks of y in reverse lane
   order, so that both run the same loop (and give the same result as
   the forward one on a reversed copy of y) */
FLINT_FORCE_INLINE void
DFV(_vec_dot_m)(DF(t) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, slong len, int rev)
{
    slong i = 0;
    DF(t) s, t;
    if (initial == NULL)
        DF(_zero_inline)(s);
    else
        DF(_set_inline)(s, initial);
#ifdef DFLOAT_HAVE_SIMD
    if (len >= DFLOAT_DOT_PLAIN_SIMD_MIN)
    {
        vec4d acc[N];
        DF(struct) a4[4];
        int k;
        for (k = 0; k < N; k++)
            acc[k] = vec4d_zero();
        for (; i + 4 <= len; i += 4)
            DFV_DOT_BLOCK(acc, i);
        DF(_vstore4)(a4->d, acc, N);
        DF(_add_impl)(a4, a4, a4 + 1);
        DF(_add_impl)(a4 + 2, a4 + 2, a4 + 3);
        DF(_add_impl)(a4, a4, a4 + 2);
        if (subtract)
            DF(_sub_impl)(s, s, a4);
        else
            DF(_add_impl)(s, s, a4);
    }
#endif
    for (; i < len; i++)
    {
        DF(_mul_impl)(t, x + i, rev ? y + len - 1 - i : y + i);
        if (subtract)
            DF(_sub_impl)(s, s, t);
        else
            DF(_add_impl)(s, s, t);
    }
    DF(_set_inline)(res, s);
}

#undef DFV_DOT_BLOCK

/* short dot products (the SIMD loop only starts at 8 terms): a lean
   scalar loop of its own, so that the many short dot products of the
   polynomial basecases do not pay for the SIMD function's frame (rev:
   y reversed) */
/* (inlined into the dispatch for N <= 2, where a call costs 5-10% of a
   short dot product) */
#if N <= 2
FLINT_FORCE_INLINE void
#else
static DFLOAT_NOINLINE void
#endif
DFV(_vec_dot_short)(DF(t) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, int rev, slong len)
{
    slong i;
    DF(t) s, t;
    if (initial == NULL)
        DF(_zero_inline)(s);
    else
        DF(_set_inline)(s, initial);
    /* (no unrolling: knowing len < DFLOAT_DOT_PLAIN_SIMD_MIN, GCC would
       otherwise unroll this fully, which costs a stack frame per call) */
    DFLOAT_NO_UNROLL
    for (i = 0; i < len; i++)
    {
        DF(_mul_impl)(t, x + i, rev ? y + len - 1 - i : y + i);
        if (subtract)
            DF(_sub_impl)(s, s, t);
        else
            DF(_add_impl)(s, s, t);
    }
    DF(_set_inline)(res, s);
}

/* the SIMD bodies return GR_SUCCESS, so that the gr wrappers reach them
   by a tail call (with a plain call after the inlined short loop, the
   wrapper would set up a stack frame on every call, short or long) */
static DFLOAT_NOINLINE int DFV(_vec_dot_long)(DF(t) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, slong len)
{ DFV(_vec_dot_m)(res, initial, subtract, x, y, len, 0); return GR_SUCCESS; }
static DFLOAT_NOINLINE int DFV(_vec_dot_rev_long)(DF(t) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, slong len)
{ DFV(_vec_dot_m)(res, initial, subtract, x, y, len, 1); return GR_SUCCESS; }

FLINT_FORCE_INLINE int
DFV(_vec_dot_dispatch)(DF(t) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, slong len, int rev)
{
    /* the short case falls through (most calls from the polynomial
       basecases are short) */
    if (FLINT_UNLIKELY(len >= DFLOAT_DOT_PLAIN_SIMD_MIN))
    {
        if (rev)
            return DFV(_vec_dot_rev_long)(res, initial, subtract, x, y, len);
        else
            return DFV(_vec_dot_long)(res, initial, subtract, x, y, len);
    }
    if (rev)
        DFV(_vec_dot_short)(res, initial, subtract, x, y, 1, len);
    else
        DFV(_vec_dot_short)(res, initial, subtract, x, y, 0, len);
    return GR_SUCCESS;
}

void DFV(vec_dot)(DF(t) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, slong len)
{ DFV(_vec_dot_dispatch)(res, initial, subtract, x, y, len, 0); }

void DFV(vec_dot_rev)(DF(t) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, slong len)
{ DFV(_vec_dot_dispatch)(res, initial, subtract, x, y, len, 1); }

/* ------------------------------------------------------------------ */
/* Ball arithmetic                                                     */
/* ------------------------------------------------------------------ */

FLINT_FORCE_INLINE void
DB(_indeterminate_inline)(DB(t) res)
{
    int i;
    for (i = 0; i < N; i++)
        res->d[i] = 0.0;
    res->rad = D_INF;
}
void DB(indeterminate)(DB(t) res) { DB(_indeterminate_inline)(res); }

FLINT_FORCE_INLINE void
DB(_set_d_inline)(DB(t) res, double x)
{
    int i;
    if (x == x && fabs(x) <= DBL_MAX)
    {
        res->d[0] = x;
        for (i = 1; i < N; i++)
            res->d[i] = 0.0;
        res->rad = 0.0;
    }
    else
        DB(_indeterminate_inline)(res);
}
void DB(set_d)(DB(t) res, double x) { DB(_set_d_inline)(res, x); }

void DB(set_dd)(DB(t) res, double x0, double x1)
{
    double t[2] = { x0, x1 };
    double err = 0.0;
    _dfloat_renorm(res->d, N, &err, t, 2);
    res->rad = _dfloat_rad_finish(err);
    DB(_fix_nonfinite)(res->d, &res->rad);
}

void DB(set_d_rad)(DB(t) res, double x, double rad)
{
    DB(_set_d_inline)(res, x);
    /* (a nonfinite x or rad, or rad < 0, gives the whole line) */
    res->rad = (rad >= 0.0 && fabs(x) <= DBL_MAX) ? rad : D_INF;
    DB(_fix_nonfinite)(res->d, &res->rad);
}

FLINT_FORCE_INLINE void DB(_zero_inline)(DB(t) res) { DB(_set_d_inline)(res, 0.0); }
void DB(zero)(DB(t) res) { DB(_zero_inline)(res); }
FLINT_FORCE_INLINE void DB(_one_inline)(DB(t) res) { DB(_set_d_inline)(res, 1.0); }
void DB(one)(DB(t) res) { DB(_one_inline)(res); }

FLINT_FORCE_INLINE void DB(_set_inline)(DB(t) res, DB(srcptr) x) { *res = *x; }
void DB(set)(DB(t) res, DB(srcptr) x) { DB(_set_inline)(res, x); }

/* [0 +/- b] */
FLINT_FORCE_INLINE void
DB(_set_zero_pm)(DB(t) res, double b)
{
    int i;
    for (i = 0; i < N; i++)
        res->d[i] = 0.0;
    res->rad = b;
}

FLINT_FORCE_INLINE void
DB(_neg_inline)(DB(t) res, DB(srcptr) x)
{
    int i;
    for (i = 0; i < N; i++)
        res->d[i] = -x->d[i];
    res->rad = x->rad;
    if (FLINT_UNLIKELY(DB(_is_W)(res)))
        DB(_indeterminate_inline)(res);
}
void DB(neg)(DB(t) res, DB(srcptr) x) { DB(_neg_inline)(res, x); }

/* | [m +/- r] | is contained in [|m| +/- r] */
FLINT_FORCE_INLINE void
DB(_abs_inline)(DB(t) res, DB(srcptr) x)
{
    if (x->d[0] < 0.0)
        DB(_neg_inline)(res, x);
    else if (FLINT_UNLIKELY(DB(_is_W)(x)))
        DB(_indeterminate_inline)(res);
    else
        DB(_set_inline)(res, x);
}
void DB(abs)(DB(t) res, DB(srcptr) x) { DB(_abs_inline)(res, x); }

FLINT_FORCE_INLINE void
DB(_canonicalise_inline)(DB(t) res, DB(srcptr) x)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x)))
    {
        DB(_indeterminate_inline)(res);
        return;
    }
    DB(_set_inline)(res, x);
    DF(_canonicalise_inplace)(res->d);
}
void DB(canonicalise)(DB(t) res, DB(srcptr) x) { DB(_canonicalise_inline)(res, x); }

DFLOAT_ARITH_IMPL void
DB(_add_impl)(DB(t) res, DB(srcptr) x, DB(srcptr) y)
{
    double err = 0.0;
    DF(_add_kernel)(res->d, &err, x->d, y->d, 1);
    res->rad = _dfloat_rad_finish((x->rad + y->rad) + err);
    DB(_fix_nonfinite)(res->d, &res->rad);
}
void DB(add)(DB(t) res, DB(srcptr) x, DB(srcptr) y) { DB(_add_impl)(res, x, y); }

DFLOAT_ARITH_IMPL void
DB(_sub_impl)(DB(t) res, DB(srcptr) x, DB(srcptr) y)
{
    double err = 0.0;
    DF(_sub_kernel)(res->d, &err, x->d, y->d, 1);
    res->rad = _dfloat_rad_finish((x->rad + y->rad) + err);
    DB(_fix_nonfinite)(res->d, &res->rad);
}
void DB(sub)(DB(t) res, DB(srcptr) x, DB(srcptr) y) { DB(_sub_impl)(res, x, y); }

/* |x' y' - x y| <= |x| ry + |y'| rx <= |x| ry + |y| rx + rx ry */
/* Error modes of the products: the rounding error of the midpoint
   product is tracked exactly (track = 1) when both inputs are exact,
   so that exact inputs give exact outputs, and otherwise bounded
   (track = 2, cheaper); DFLOAT_MODE_FAST always bounds.  The mode is
   a compile-time constant in every caller. */

FLINT_FORCE_INLINE void
DB(_mul_m)(DB(t) res, DB(srcptr) x, DB(srcptr) y, int mode)
{
    double err = 0.0, ax, ay, rad, rx = x->rad, ry = y->rad;

    if (FLINT_UNLIKELY(DF(_has_tiny)(x->d) | DF(_has_tiny)(y->d)))
    {
        _dfloat_ball_mul_scaled(res->d, &res->rad, N, x->d, x->rad, y->d, y->rad);
        return;
    }

    ax = _dfloat_abs_sum(x->d, N);
    ay = _dfloat_abs_sum(y->d, N);
    rad = _dfloat_rad_mul(ax, ry) + _dfloat_rad_mul(ay, rx);
    rad += _dfloat_rad_mul(rx, ry);
    if (TRACK_MODE(mode, rad == 0.0) == 1)
        DF(_mul_kernel)(res->d, &err, x->d, y->d, 1);
    else
        DF(_mul_kernel)(res->d, &err, x->d, y->d, 2);
    res->rad = _dfloat_rad_finish(rad + err);
    if (FLINT_UNLIKELY(!(_dfloat_abs_sum(res->d, N) + res->rad <= DBL_MAX)))
    {
        /* the whole line, except that an exact zero (ax = 0 only for
           zero components) times anything is zero */
        if ((ax == 0.0 && rx == 0.0) || (ay == 0.0 && ry == 0.0))
            DB(_zero_inline)(res);
        else
            DB(_indeterminate_inline)(res);
    }
}

FLINT_FORCE_INLINE void
DB(_sqr_m)(DB(t) res, DB(srcptr) x, int mode)
{
    double err = 0.0, ax, rad;

    if (FLINT_UNLIKELY(DF(_has_tiny)(x->d)))
    {
        _dfloat_ball_mul_scaled(res->d, &res->rad, N, x->d, x->rad, x->d, x->rad);
        return;
    }

    ax = _dfloat_abs_sum(x->d, N);
    rad = _dfloat_rad_mul(ax + ax, x->rad) + _dfloat_rad_mul(x->rad, x->rad);
    if (TRACK_MODE(mode, rad == 0.0) == 1)
        DF(_sqr_kernel)(res->d, &err, x->d, 1);
    else
        DF(_sqr_kernel)(res->d, &err, x->d, 2);
    res->rad = _dfloat_rad_finish(rad + err);
    DB(_fix_nonfinite)(res->d, &res->rad);
}

/* out-of-line versions with a run-time mode, for the scalar fallbacks
   of the vector functions, the complex products and, for N >= 2, all
   the products outside the dot products (in either mode: the automatic
   mode's body contains both kernels, so one copy serves both) */
static DFLOAT_NOINLINE void DB(_mul_ool)(DB(t) res, DB(srcptr) x, DB(srcptr) y, int mode) { DB(_mul_m)(res, x, y, mode); }
static DFLOAT_NOINLINE void DB(_sqr_ool)(DB(t) res, DB(srcptr) x, int mode) { DB(_sqr_m)(res, x, mode); }

/* the products in a given mode: inlined for N = 1, where they are
   small, and otherwise the out-of-line versions */
FLINT_FORCE_INLINE void
DB(_mul_mode)(DB(t) res, DB(srcptr) x, DB(srcptr) y, int mode)
{
    if (N == 1)
        DB(_mul_m)(res, x, y, mode);
    else
        DB(_mul_ool)(res, x, y, mode);
}

FLINT_FORCE_INLINE void
DB(_sqr_mode)(DB(t) res, DB(srcptr) x, int mode)
{
    if (N == 1)
        DB(_sqr_m)(res, x, mode);
    else
        DB(_sqr_ool)(res, x, mode);
}

FLINT_FORCE_INLINE void
DB(_mul_impl)(DB(t) res, DB(srcptr) x, DB(srcptr) y) { DB(_mul_mode)(res, x, y, DFLOAT_MODE_AUTO); }
void DB(mul)(DB(t) res, DB(srcptr) x, DB(srcptr) y) { DB(_mul_impl)(res, x, y); }
FLINT_FORCE_INLINE void
DB(_sqr_impl)(DB(t) res, DB(srcptr) x) { DB(_sqr_mode)(res, x, DFLOAT_MODE_AUTO); }
void DB(sqr)(DB(t) res, DB(srcptr) x) { DB(_sqr_impl)(res, x); }

void DB(mul_d)(DB(t) res, DB(srcptr) x, double y)
{
    double err = 0.0, rad;

    if (FLINT_UNLIKELY(DF(_has_tiny)(x->d) || (fabs(y) < TINY_THRESHOLD && y != 0.0)
        || !(fabs(y) <= DBL_MAX)))
    {
        double yy[N];
        int i;
        yy[0] = (fabs(y) <= DBL_MAX) ? y : 0.0;
        for (i = 1; i < N; i++)
            yy[i] = 0.0;
        _dfloat_ball_mul_scaled(res->d, &res->rad, N, x->d, x->rad, yy, (fabs(y) <= DBL_MAX) ? 0.0 : D_INF);
        return;
    }

    rad = _dfloat_rad_mul(fabs(y), x->rad);
    DF(_mul_d_kernel)(res->d, &err, x->d, y, 1);
    res->rad = _dfloat_rad_finish(rad + err);
    if (FLINT_UNLIKELY(!(_dfloat_abs_sum(res->d, N) + res->rad <= DBL_MAX)))
    {
        if (y == 0.0)
            DB(_zero_inline)(res);      /* anything times an exact zero */
        else
            DB(_indeterminate_inline)(res);
    }
}

DFLOAT_ARITH_IMPL void
DB(_mul_2exp_si_impl)(DB(t) res, DB(srcptr) x, slong e)
{
    double err = 0.0, r[N], rad;
    int i;

    if (e >= -1000 && e <= 1000)
    {
        double s = d_mul_2exp(1.0, e);
        for (i = 0; i < N; i++)
            r[i] = x->d[i] * s;
        rad = _dfloat_rad_mul(x->rad, s);
    }
    else
    {
        int ee = (int) FLINT_MAX(FLINT_MIN(e, 5000), -5000);
        for (i = 0; i < N; i++)
            r[i] = ldexp(x->d[i], ee);
        rad = ldexp(x->rad, ee);
        if (x->rad != 0.0 && rad < DFLOAT_RAD_UNDERFLOW_THRESHOLD)
            rad += DFLOAT_RAD_UNDERFLOW_BUMP;
    }

    /* components scaled into the subnormal range may have been rounded
       (possibly up to DBL_MIN itself) */
    for (i = 0; i < N; i++)
        if (x->d[i] != 0.0 && fabs(r[i]) <= DBL_MIN)
            err += DFLOAT_RAD_UNDERFLOW_BUMP;

    for (i = 0; i < N; i++)
        res->d[i] = r[i];
    res->rad = _dfloat_rad_finish(rad + err);
    DB(_fix_nonfinite)(res->d, &res->rad);
}
void DB(mul_2exp_si)(DB(t) res, DB(srcptr) x, slong e) { DB(_mul_2exp_si_impl)(res, x, e); }

void DB(add_error_d)(DB(t) res, double err)
{
    res->rad = _dfloat_rad_finish(res->rad + fabs(err));
    DB(_fix_nonfinite)(res->d, &res->rad);
}

/* the midpoint (zero for the whole line) */
void DB(get_mid_d)(double * mid, DB(srcptr) x)
{
    int i, W = DB(_is_W)(x);
    for (i = 0; i < N; i++)
        mid[i] = W ? 0.0 : x->d[i];
}

int DB(is_exact)(DB(srcptr) x) { return x->rad == 0.0 && !DB(_is_W)(x); }
int DB(is_finite)(DB(srcptr) x) { return !DB(_is_W)(x); }
/* (the sum of the absolute values is 0 only for finite zeros) */
FLINT_FORCE_INLINE int DB(_is_zero_inline)(DB(srcptr) x) { return x->rad == 0.0 && _dfloat_abs_sum(x->d, N) == 0.0; }
int DB(is_zero)(DB(srcptr) x) { return DB(_is_zero_inline)(x); }

/* lower bound on |mid| */
FLINT_FORCE_INLINE double
DB(_abs_mid_lower)(DB(srcptr) x)
{
    return fabs(x->d[0]) - _dfloat_rad_finish(DF(_tail_abs_sum)(x->d));
}

/* the comparisons below are made with a slack factor on the radius,
   so that the rounding of the (non-cancelling) subtraction in the
   lower bound cannot flip a borderline verdict */
/* (every form of the whole line fails the test: a nonfinite tail or
   radius makes the comparison false, and the head is checked) */
FLINT_FORCE_INLINE int
DB(_contains_zero_inline)(DB(srcptr) x)
{
    return !(DB(_abs_mid_lower)(x) > _dfloat_rad_finish(x->rad) && fabs(x->d[0]) <= DBL_MAX);
}
int DB(contains_zero)(DB(srcptr) x) { return DB(_contains_zero_inline)(x); }

/* c in [mid +/- rad]?  (certainly: |mid - c| <= rad up to slack);
   the whole line contains every finite c */
int DB(contains_d)(DB(srcptr) x, double c)
{
    double lo, hi, d0;
    if (!(fabs(c) <= DBL_MAX))
        return 0;
    if (DB(_is_W)(x))
        return 1;
    DF(_diff_bounds_d)(&lo, &hi, &d0, x->d, c);
    return _dfloat_rad_finish(hi) <= x->rad;
}

/* certain comparisons of every point of the ball with c (all 0 for
   the whole line): |mid - c| >= lo > rad (or >= rad) with the sign of
   mid - c, or mid = c exactly with rad = 0 for the non-strict ones */
static int
DB(_cmp_all_d)(DB(srcptr) x, double c, int sign, int strict)
{
    double lo, hi, d0;
    int zero;
    if (DB(_is_W)(x))
        return 0;
    zero = DF(_diff_bounds_d)(&lo, &hi, &d0, x->d, c);
    if (zero)
        return !strict && x->rad == 0.0;
    if ((sign > 0) ? !(d0 > 0.0) : !(d0 < 0.0))
        return 0;
    return strict ? (lo > _dfloat_rad_finish(x->rad)) : (lo >= _dfloat_rad_finish(x->rad));
}
#define DB_ALL_GT(x, c) DB(_cmp_all_d)(x, c, 1, 1)
#define DB_ALL_GE(x, c) DB(_cmp_all_d)(x, c, 1, 0)
#define DB_ALL_LT(x, c) DB(_cmp_all_d)(x, c, -1, 1)
#define DB_ALL_LE(x, c) DB(_cmp_all_d)(x, c, -1, 0)

/* y contained in x: |mx - my| + ry <= rx */
int DB(contains)(DB(srcptr) x, DB(srcptr) y)
{
    double lo, hi, d0;
    if (DB(_is_W)(x))
        return 1;
    if (DB(_is_W)(y))
        return 0;
    DF(_diff_bounds)(&lo, &hi, &d0, x->d, y->d);
    return _dfloat_rad_finish(hi + y->rad) <= x->rad;
}

/* certainly overlapping: |mx - my| <= rx + ry; answers 1 in doubtful
   cases (the negation "certainly disjoint" is what the predicates use) */
FLINT_FORCE_INLINE int
DB(_overlaps_inline)(DB(srcptr) x, DB(srcptr) y)
{
    double lo, hi, d0;
    if (DB(_is_W)(x) || DB(_is_W)(y))
        return 1;
    DF(_diff_bounds)(&lo, &hi, &d0, x->d, y->d);
    return !(lo > _dfloat_rad_finish(x->rad + y->rad));
}
int DB(overlaps)(DB(srcptr) x, DB(srcptr) y) { return DB(_overlaps_inline)(x, y); }

FLINT_FORCE_INLINE int
DB(_equal_inline)(DB(srcptr) x, DB(srcptr) y)
{
    double lo, hi, d0;
    if (x->rad != 0.0 || y->rad != 0.0 || DB(_is_W)(x) || DB(_is_W)(y))
        return 0;
    return DF(_diff_bounds)(&lo, &hi, &d0, x->d, y->d);
}
int DB(equal)(DB(srcptr) x, DB(srcptr) y) { return DB(_equal_inline)(x, y); }

/* -1, 0, 1 with success when decidable */
int DB(cmp)(int * res, DB(srcptr) x, DB(srcptr) y)
{
    double lo, hi, d0;
    int zero;
    if (DB(_is_W)(x) || DB(_is_W)(y))
    {
        *res = 0;
        return GR_UNABLE;
    }
    zero = DF(_diff_bounds)(&lo, &hi, &d0, x->d, y->d);
    if (zero && x->rad == 0.0 && y->rad == 0.0)
    {
        *res = 0;
        return GR_SUCCESS;
    }
    if (lo > _dfloat_rad_finish(x->rad + y->rad))
    {
        *res = (d0 > 0.0) ? 1 : -1;
        return GR_SUCCESS;
    }
    *res = 0;
    return GR_UNABLE;
}

/* the head of the midpoint (0 for the whole line) */
double DB(get_d)(DB(srcptr) x) { return DB(_is_W)(x) ? 0.0 : x->d[0]; }

void DB(get_arb)(arb_t res, DB(srcptr) x)
{
    if (DB(_is_W)(x))
        arb_zero_pm_inf(res);
    else
        _dfloat_get_arb(res, x->d, N, x->rad);
}
void DB(set_arb)(DB(t) res, const arb_t x) { _dfloat_set_arb(res->d, &res->rad, N, x); }

void DB(set_arf)(DB(t) res, const arf_t x)
{
    double err = 0.0;
    if (!arf_is_finite(x))
    {
        DB(_indeterminate_inline)(res);
        return;
    }
    _dfloat_set_arf(res->d, N, &err, x);
    res->rad = _dfloat_rad_finish(err);
}

void DB(set_fmpz)(DB(t) res, const fmpz_t x)
{
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, x);
    DB(set_arf)(res, t);
    arf_clear(t);
}

void DB(set_fmpq)(DB(t) res, const fmpq_t x)
{
    arb_t t;
    arb_init(t);
    arb_set_fmpq(t, x, 53 * N + 64);
    DB(set_arb)(res, t);
    arb_clear(t);
}

int DB(set_str)(DB(t) res, const char * s)
{
    return _dfloat_set_str(res->d, &res->rad, N, s);
}

char * DB(get_str)(DB(srcptr) x, slong digits)
{
    return _dfloat_get_str(x->d, N, x->rad, digits);
}

void DB(print)(DB(srcptr) x)
{
    char * s = DB(get_str)(x, 0);
    flint_printf("%s", s);
    flint_free(s);
}

void DB(randtest)(DB(t) res, flint_rand_t state)
{
    _dfloat_randtest(res->d, N, state);
    res->rad = _dfloat_randtest_rad(state) * fabs(res->d[0]);
    if (res->rad != res->rad)   /* 0 * inf */
        res->rad = D_INF;
    if (n_randint(state, 8) == 0)
        res->rad = _dfloat_randtest_rad(state);
    DB(_fix_nonfinite)(res->d, &res->rad);
}

void DB(randtest_special)(DB(t) res, flint_rand_t state)
{
    _dfloat_randtest_special(res->d, N, state);
    res->rad = _dfloat_randtest_rad(state) * fabs(res->d[0]);
    if (res->rad != res->rad)
        res->rad = D_INF;
    if (n_randint(state, 4) == 0)
        res->rad = _dfloat_randtest_rad(state);
    if (n_randint(state, 4) == 0)
        res->rad = fabs(d_randtest_signed(state, -1074, 1023));
    DB(_fix_nonfinite)(res->d, &res->rad);
}

/* Division and square root.  Fast path: heads within 2^-150 .. 2^150
   and no tiny components, so that no product in the residual
   computation can underflow and no scaling is needed; anything else
   goes through the scaled generic code (scaled.c), which also holds
   the derivation of the bounds. */
#define FAST_RANGE(a) ((a) >= 0x1p-150 && (a) <= 0x1p150)

FLINT_FORCE_INLINE int
DB(_div_m)(DB(t) res, DB(srcptr) x, DB(srcptr) y, int mode)
{
    double q[N], P[N], D[N], Perr = 0.0, Derr = 0.0, E, ylo, aq, rad;
    int i;

    /* (a zero head with a nonzero tail only occurs in the whole line) */
    if (FLINT_UNLIKELY(y->d[0] == 0.0 && y->rad == 0.0))
        return DB(_is_zero_inline)(y) ? GR_DOMAIN : GR_UNABLE;
    /* contains zero: the test of contains_zero without its check of the
       head, which the range check below does (a whole line y that
       passes this test has an infinite head) */
    if (!(DB(_abs_mid_lower)(y) > _dfloat_rad_finish(y->rad)))
        return GR_UNABLE;

    if (FLINT_UNLIKELY(!FAST_RANGE(fabs(x->d[0])) || !FAST_RANGE(fabs(y->d[0]))
        || DF(_has_tiny)(x->d) || DF(_has_tiny)(y->d) || !(x->rad <= DBL_MAX)))
    {
        if (DB(_is_W)(y))
            return GR_UNABLE;
        _dfloat_ball_div(res->d, &res->rad, N, x->d, x->rad, y->d, y->rad);
        return GR_SUCCESS;
    }

    /* |y| - ry > 0 was checked by contains_zero */
    ylo = _dfloat_lower(fabs(y->d[0]) - _dfloat_rad_finish(DF(_tail_abs_sum)(y->d) + y->rad));

    /* (the kernel inlined: as a call, the ball division is 5-15% slower
       for N >= 3) */
    DF(_div_level)(q, x->d, y->d, N, N);
    UNROLL
    for (i = 1; i < N; i++)
        if (fabs(q[i]) < TINY_THRESHOLD)
            q[i] = 0.0;

    /* |x' - q y'| <= |x - q y| + rx + |q| ry, with x - q y in
       D +/- (Derr + Perr) */
    if (TRACK_MODE(mode, x->rad == 0.0 && y->rad == 0.0) == 1)
        DF(_mul_kernel)(P, &Perr, q, y->d, 1);
    else
        DF(_mul_kernel)(P, &Perr, q, y->d, 2);
    DF(_sub_kernel)(D, &Derr, x->d, P, 1);
    E = _dfloat_abs_sum(D, N) + Derr + Perr;
    aq = _dfloat_abs_sum(q, N);
    rad = E + x->rad + _dfloat_rad_mul(aq, y->rad);
    rad = _dfloat_rad_div(rad, ylo);

    UNROLL
    for (i = 0; i < N; i++)
        res->d[i] = q[i];
    res->rad = _dfloat_rad_finish(rad);
    DB(_fix_nonfinite)(res->d, &res->rad);
    return GR_SUCCESS;
}

/* one out-of-line body with a run-time mode for div, the fast
   contexts and the vector fallbacks; inv keeps an inlined copy, which
   the constant numerator simplifies */
static DFLOAT_NOINLINE int DB(_div_ool)(DB(t) res, DB(srcptr) x, DB(srcptr) y, int mode) { return DB(_div_m)(res, x, y, mode); }

int DB(div)(DB(t) res, DB(srcptr) x, DB(srcptr) y) { return DB(_div_ool)(res, x, y, DFLOAT_MODE_AUTO); }

static DFLOAT_NOINLINE int
DB(_inv_m)(DB(t) res, DB(srcptr) x, int mode)
{
    DB(t) one;
    DB(_one_inline)(one);
    return DB(_div_m)(res, one, x, mode);
}

int DB(inv)(DB(t) res, DB(srcptr) x) { return DB(_inv_m)(res, x, DFLOAT_MODE_AUTO); }

static DFLOAT_NOINLINE int
DB(_sqrt_m)(DB(t) res, DB(srcptr) x, int mode)
{
    double s[N], P[N], D[N], Perr = 0.0, Derr = 0.0, E, xlo, slo, rad;
    int i;

    if (FLINT_UNLIKELY(!FAST_RANGE(x->d[0]) || DF(_has_tiny)(x->d) || !(x->rad <= DBL_MAX)))
        return _dfloat_ball_sqrt(res->d, &res->rad, N, x->d, x->rad);

    xlo = x->d[0] - _dfloat_rad_finish(DF(_tail_abs_sum)(x->d) + x->rad);
    if (!(xlo >= 0.0))
        return _dfloat_ball_sqrt(res->d, &res->rad, N, x->d, x->rad);

    DF(_sqrt_kernel)(s, x->d);
    UNROLL
    for (i = 1; i < N; i++)
        if (fabs(s[i]) < TINY_THRESHOLD)
            s[i] = 0.0;

    /* |sqrt(x') - s| <= (|x - s^2| + rx) / s */
    if (TRACK_MODE(mode, x->rad == 0.0) == 1)
        DF(_sqr_kernel)(P, &Perr, s, 1);
    else
        DF(_sqr_kernel)(P, &Perr, s, 2);
    DF(_sub_kernel)(D, &Derr, x->d, P, 1);
    E = _dfloat_abs_sum(D, N) + Derr + Perr;
    slo = _dfloat_lower(fabs(s[0]) - _dfloat_rad_finish(DF(_tail_abs_sum)(s)));
    if (!(slo > 0.0))
        return _dfloat_ball_sqrt(res->d, &res->rad, N, x->d, x->rad);
    rad = _dfloat_rad_div(E + x->rad, slo);

    UNROLL
    for (i = 0; i < N; i++)
        res->d[i] = s[i];
    res->rad = _dfloat_rad_finish(rad);
    DB(_fix_nonfinite)(res->d, &res->rad);
    return GR_SUCCESS;
}

int DB(sqrt)(DB(t) res, DB(srcptr) x) { return DB(_sqrt_m)(res, x, DFLOAT_MODE_AUTO); }

/* For x' in [x +/- rx], all > 0, and r the computed reciprocal root:
   1/sqrt(x') = r (1 - e')^(-1/2) with e' = 1 - x' r^2, so that
   |1/sqrt(x') - r| <= |r| (|e'|/2) (1 + |e'|) for |e'| <= 1/4, and
   |e'| <= |1 - x r^2| + rx r^2. */
FLINT_FORCE_INLINE int
DB(_rsqrt_m)(DB(t) res, DB(srcptr) x, int mode)
{
    double r[N], P[N], Q[N], D[N], one[N], Perr = 0.0, Qerr = 0.0, Derr = 0.0, E, xlo, xhi, rhi, rad;
    int i, track;

    if (FLINT_UNLIKELY(!FAST_RANGE(x->d[0]) || DF(_has_tiny)(x->d) || !(x->rad <= DBL_MAX)))
        return _dfloat_ball_rsqrt(res->d, &res->rad, N, x->d, x->rad);

    xlo = x->d[0] - _dfloat_rad_finish(DF(_tail_abs_sum)(x->d) + x->rad);
    if (!(xlo > 0.0))
    {
        xhi = x->d[0] + _dfloat_rad_finish(DF(_tail_abs_sum)(x->d) + x->rad);
        return (xhi <= 0.0) ? GR_DOMAIN : GR_UNABLE;
    }

    DF(_rsqrt_kernel)(r, x->d);
    UNROLL
    for (i = 1; i < N; i++)
        if (fabs(r[i]) < TINY_THRESHOLD)
            r[i] = 0.0;

    /* (a run-time track in the kernels compiles better here than two
       specialised copies) */
    track = TRACK_MODE(mode, x->rad == 0.0);
    DF(_sqr_kernel)(P, &Perr, r, track);
    DF(_mul_kernel)(Q, &Qerr, x->d, P, track);
    one[0] = 1.0;
    for (i = 1; i < N; i++)
        one[i] = 0.0;
    DF(_sub_kernel)(D, &Derr, one, Q, 1);
    rhi = _dfloat_rad_finish(_dfloat_abs_sum(r, N));
    xhi = x->d[0] + _dfloat_rad_finish(DF(_tail_abs_sum)(x->d) + x->rad);
    E = _dfloat_abs_sum(D, N) + Derr + Qerr + _dfloat_rad_mul(xhi, Perr);
    E = _dfloat_rad_finish(E + _dfloat_rad_mul(x->rad, _dfloat_rad_mul(rhi, rhi)));
    if (!(E <= 0.25))
        return _dfloat_ball_rsqrt(res->d, &res->rad, N, x->d, x->rad);
    rad = _dfloat_rad_mul(rhi, 0.5 * E * (1.0 + E));

    UNROLL
    for (i = 0; i < N; i++)
        res->d[i] = r[i];
    res->rad = _dfloat_rad_finish(rad);
    DB(_fix_nonfinite)(res->d, &res->rad);
    return GR_SUCCESS;
}

/* one out-of-line body with a run-time mode */
static DFLOAT_NOINLINE int DB(_rsqrt_ool)(DB(t) res, DB(srcptr) x, int mode) { return DB(_rsqrt_m)(res, x, mode); }
int DB(rsqrt)(DB(t) res, DB(srcptr) x) { return DB(_rsqrt_ool)(res, x, DFLOAT_MODE_AUTO); }

/* an upper bound on exp(r) - 1 for r >= 0 (r may be +inf) */
static double
DB(_expm1_upper)(double r)
{
    /* for tiny r, r + r^2 rounds to r (and r (1 + 2^-45) may too, if
       r is subnormal): add 2^-1074 */
    if (r < 0x1p-500)
        return (r == 0.0) ? 0.0 : _dfloat_rad_finish(r) + 0x1p-1074;
    if (r <= 1.0)
        return _dfloat_rad_finish(r + r * r);
    if (!(r < 700.0))
        return D_INF;
    return ldexp(1.0, (int) ceil(r * 1.4426950408889634) + 1);
}

/* an upper bound on exp(x') over the ball with the canonical midpoint xc
   and radius rx (+inf where it may overflow): exp of an upper bound hi
   of the ball, for wide balls, where the bound exp(m) (exp(r) - 1) of
   the radius overflows although exp(m + r) does not */
static double
DB(_exp_hi_upper)(const double * xc, double rx)
{
    double hi, t;
    DF(t) h, e;
    /* the tail of xc is below 2^-52 |xc[0]|, and the margin covers it
       and the rounding of the sum */
    hi = xc[0] + rx;
    hi += (fabs(hi) + fabs(xc[0]) + rx) * 0x1p-50;
    if (hi < -745.0)
        return 0x1p-1074;           /* exp(-745) < 2^-1074 */
    if (!(hi <= 709.0))
        return D_INF;
    DF(_set_d_inline)(h, hi);
    DF(exp)(e, h);
    t = _dfloat_rad_finish(_dfloat_abs_sum(e->d, N) * (1.0 + 0x1p-40));
    return t + 0x1p-1069;           /* (a subnormal result) */
}

DFLOAT_FUNC_NOINLINE void DB(exp)(DB(t) res, DB(srcptr) x)
{
    double xc[N], m0, rx, E, ay, rad;
    const double eps = DFLOAT_EXP_EPS(N);
    int idx;

    rx = x->rad;
    if (DB(_is_W)(x))
    {
        DB(_indeterminate_inline)(res);
        return;
    }

    /* the canonical midpoint (exact): its head is the midpoint to
       within a relative 2^-52 */
    UNROLL
    for (idx = 0; idx < N; idx++)
        xc[idx] = x->d[idx];
    DF(_canonicalise_inplace)(xc);
    m0 = xc[0];

    if (m0 > 709.9)
    {
        DB(_indeterminate_inline)(res);
        return;
    }
    if (m0 < -746.0)
    {
        /* exp(m) < 2^-1075: [0 +/- exp(m + r)] */
        DB(_zero_inline)(res);
        res->rad = DB(_exp_hi_upper)(xc, rx);
        return;
    }
    if (m0 == 0.0)
    {
        /* exp(0) = 1 exactly; the radius covers exp([-r, r]) - 1 */
        DF(_one_inline)((DF(ptr)) res->d);
        res->rad = (rx == 0.0) ? 0.0 : DB(_expm1_upper)(rx);
        DB(_fix_nonfinite)(res->d, &res->rad);
        return;
    }

    if (FLINT_UNLIKELY(!DF(_is_canonical)(xc) || DF(_exp_core)(res->d, xc)))
        _dfloat_exp_slow(res->d, N, xc);

    /* exp(m) = y / (1 + delta), |delta| <= eps, and for |x' - m| <= r,
       |exp(x') - exp(m)| <= exp(m) (exp(r) - 1); together
       |exp(x') - y| <= |y| (1 + 2 eps) (E + 2 eps) with E >= exp(r) - 1.
       A subnormal y (or y = 0) can carry a further absolute error of
       2^-1070, which enters the bound for exp(m) in the E term too. */
    E = DB(_expm1_upper)(rx);
    ay = _dfloat_abs_sum(res->d, N) * (1.0 + 2.0 * eps);
    if (FLINT_UNLIKELY(ay < DFLOAT_RAD_UNDERFLOW_THRESHOLD))
        ay += 0x1p-1069;
    rad = _dfloat_rad_finish(_dfloat_rad_mul(ay, E) + ay * (2.0 * eps));
    if (FLINT_UNLIKELY(rad < DFLOAT_RAD_UNDERFLOW_THRESHOLD))
        rad += 0x1p-1069;
    if (FLINT_UNLIKELY(!(rad <= DBL_MAX)))
    {
        /* a wide ball: [0 +/- exp(m + r)] unless that overflows too */
        rad = DB(_exp_hi_upper)(xc, rx);
        if (rad <= DBL_MAX)
            DB(_zero_inline)(res);
    }
    res->rad = rad;
    DB(_fix_nonfinite)(res->d, &res->rad);
}

/* expm1 on balls: the kernel on m, and |expm1(x') - expm1(m)| <=
   exp(m) (exp(r) - 1) with exp(m) = 1 + expm1(m) <= 1 + max(y, 0)
   (1 + 2 eps) */
DFLOAT_FUNC_NOINLINE void DB(expm1)(DB(t) res, DB(srcptr) x)
{
    double xc[N], m0, rx, E, ay, em, rad;
    const double eps = DFLOAT_EXPM1_EPS(N);
    int idx;

    rx = x->rad;
    if (DB(_is_W)(x))
    {
        DB(_indeterminate_inline)(res);
        return;
    }
    UNROLL
    for (idx = 0; idx < N; idx++)
        xc[idx] = x->d[idx];
    DF(_canonicalise_inplace)(xc);
    m0 = xc[0];
    if (m0 > 709.9)
    {
        DB(_indeterminate_inline)(res);
        return;
    }
    if (m0 < -746.0)
    {
        /* -1 + exp(m) with exp(m) < 2^-1075: [-1 +/- exp(m + r)] */
        DF(_set_d_inline)((DF(ptr)) res->d, -1.0);
        res->rad = DB(_exp_hi_upper)(xc, rx);
        DB(_fix_nonfinite)(res->d, &res->rad);
        return;
    }
    if (m0 == 0.0 && _dfloat_abs_sum(xc, N) == 0.0)
    {
        DB(_zero_inline)(res);
        res->rad = (rx == 0.0) ? 0.0 : DB(_expm1_upper)(rx);
        return;
    }
    if (FLINT_UNLIKELY(!DF(_is_canonical)(xc) || DF(_expm1_core)(res->d, xc)))
        _dfloat_expm1_slow(res->d, N, xc);
    E = DB(_expm1_upper)(rx);
    ay = _dfloat_abs_sum(res->d, N);
    em = 1.0 + ((res->d[0] > 0.0) ? ay * (1.0 + 2.0 * eps) : 0.0);
    rad = _dfloat_rad_finish(_dfloat_rad_mul(em, E) + ay * eps);
    if (FLINT_UNLIKELY(rad < DFLOAT_RAD_UNDERFLOW_THRESHOLD))
        rad += 0x1p-1069;
    if (FLINT_UNLIKELY(!(rad <= DBL_MAX)))
    {
        /* a wide ball: [-1 +/- exp(m + r)] unless that overflows too */
        rad = DB(_exp_hi_upper)(xc, rx);
        if (rad <= DBL_MAX)
            DF(_set_d_inline)((DF(ptr)) res->d, -1.0);
    }
    res->rad = rad;
    DB(_fix_nonfinite)(res->d, &res->rad);
}

/* vectors */

void DBV(vec_add)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
        if (FLINT_UNLIKELY(DB(_vaddsub4)(res + i, x + i, y + i, 0)))
        {
            DB(_add_impl)(res + i, x + i, y + i);
            DB(_add_impl)(res + i + 1, x + i + 1, y + i + 1);
            DB(_add_impl)(res + i + 2, x + i + 2, y + i + 2);
            DB(_add_impl)(res + i + 3, x + i + 3, y + i + 3);
        }
#endif
    for (; i < len; i++)
        DB(_add_impl)(res + i, x + i, y + i);
}

void DBV(vec_sub)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
        if (FLINT_UNLIKELY(DB(_vaddsub4)(res + i, x + i, y + i, 1)))
        {
            DB(_sub_impl)(res + i, x + i, y + i);
            DB(_sub_impl)(res + i + 1, x + i + 1, y + i + 1);
            DB(_sub_impl)(res + i + 2, x + i + 2, y + i + 2);
            DB(_sub_impl)(res + i + 3, x + i + 3, y + i + 3);
        }
#endif
    for (; i < len; i++)
        DB(_sub_impl)(res + i, x + i, y + i);
}

/* res[i] = x[i] y[i] (ystep = 1), or x[i] y[0] (ystep = 0, in which
   case y must point to 4 copies of the scalar, for the SIMD loads) */
static DFLOAT_NOINLINE void
DBV(_vec_mul_m)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong ystep, slong len, int mode)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DB(srcptr) yb = y + i * ystep;
        if (FLINT_UNLIKELY(DB(_vmul4)(res + i, x + i, yb, mode)))
        {
            DB(_mul_ool)(res + i, x + i, yb, mode);
            DB(_mul_ool)(res + i + 1, x + i + 1, yb + 1, mode);
            DB(_mul_ool)(res + i + 2, x + i + 2, yb + 2, mode);
            DB(_mul_ool)(res + i + 3, x + i + 3, yb + 3, mode);
        }
    }
#endif
    for (; i < len; i++)
        DB(_mul_ool)(res + i, x + i, y + i * ystep, mode);
}

void DBV(vec_mul)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len) { DBV(_vec_mul_m)(res, x, y, 1, len, DFLOAT_MODE_AUTO); }

/* the products by a scalar go through the SIMD product with the scalar
   broadcast (bitwise the same results as the scalar products) */
static void
DBV(_vec_mul_scalar_m)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c, int mode)
{
    DB(struct) c4[4];
    c4[0] = c4[1] = c4[2] = c4[3] = *c;
    DBV(_vec_mul_m)(res, x, c4, 0, len, mode);
}

#define DFLOAT_ADDMUL_CHUNK 64

static void
DBV(_vec_addmul_scalar_m)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c, int negate, int mode)
{
    DB(struct) c4[4], t[DFLOAT_ADDMUL_CHUNK];
    slong i, m;
    c4[0] = c4[1] = c4[2] = c4[3] = *c;
    for (i = 0; i < len; i += m)
    {
        m = FLINT_MIN(len - i, DFLOAT_ADDMUL_CHUNK);
        DBV(_vec_mul_m)(t, x + i, c4, 0, m, mode);
        if (negate)
            DBV(vec_sub)(res + i, res + i, t, m);
        else
            DBV(vec_add)(res + i, res + i, t, m);
    }
}

void DBV(vec_mul_scalar)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c) { DBV(_vec_mul_scalar_m)(res, x, len, c, DFLOAT_MODE_AUTO); }
void DBV(vec_addmul_scalar)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c) { DBV(_vec_addmul_scalar_m)(res, x, len, c, 0, DFLOAT_MODE_AUTO); }
void DBV(vec_submul_scalar)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c) { DBV(_vec_addmul_scalar_m)(res, x, len, c, 1, DFLOAT_MODE_AUTO); }

FLINT_FORCE_INLINE void
DBV(_vec_dot_m)(DB(t) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len, int rev, int mode)
{
    slong i = 0;
    DB(t) s, t;
    if (initial == NULL)
        DB(_zero_inline)(s);
    else
        DB(_set_inline)(s, initial);
#ifdef DFLOAT_HAVE_SIMD
    if (len >= DFLOAT_DOT_SIMD_MIN)
    {
        /* four lane accumulators kept in registers (SoA), combined at
           the end; a block with a special case is redone in scalar */
        vec4d acc[N], accrad = vec4d_zero(), mx[N], my[N], rx, ry, prod[N], prad, sum[N], srad;
        DB(struct) a4[4];
        int k;
        for (k = 0; k < N; k++)
            acc[k] = vec4d_zero();
        for (; i + 4 <= len; i += 4)
        {
            DF(_vload4)(mx, x[i].d, N + 1);
            DB(_vload_rad)(&rx, x + i);
            if (rev)
            {
                DF(_vload4_rev)(my, y[len - 4 - i].d, N + 1);
                DB(_vload_rad_rev)(&ry, y + len - 4 - i);
            }
            else
            {
                DF(_vload4)(my, y[i].d, N + 1);
                DB(_vload_rad)(&ry, y + i);
            }
            if (FLINT_LIKELY(!DB(_vmul4_soa)(prod, &prad, mx, rx, my, ry, mode)
                && !DB(_vaddsub4_soa)(sum, &srad, acc, accrad, prod, prad, 0)))
            {
                for (k = 0; k < N; k++)
                    acc[k] = sum[k];
                accrad = srad;
            }
            else
            {
                DB(_vstore4)(a4, acc, accrad);
                for (k = 0; k < 4; k++)
                {
                    DB(_mul_m)(t, x + i + k, rev ? y + len - 1 - i - k : y + i + k, mode);
                    DB(_add_impl)(a4 + k, a4 + k, t);
                }
                DF(_vload4)(acc, a4->d, N + 1);
                DB(_vload_rad)(&accrad, a4);
            }
        }
        DB(_vstore4)(a4, acc, accrad);
        DB(_add_impl)(a4, a4, a4 + 1);
        DB(_add_impl)(a4 + 2, a4 + 2, a4 + 3);
        DB(_add_impl)(a4, a4, a4 + 2);
        if (subtract)
            DB(_sub_impl)(s, s, a4);
        else
            DB(_add_impl)(s, s, a4);
    }
#endif
    for (; i < len; i++)
    {
        DB(_mul_m)(t, x + i, rev ? y + len - 1 - i : y + i, mode);
        if (subtract)
            DB(_sub_impl)(s, s, t);
        else
            DB(_add_impl)(s, s, t);
    }
    DB(_set_inline)(res, s);
    if (FLINT_UNLIKELY(DB(_is_W)(res)))    /* (an initial value of any form) */
        DB(_indeterminate_inline)(res);
}

/* short dot products: a lean scalar loop (see DFV(_vec_dot_short)) */
FLINT_FORCE_INLINE void
DBV(_vec_dot_short_m)(DB(t) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong ystep, slong len, int mode)
{
    slong i;
    DB(t) s, t;
    if (initial == NULL)
        DB(_zero_inline)(s);
    else
        DB(_set_inline)(s, initial);
    for (i = 0; i < len; i++)
    {
        DB(_mul_m)(t, x + i, y + i * ystep, mode);
        if (subtract)
            DB(_sub_impl)(s, s, t);
        else
            DB(_add_impl)(s, s, t);
    }
    DB(_set_inline)(res, s);
    if (FLINT_UNLIKELY(DB(_is_W)(res)))    /* (an initial value of any form) */
        DB(_indeterminate_inline)(res);
}

static DFLOAT_NOINLINE int DBV(_vec_dot_short_auto)(DB(t) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong ystep, slong len)
{ DBV(_vec_dot_short_m)(res, initial, subtract, x, y, ystep, len, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static DFLOAT_NOINLINE int DBV(_vec_dot_short_fast)(DB(t) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong ystep, slong len)
{ DBV(_vec_dot_short_m)(res, initial, subtract, x, y, ystep, len, DFLOAT_MODE_FAST); return GR_SUCCESS; }

/* the dot products: one out-of-line body per mode, which must be a
   compile-time constant in the loop (a run-time mode costs up to 30%
   here), with the direction as a run-time argument */
static DFLOAT_NOINLINE int DBV(_vec_dot_auto)(DB(t) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len, int rev)
{ DBV(_vec_dot_m)(res, initial, subtract, x, y, len, rev, DFLOAT_MODE_AUTO); return GR_SUCCESS; }
static DFLOAT_NOINLINE int DBV(_vec_dot_fast)(DB(t) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len, int rev)
{ DBV(_vec_dot_m)(res, initial, subtract, x, y, len, rev, DFLOAT_MODE_FAST); return GR_SUCCESS; }
#define DBV_vec_dot_fwd(res, initial, subtract, x, y, len, mode) \
    (((len) < DFLOAT_DOT_SIMD_MIN) ? \
        (((mode) == DFLOAT_MODE_FAST) ? DBV(_vec_dot_short_fast)(res, initial, subtract, x, y, 1, len) : DBV(_vec_dot_short_auto)(res, initial, subtract, x, y, 1, len)) : \
        (((mode) == DFLOAT_MODE_FAST) ? DBV(_vec_dot_fast)(res, initial, subtract, x, y, len, 0) : DBV(_vec_dot_auto)(res, initial, subtract, x, y, len, 0)))
#define DBV_vec_dot_bwd(res, initial, subtract, x, y, len, mode) \
    (((len) < DFLOAT_DOT_SIMD_MIN) ? \
        (((mode) == DFLOAT_MODE_FAST) ? DBV(_vec_dot_short_fast)(res, initial, subtract, x, y + (len) - 1, -1, len) : DBV(_vec_dot_short_auto)(res, initial, subtract, x, y + (len) - 1, -1, len)) : \
        (((mode) == DFLOAT_MODE_FAST) ? DBV(_vec_dot_fast)(res, initial, subtract, x, y, len, 1) : DBV(_vec_dot_auto)(res, initial, subtract, x, y, len, 1)))

void DBV(vec_dot)(DB(t) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len)
{ (void) DBV_vec_dot_fwd(res, initial, subtract, x, y, len, DFLOAT_MODE_AUTO); }
void DBV(vec_dot_rev)(DB(t) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len)
{ (void) DBV_vec_dot_bwd(res, initial, subtract, x, y, len, DFLOAT_MODE_AUTO); }

/* ---- SIMD division ---- */

#ifdef DFLOAT_HAVE_SIMD

/* top-down TwoSum chain (the fast path of _canonicalise) followed by
   the canonical check; returns the mask of lanes that fail it */
FLINT_FORCE_INLINE vec4d
DF(_v_canonicalise)(vec4d * t)
{
    vec4d bad = vec4d_zero();
    int i;
    if (N > 2)
    {
        UNROLL
        for (i = 0; i < N - 1; i++)
            V_TWO_SUM(t[i], t[i + 1], t[i], t[i + 1]);
    }
    UNROLL
    for (i = 0; i + 1 < N; i++)
        bad = vec4d_bit_or(bad, vec4d_cmp_ne(vec4d_add(t[i], t[i + 1]), t[i]));
    return bad;
}


/* the level-M vector kernels with the output zeroed above M (a zero
   head is left in place: only the value matters here) */
#define VECTOR_LEVEL_OP(fn, M, ...) \
    do { \
        switch (M) \
        { \
            case 1: DFV1(fn)(__VA_ARGS__); break; \
            case 2: DFV2(fn)(__VA_ARGS__); break; \
            case 3: DFV3(fn)(__VA_ARGS__); break; \
            default: DF(_v##fn)(__VA_ARGS__); break; \
        } \
    } while (0)

/* the long division of _div_level on 4 lanes; y[0] must be nonzero */
FLINT_FORCE_INLINE void
DF(_v_div_level)(vec4d * q, const vec4d * x, const vec4d * y, int L, int D)
{
    /* (for N = 1, X[1] appears in code that is dead) */
    vec4d r[N], X[N == 1 ? 2 : N], Z[N], qk[N], y1[N], yinv, p, e, d0;
    int i, k;

    for (i = 0; i < N; i++)
        r[i] = (i < L) ? x[i] : vec4d_zero();
    for (i = 0; i < N; i++)
        y1[i] = (i + 1 < N) ? y[i + 1] : vec4d_zero();
    yinv = vec4d_div(vec4d_set_d(1.0), y[0]);

    UNROLL
    for (k = 0; k < D; k++)
    {
        const int M = L - k;
        qk[k] = (M == 1) ? vec4d_div(r[0], y[0]) : vec4d_mul(r[0], yinv);
        if (k + 1 < D)
        {
            p = vec4d_mul(qk[k], y[0]);
            e = vec4d_fmsub_fused(qk[k], y[0], p);
            d0 = vec4d_sub(r[0], p);
            V_TWO_SUM(X[0], X[1], d0, vec4d_neg(e));
            for (i = 2; i < N; i++)
                X[i] = vec4d_zero();
            for (i = 0; i + 1 < M; i++)
                r[i] = r[i + 1];
            for (i = 0; i < N; i++)
                if (i >= M - 1)
                    r[i] = vec4d_zero();
            VECTOR_LEVEL_OP(_mul_d_kernel, M - 1, Z, NULL, y1, vec4d_neg(qk[k]), 0, 0);
            for (i = 0; i < N; i++)
                if (i >= M - 1)
                    Z[i] = vec4d_zero();
            VECTOR_LEVEL_OP(_add3_kernel, M - 1, r, NULL, X, r, Z, 0, 0);
            for (i = 0; i < N; i++)
                if (i >= M - 1)
                    r[i] = vec4d_zero();
        }
    }
    for (i = 0; i < N; i++)
        if (i >= D)
            qk[i] = vec4d_zero();
    VECTOR_LEVEL_OP(_renorm_orders, D, q, qk, 0);
    for (i = 0; i < N; i++)
        if (i >= D)
            q[i] = vec4d_zero();
}

#endif

void DFV(vec_div)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DF(struct) xs[4], ys[4];
        vec4d mx[N], my[N], mq[N], bad, chk;
        int mask, l;
        DF(_vload4)(mx, x[i].d, N);
        DF(_vload4)(my, y[i].d, N);
        /* zero or nonfinite divisors and nonfinite results: scalar */
        bad = vec4d_cmp_eq(my[0], vec4d_zero());
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(vec4d_abs(my[0]), vec4d_set_d(DBL_MAX)));
        DF(_v_div_level)(mq, mx, my, N, N);
        chk = DF(_vabs_sum)(mq);
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(chk, vec4d_set_d(DBL_MAX)));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
            memcpy(ys, y + i, 4 * sizeof(ys[0]));
        }
        DF(_vstore4)(res[i].d, mq, N);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DF(_div_inline)(res + i + l, xs + l, ys + l);
    }
#endif
    for (; i < len; i++)
        DF(_div_inline)(res + i, x + i, y + i);
}

static DFLOAT_NOINLINE void
DBV(_vec_div_m)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len, int mode)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DB(struct) xs[4], ys[4];
        vec4d mx[N], my[N], mq[N], P[N], D[N], rx, ry, bad, ylo, E, aq, rad, err, ax, ay, rr;
        int mask, l, exact, k;
        DF(_vload4)(mx, x[i].d, N + 1);
        DF(_vload4)(my, y[i].d, N + 1);
        DB(_vload_rad)(&rx, x + i);
        DB(_vload_rad)(&ry, y + i);
        ax = vec4d_abs(mx[0]);
        ay = vec4d_abs(my[0]);
        /* the fast range of DB(div), tiny components, infinite input
           radius, and divisors not bounded away from zero: scalar */
        bad = vec4d_bit_or(vec4d_cmp_nge(ax, vec4d_set_d(0x1p-150)), vec4d_cmp_nle(ax, vec4d_set_d(0x1p150)));
        bad = vec4d_bit_or(bad, vec4d_bit_or(vec4d_cmp_nge(ay, vec4d_set_d(0x1p-150)), vec4d_cmp_nle(ay, vec4d_set_d(0x1p150))));
        bad = vec4d_bit_or(bad, vec4d_bit_or(DF(_vtiny_mask)(mx), DF(_vtiny_mask)(my)));
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(rx, vec4d_set_d(DBL_MAX)));
        /* ylo = lower(|y0| - finish(tail + ry)) */
        rr = vec4d_zero();
        for (k = 1; k < N; k++)
            rr = vec4d_add(rr, vec4d_abs(my[k]));
        rr = vec4d_mul(vec4d_add(rr, ry), vec4d_set_d(DFLOAT_RAD_SLACK));
        ylo = vec4d_mul(vec4d_sub(ay, rr), vec4d_set_d(0x1.fffffffffff00p-1));
        bad = vec4d_bit_or(bad, vec4d_cmp_nge(ylo, vec4d_set_d(0x1p-1000)));

        DF(_v_div_level)(mq, mx, my, N, N);
        /* |x' - q y'| <= |x - q y| + rx + |q| ry, with x - q y in
           D +/- err */
        err = vec4d_zero();
        if (mode == DFLOAT_MODE_FAST)
            exact = 0;
        else
        {
            exact = vec4d_movemask(vec4d_bit_and(vec4d_cmp_eq(rx, vec4d_zero()), vec4d_cmp_eq(ry, vec4d_zero())));
            /* a block with both exact and inexact lanes: all in scalar */
            if (FLINT_UNLIKELY(exact != 0 && exact != 15))
                bad = vec4d_cmp_eq(vec4d_zero(), vec4d_zero());
        }
        if (exact)
            DF(_v_mul_kernel)(P, &err, mq, my, 1, 0);
        else
            DF(_v_mul_kernel)(P, &err, mq, my, 2, 0);
        DF(_v_addsub_kernel)(D, &err, mx, P, 1, 1, 0);
        E = vec4d_add(DF(_vabs_sum)(D), err);
        aq = DF(_vabs_sum)(mq);
        rad = vec4d_add(vec4d_add(E, rx), vec4d_mul(aq, ry));
        rad = vec4d_mul(vec4d_div(rad, ylo), vec4d_set_d(DFLOAT_RAD_SLACK));
        /* radius underflow (including to exactly zero from nonzero
           ingredients: the plain vector products and quotient here do
           not bump) and nonfinite results: scalar */
        bad = vec4d_bit_or(bad, DB(_v_rad_underflow)(rad, vec4d_bit_or(vec4d_bit_or(V_NONZERO(E), V_NONZERO(rx)), V_NONZERO(ry))));
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(vec4d_add(aq, rad), vec4d_set_d(DBL_MAX)));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
            memcpy(ys, y + i, 4 * sizeof(ys[0]));
        }
        DB(_vstore4)(res + i, mq, rad);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    GR_IGNORE(DB(_div_ool)(res + i + l, xs + l, ys + l, mode));
    }
#endif
    for (; i < len; i++)
        GR_IGNORE(DB(_div_ool)(res + i, x + i, y + i, mode));
}

void DBV(vec_div)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len) { DBV(_vec_div_m)(res, x, y, len, DFLOAT_MODE_AUTO); }

/* ---- SIMD square root and reciprocal square root ---- */

#ifdef DFLOAT_HAVE_SIMD

/* one pass of the zero-head repair on the first L components: a zero
   head followed by a nonzero component is shifted down; returns the
   mask of lanes still left with a zero head (repeated zeros) */
FLINT_FORCE_INLINE vec4d
DF(_v_fix_zero_head)(vec4d * t, int L)
{
    vec4d zero = vec4d_zero(), bad = zero;
    int i;
    for (i = 0; i + 1 < L; i++)
    {
        vec4d z = vec4d_bit_and(vec4d_cmp_eq(t[i], zero), V_NONZERO(t[i + 1]));
        vec4d a = vec4d_blendv(t[i], t[i + 1], z);
        t[i + 1] = vec4d_blendv(t[i + 1], t[i], z);
        t[i] = a;
    }
    for (i = 1; i < L; i++)
        bad = vec4d_bit_or(bad, V_NONZERO(t[i]));
    return vec4d_bit_and(bad, vec4d_cmp_eq(t[0], zero));
}

/* _sqrt_kernel on 4 lanes; x > 0 in every lane.  Returns the mask of
   lanes whose intermediate remainder was left with a zero head (to be
   redone in scalar) */
FLINT_FORCE_INLINE vec4d
DF(_v_sqrt_kernel)(vec4d * s, const vec4d * x)
{
    /* (for N = 1, X[1] appears in code that is dead) */
    vec4d t[N], u[N], X[N == 1 ? 2 : N], A[N], B[N], P[N], two_s[N], p, e, d0, zero = vec4d_zero(), bad = zero;
    int i, M, Mn;

    s[0] = vec4d_sqrt(x[0]);
    for (i = 1; i < N; i++)
        s[i] = zero;

    UNROLL
    for (M = 1; M < N; M = Mn)
    {
        Mn = FLINT_MIN(2 * M, N);
        p = vec4d_mul(s[0], s[0]);
        e = vec4d_fmsub_fused(s[0], s[0], p);
        d0 = vec4d_sub(x[0], p);
        V_TWO_SUM(X[0], X[1], d0, vec4d_neg(e));
        for (i = 2; i < N; i++)
            X[i] = zero;
        if (M == 1)
        {
            for (i = 0; i < N; i++)
                A[i] = (i + 1 < N && i + 1 < Mn) ? x[i + 1] : zero;
            /* Mn - 1 == 1 */
            t[0] = vec4d_add(X[0], A[0]);
            for (i = 1; i < N; i++)
                t[i] = zero;
            u[0] = vec4d_div(t[0], vec4d_add(s[0], s[0]));
            for (i = 1; i < N; i++)
                u[i] = zero;
        }
        else
        {
            for (i = 0; i < N; i++)
            {
                A[i] = (i + 1 < N && i + 1 < M) ? s[i + 1] : zero;
                B[i] = (i == 0) ? vec4d_add(s[0], s[0]) : ((i < M) ? s[i] : zero);
            }
            VECTOR_LEVEL_OP(_mul_kernel, Mn - 1, P, NULL, A, B, 0, 0);
            for (i = 0; i < N; i++)
            {
                P[i] = (i < Mn - 1) ? vec4d_neg(P[i]) : zero;
                A[i] = (i + 1 < N && i + 1 < Mn) ? x[i + 1] : zero;
            }
            VECTOR_LEVEL_OP(_add3_kernel, Mn - 1, t, NULL, X, A, P, 0, 0);
            for (i = 0; i < N; i++)
                if (i >= Mn - 1)
                    t[i] = zero;
            bad = vec4d_bit_or(bad, DF(_v_fix_zero_head)(t, Mn - 1));
            for (i = 0; i < N; i++)
                two_s[i] = (i < M) ? vec4d_add(s[i], s[i]) : zero;
            DF(_v_div_level)(u, t, two_s, Mn - 1, M);
        }
        VECTOR_LEVEL_OP(_addsub_kernel, Mn, s, NULL, s, u, 0, 0, 0);
        for (i = 0; i < N; i++)
            if (i >= Mn)
                s[i] = zero;
    }
    return bad;
}

/* _rsqrt_kernel on 4 lanes; x > 0 in every lane.  The intermediate
   zero heads do not matter here (every later operation reads all the
   components); the result may come out non-canonical, which the
   callers check. */
FLINT_FORCE_INLINE void
DF(_v_rsqrt_kernel)(vec4d * r, const vec4d * x)
{
    vec4d P[N], Q[N], X[N], Y[N], t[N], U[N], zero = vec4d_zero();
    int i, M, Mn;

    r[0] = vec4d_div(vec4d_set_d(1.0), vec4d_sqrt(x[0]));
    for (i = 1; i < N; i++)
        r[i] = zero;
    if (N == 1)
        return;

    {
        vec4d p, e, q, f, t0, u0;
        p = vec4d_mul(r[0], r[0]);
        e = vec4d_fmsub_fused(r[0], r[0], p);
        q = vec4d_mul(x[0], p);
        f = vec4d_fmsub_fused(x[0], p, q);
        t0 = vec4d_sub(vec4d_sub(vec4d_set_d(1.0), q), f);
        t0 = vec4d_fnmadd_fused(x[0], e, t0);
        t0 = vec4d_fnmadd_fused(x[1], vec4d_add(p, e), t0);
        u0 = vec4d_mul(vec4d_mul(vec4d_set_d(0.5), r[0]), t0);
        V_TWO_SUM(r[0], r[1], r[0], u0);
    }

    UNROLL
    for (M = 2; M < N; M = Mn)
    {
        Mn = FLINT_MIN(2 * M, N);
        VECTOR_LEVEL_OP(_sqr_kernel, Mn, P, NULL, r, 0, 0);
        for (i = 0; i < N; i++)
            if (i >= Mn)
                P[i] = zero;
        VECTOR_LEVEL_OP(_mul_kernel, Mn, Q, NULL, x, P, 0, 0);
        X[0] = vec4d_sub(vec4d_set_d(1.0), Q[0]);
        for (i = 1; i < N; i++)
            X[i] = zero;
        for (i = 0; i < N; i++)
            Y[i] = (i + 1 < N && i + 1 < Mn) ? Q[i + 1] : zero;
        VECTOR_LEVEL_OP(_addsub_kernel, Mn - 1, t, NULL, X, Y, 1, 0, 0);
        for (i = 0; i < N; i++)
            if (i >= Mn - 1)
                t[i] = zero;
        VECTOR_LEVEL_OP(_mul_kernel, Mn - 1, U, NULL, r, t, 0, 0);
        for (i = 0; i < N; i++)
            U[i] = (i < Mn - 1) ? vec4d_mul(U[i], vec4d_set_d(0.5)) : zero;
        VECTOR_LEVEL_OP(_addsub_kernel, Mn, r, NULL, r, U, 0, 0, 0);
        for (i = 0; i < N; i++)
            if (i >= Mn)
                r[i] = zero;
    }
}

#endif

void DFV(vec_sqrt)(DF(ptr) res, DF(srcptr) x, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DF(struct) xs[4];
        vec4d mx[N], ms[N], bad, chk;
        int mask, l;
        DF(_vload4)(mx, x[i].d, N);
        /* nonpositive or nonfinite arguments: scalar */
        bad = vec4d_cmp_ngt(mx[0], vec4d_zero());
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(mx[0], vec4d_set_d(DBL_MAX)));
        bad = vec4d_bit_or(bad, DF(_v_sqrt_kernel)(ms, mx));
        chk = DF(_vabs_sum)(ms);
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(chk, vec4d_set_d(DBL_MAX)));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DF(_vstore4)(res[i].d, ms, N);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DF(_sqrt_inline)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DF(_sqrt_inline)(res + i, x + i);
}

void DFV(vec_rsqrt)(DF(ptr) res, DF(srcptr) x, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DF(struct) xs[4];
        vec4d mx[N], mr[N], bad, chk;
        int mask, l;
        DF(_vload4)(mx, x[i].d, N);
        bad = vec4d_cmp_ngt(mx[0], vec4d_zero());
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(mx[0], vec4d_set_d(DBL_MAX)));
        DF(_v_rsqrt_kernel)(mr, mx);
        if (N > 1)
            bad = vec4d_bit_or(bad, DF(_v_canonicalise)(mr));
        chk = DF(_vabs_sum)(mr);
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(chk, vec4d_set_d(DBL_MAX)));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DF(_vstore4)(res[i].d, mr, N);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DF(_rsqrt_inline)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DF(_rsqrt_inline)(res + i, x + i);
}

/* the ball versions: the fast range of the scalar functions, a
   posteriori bounds as there (the residuals with exactly tracked
   errors when the block is exact, bounded otherwise), the rest in
   scalar; the status is that of the scalar calls */
FLINT_FORCE_INLINE int
DBV(_vec_sqrt_m)(DB(ptr) res, DB(srcptr) x, slong len, int mode)
{
    slong i = 0;
    int status = GR_SUCCESS;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DB(struct) xs[4];
        vec4d mx[N], ms[N], P[N], D[N], rx, bad, xlo, slo, E, rad, err, rr, ubad;
        int mask, l, exact, k;
        DF(_vload4)(mx, x[i].d, N + 1);
        DB(_vload_rad)(&rx, x + i);
        bad = vec4d_bit_or(vec4d_cmp_nge(mx[0], vec4d_set_d(0x1p-150)), vec4d_cmp_nle(mx[0], vec4d_set_d(0x1p150)));
        bad = vec4d_bit_or(bad, DF(_vtiny_mask)(mx));
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(rx, vec4d_set_d(DBL_MAX)));
        /* xlo = x0 - finish(tail + rx) must be >= 0 */
        rr = vec4d_zero();
        for (k = 1; k < N; k++)
            rr = vec4d_add(rr, vec4d_abs(mx[k]));
        rr = vec4d_mul(vec4d_add(rr, rx), vec4d_set_d(DFLOAT_RAD_SLACK));
        xlo = vec4d_sub(mx[0], rr);
        bad = vec4d_bit_or(bad, vec4d_cmp_nge(xlo, vec4d_zero()));

        bad = vec4d_bit_or(bad, DF(_v_sqrt_kernel)(ms, mx));
        for (k = 1; k < N; k++)
            ms[k] = vec4d_blendv(ms[k], vec4d_zero(), vec4d_cmp_lt(vec4d_abs(ms[k]), vec4d_set_d(TINY_THRESHOLD)));
        /* |sqrt(x') - s| <= (|x - s^2| + rx) / s */
        err = vec4d_zero();
        if (mode == DFLOAT_MODE_FAST)
            exact = 0;
        else
        {
            exact = vec4d_movemask(vec4d_cmp_eq(rx, vec4d_zero()));
            if (FLINT_UNLIKELY(exact != 0 && exact != 15))
                bad = vec4d_cmp_eq(vec4d_zero(), vec4d_zero());
        }
        if (exact)
            DF(_v_sqr_kernel)(P, &err, ms, 1, 0);
        else
            DF(_v_sqr_kernel)(P, &err, ms, 2, 0);
        DF(_v_addsub_kernel)(D, &err, mx, P, 1, 1, 0);
        E = vec4d_add(DF(_vabs_sum)(D), err);
        rr = vec4d_zero();
        for (k = 1; k < N; k++)
            rr = vec4d_add(rr, vec4d_abs(ms[k]));
        slo = vec4d_mul(vec4d_sub(vec4d_abs(ms[0]), vec4d_mul(rr, vec4d_set_d(DFLOAT_RAD_SLACK))), vec4d_set_d(0x1.fffffffffff00p-1));
        bad = vec4d_bit_or(bad, vec4d_cmp_nge(slo, vec4d_set_d(0x1p-1000)));
        rad = vec4d_mul(vec4d_div(vec4d_add(E, rx), slo), vec4d_set_d(DFLOAT_RAD_SLACK));
        ubad = DB(_v_rad_underflow)(rad, vec4d_bit_or(V_NONZERO(E), V_NONZERO(rx)));
        bad = vec4d_bit_or(bad, ubad);
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(vec4d_add(DF(_vabs_sum)(ms), rad), vec4d_set_d(DBL_MAX)));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DB(_vstore4)(res + i, ms, rad);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    status |= DB(_sqrt_m)(res + i + l, xs + l, mode);
    }
#endif
    for (; i < len; i++)
        status |= DB(_sqrt_m)(res + i, x + i, mode);
    return status;
}

/* one out-of-line body with a run-time mode, except for N = 1, where
   specialising the mode is worth 10% and the code is small */
static DFLOAT_NOINLINE int DBV(_vec_sqrt_ool)(DB(ptr) res, DB(srcptr) x, slong len, int mode) { return DBV(_vec_sqrt_m)(res, x, len, mode); }
#define DBV_vec_sqrt_mode(res, x, len, mode) ((N == 1) ? DBV(_vec_sqrt_m)(res, x, len, mode) : DBV(_vec_sqrt_ool)(res, x, len, mode))

int DBV(vec_sqrt)(DB(ptr) res, DB(srcptr) x, slong len) { return DBV_vec_sqrt_mode(res, x, len, DFLOAT_MODE_AUTO); }

static DFLOAT_NOINLINE int
DBV(_vec_rsqrt_m)(DB(ptr) res, DB(srcptr) x, slong len, int mode)
{
    slong i = 0;
    int status = GR_SUCCESS;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DB(struct) xs[4];
        vec4d mx[N], mr[N], P[N], Q[N], D[N], one[N], rx, bad, xlo, xhi, rhi, E, rad, perr, qerr, derr, rr, ubad, zero = vec4d_zero();
        int mask, l, exact, k;
        DF(_vload4)(mx, x[i].d, N + 1);
        DB(_vload_rad)(&rx, x + i);
        bad = vec4d_bit_or(vec4d_cmp_nge(mx[0], vec4d_set_d(0x1p-150)), vec4d_cmp_nle(mx[0], vec4d_set_d(0x1p150)));
        bad = vec4d_bit_or(bad, DF(_vtiny_mask)(mx));
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(rx, vec4d_set_d(DBL_MAX)));
        rr = zero;
        for (k = 1; k < N; k++)
            rr = vec4d_add(rr, vec4d_abs(mx[k]));
        rr = vec4d_mul(vec4d_add(rr, rx), vec4d_set_d(DFLOAT_RAD_SLACK));
        xlo = vec4d_sub(mx[0], rr);
        xhi = vec4d_add(mx[0], rr);
        bad = vec4d_bit_or(bad, vec4d_cmp_ngt(xlo, zero));

        DF(_v_rsqrt_kernel)(mr, mx);
        if (N > 1)
            bad = vec4d_bit_or(bad, DF(_v_canonicalise)(mr));
        for (k = 1; k < N; k++)
            mr[k] = vec4d_blendv(mr[k], zero, vec4d_cmp_lt(vec4d_abs(mr[k]), vec4d_set_d(TINY_THRESHOLD)));
        /* E = |1 - x r^2| + rx r^2, rad = r (E/2)(1 + E) for E <= 1/4 */
        perr = qerr = derr = zero;
        if (mode == DFLOAT_MODE_FAST)
            exact = 0;
        else
        {
            exact = vec4d_movemask(vec4d_cmp_eq(rx, zero));
            if (FLINT_UNLIKELY(exact != 0 && exact != 15))
                bad = vec4d_cmp_eq(zero, zero);
        }
        if (exact)
        {
            DF(_v_sqr_kernel)(P, &perr, mr, 1, 0);
            DF(_v_mul_kernel)(Q, &qerr, mx, P, 1, 0);
        }
        else
        {
            DF(_v_sqr_kernel)(P, &perr, mr, 2, 0);
            DF(_v_mul_kernel)(Q, &qerr, mx, P, 2, 0);
        }
        one[0] = vec4d_set_d(1.0);
        for (k = 1; k < N; k++)
            one[k] = zero;
        DF(_v_addsub_kernel)(D, &derr, one, Q, 1, 1, 0);
        rhi = vec4d_mul(DF(_vabs_sum)(mr), vec4d_set_d(DFLOAT_RAD_SLACK));
        E = vec4d_add(vec4d_add(DF(_vabs_sum)(D), derr), vec4d_add(qerr, vec4d_mul(xhi, perr)));
        E = vec4d_add(E, vec4d_mul(rx, vec4d_mul(rhi, rhi)));
        E = vec4d_mul(E, vec4d_set_d(DFLOAT_RAD_SLACK));
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(E, vec4d_set_d(0.25)));
        rad = vec4d_mul(rhi, vec4d_mul(vec4d_mul(E, vec4d_set_d(0.5)), vec4d_add(vec4d_set_d(1.0), E)));
        rad = vec4d_mul(rad, vec4d_set_d(DFLOAT_RAD_SLACK));
        ubad = DB(_v_rad_underflow)(rad, vec4d_bit_or(V_NONZERO(E), V_NONZERO(rx)));
        bad = vec4d_bit_or(bad, ubad);
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(vec4d_add(DF(_vabs_sum)(mr), rad), vec4d_set_d(DBL_MAX)));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DB(_vstore4)(res + i, mr, rad);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    status |= DB(_rsqrt_ool)(res + i + l, xs + l, mode);
    }
#endif
    for (; i < len; i++)
        status |= DB(_rsqrt_ool)(res + i, x + i, mode);
    return status;
}

int DBV(vec_rsqrt)(DB(ptr) res, DB(srcptr) x, slong len) { return DBV(_vec_rsqrt_m)(res, x, len, DFLOAT_MODE_AUTO); }

/* ---- the exponential function ----

   See exp.c for the algorithm.  The polynomial p(s) = expm1 Taylor
   polynomial of degree EXP_DEG is evaluated through its even and odd
   parts in w = s^2, E(w) = sum_j w^j/(2j+1)! and O(w) = sum_j
   w^j/(2j+2)!, each by Horner's rule; the result is assembled as
   y = T + (T s) E + (T w) O in one three-term sum, with
   T = exp(i/256) exp(j/65536) and the products T s, T w computed
   ahead of the chains (they do not depend on them).  The two chains have the same shape, so the scalar
   kernel runs them in two SIMD lanes when available, which halves
   the latency of the dependent chain of expansion operations.
   DF(_exp_level)[j] is the number of components used at the step
   with coefficients 1/(2j+1)! and 1/(2j+2)! (tapered precision; the
   O chain is padded with leading zero coefficients when it is
   shorter), chosen in dev/dfloat_exp_bound.py.

   The scalar kernel splits further (4-way): with v = w^2,
   E(w) = E0(v) + w E1(v) and O(w) = O0(v) + w O1(v), so that
   y = T + (T s) E0 + (T s w) E1 + (T w) O0 + (T w^2) O1, the four
   chains in the four lanes at EXP_J4 + 1 coefficients each (schedule
   DF(_exp_level4)) and the four products in one vector product; the
   vector kernel (four arguments at a time) keeps the 2-way split. */
#if N == 1
#define EXP_DEG 3
#define EXP_JE 1
#define EXP_JO 0
#define EXP_J4 0
static const int DF(_exp_level)[EXP_JE + 1] = {1, 1};
static const int DF(_exp_level4)[EXP_J4 + 1] = {1};
#elif N == 2
#define EXP_DEG 6
#define EXP_JE 2
#define EXP_JO 2
#define EXP_J4 1
static const int DF(_exp_level)[EXP_JE + 1] = {2, 1, 1};
static const int DF(_exp_level4)[EXP_J4 + 1] = {2, 1};
#elif N == 3
#define EXP_DEG 8
#define EXP_JE 3
#define EXP_JO 3
#define EXP_J4 1
static const int DF(_exp_level)[EXP_JE + 1] = {3, 2, 2, 1};
static const int DF(_exp_level4)[EXP_J4 + 1] = {3, 2};
#else
#define EXP_DEG 11
#define EXP_JE 5
#define EXP_JO 4
#define EXP_J4 2
static const int DF(_exp_level)[EXP_JE + 1] = {4, 3, 3, 2, 1, 1};
static const int DF(_exp_level4)[EXP_J4 + 1] = {4, 3, 1};
#endif

/* the coefficients of step j: E lane 1/(2j+1)!, O lane 1/(2j+2)!
   (zero for the padding steps j > JO), component idx */
#define EXP_CE(j, idx) (_dfloat_exp_invfac[2 * (j) + 1][idx])
#define EXP_CO(j, idx) (((j) <= EXP_JO) ? _dfloat_exp_invfac[2 * (j) + 2][idx] : 0.0)

/* The head of the argument reduction, for a double x0 with
   |x0| <= 746: k = round(x0 / log 2), h = x0 - k L0, i = round(256 h),
   d = h - i/256, j = round(65536 d), s0 = d - j/65536.  All three
   subtractions are exact: k L0 is exact since L0 has 41 bits and
   |k| <= 1078 < 2^12; for |k| >= 2, x0 and k L0 are within a factor 2
   of each other (|x0 - k log 2| <= 0.35 and k log 2 >= 1.38) and
   Sterbenz's lemma applies, for |k| = 1 both are multiples of 2^-54
   with a difference below 2^-1, and for k = 0 there is nothing to
   subtract; d and s0 are differences of two multiples of ulp(h),
   resp. ulp(d), of magnitude at most 1/512, resp. 2^-17, which is at
   most 2^53 ulps when |h| >= 2^-10, resp. |d| >= 2^-18, while i,
   resp. j, is zero otherwise.  (The clamps never act: |h| <= 0.3466 +
   2^-31 gives |i| <= 89, and |d| <= 1/512 gives |j| <= 128.) */
#define EXP_REDUCE_HEAD(x0, k, kd, h, i, d, j, s0) \
    do { \
        (k) = (slong) rint((x0) * 1.4426950408889634); \
        (kd) = (double) (k); \
        (h) = (x0) - (kd) * _dfloat_exp_ln2_42[0]; \
        (i) = (int) rint((h) * 256.0); \
        (i) = FLINT_MAX(FLINT_MIN((i), 89), -89); \
        (d) = (h) - (i) / 256.0; \
        (j) = (int) rint((d) * 65536.0); \
        (j) = FLINT_MAX(FLINT_MIN((j), 128), -128); \
        (s0) = (d) - (j) / 65536.0; \
    } while (0)

#ifdef DFLOAT_HAVE_SIMD
#define V_EXP_REDUCE_HEAD(x0, kd, h, iv, d, jv, s0) \
    do { \
        (kd) = vec4d_round(vec4d_mul((x0), vec4d_set_d(1.4426950408889634))); \
        (h) = vec4d_sub((x0), vec4d_mul((kd), vec4d_set_d(_dfloat_exp_ln2_42[0]))); \
        (iv) = vec4d_round(vec4d_mul((h), vec4d_set_d(256.0))); \
        (iv) = vec4d_max(vec4d_min((iv), vec4d_set_d(89.0)), vec4d_set_d(-89.0)); \
        (d) = vec4d_sub((h), vec4d_mul((iv), vec4d_set_d(1.0 / 256.0))); \
        (jv) = vec4d_round(vec4d_mul((d), vec4d_set_d(65536.0))); \
        (jv) = vec4d_max(vec4d_min((jv), vec4d_set_d(128.0)), vec4d_set_d(-128.0)); \
        (s0) = vec4d_sub((d), vec4d_mul((jv), vec4d_set_d(1.0 / 65536.0))); \
    } while (0)
#endif

/* ---- SIMD exp: the kernel on 4 lanes ---- */

#ifdef DFLOAT_HAVE_SIMD


/* mul / add at level M <= N on vec4d components (M constant after
   unrolling); the higher components of the output are zeroed.  A zero
   head is left in place: the value is exact and the later operations
   do not need the canonical position of the components. */
FLINT_FORCE_INLINE void
DF(_vexp_mul_level)(vec4d * res, const vec4d * x, const vec4d * y, int M)
{
    int i, flag;
    switch (M)
    {
        case 1: flag = DFV1(_mul_kernel)(res, NULL, x, y, 0, 0); break;
        case 2: flag = DFV2(_mul_kernel)(res, NULL, x, y, 0, 0); break;
        case 3: flag = DFV3(_mul_kernel)(res, NULL, x, y, 0, 0); break;
        default: flag = DF(_v_mul_kernel)(res, NULL, x, y, 0, 0); break;
    }
    (void) flag;
    for (i = 0; i < N; i++)
        if (i >= M)
            res[i] = vec4d_zero();
}

FLINT_FORCE_INLINE void
DF(_vexp_add_level)(vec4d * res, const vec4d * x, const vec4d * y, int M)
{
    int i, flag;
    switch (M)
    {
        case 1: flag = DFV1(_addsub_kernel)(res, NULL, x, y, 0, 0, 0); break;
        case 2: flag = DFV2(_addsub_kernel)(res, NULL, x, y, 0, 0, 0); break;
        case 3: flag = DFV3(_addsub_kernel)(res, NULL, x, y, 0, 0, 0); break;
        default: flag = DF(_v_addsub_kernel)(res, NULL, x, y, 0, 0, 0); break;
    }
    (void) flag;
    for (i = 0; i < N; i++)
        if (i >= M)
            res[i] = vec4d_zero();
}

/* Horner's rule in the four lanes of a vec4d: lane l evaluates
   sum_{j=0}^{J} c[l][j] v^j from the top, q = c[J], q <- c[j] + v q,
   at lev[j] components for step j (the mul and the add); v is the
   same in every lane, and the coefficients c[l][j] are the first
   lev[j] components of 5-term expansions (ct[l * (J + 1) + j], NULL
   for a zero coefficient, used for padding a shorter lane).  The
   output components above lev[0] are zero. */
FLINT_FORCE_INLINE void
DF(_horner4)(vec4d * q, const vec4d * v, const double * const * ct, int J, const int * lev)
{
    vec4d c[N], u[N];
    int j, idx, M;
#define EXP_COEF(l, j, idx) ((ct[(l) * (J + 1) + (j)] == NULL) ? 0.0 : ct[(l) * (J + 1) + (j)][idx])
    M = lev[J];
    for (idx = 0; idx < N; idx++)
        q[idx] = (idx < M) ? vec4d_set_d4(EXP_COEF(0, J, idx), EXP_COEF(1, J, idx), EXP_COEF(2, J, idx), EXP_COEF(3, J, idx)) : vec4d_zero();
    for (j = J - 1; j >= 0; j--)
    {
        M = lev[j];
        for (idx = 0; idx < M; idx++)
            c[idx] = vec4d_set_d4(EXP_COEF(0, j, idx), EXP_COEF(1, j, idx), EXP_COEF(2, j, idx), EXP_COEF(3, j, idx));
        DF(_vexp_mul_level)(u, v, q, M);
        DF(_vexp_add_level)(q, c, u, M);
    }
#undef EXP_COEF
}

/* 2^k for integer-valued k in [-1074, 1024] on each lane, as a
   product of two representable powers */
FLINT_FORCE_INLINE void
DF(_v_pow2_split)(vec4d * p1, vec4d * p2, vec4d kd)
{
    vec4d k1 = vec4d_round(vec4d_mul(kd, vec4d_set_d(0.5)));   /* |k1| <= 537 */
    vec4d k2 = vec4d_sub(kd, k1);
    /* the low 12 bits of 2^52 + 2^51 + k + 1023 are k + 1023 (in [0, 2047]) */
    vec4d magic = vec4d_set_d(0x1.8p52 + 1023.0);
    *p1 = vec4d_bits_shl(vec4d_add(k1, magic), 52);
    *p2 = vec4d_bits_shl(vec4d_add(k2, magic), 52);
}

/* gather 4 table rows (5 doubles each), the first N components; row
   holds integers in [0, 2^31 / 5) */
FLINT_FORCE_INLINE void
DF(_v_gather)(vec4d * res, const double * table, vec4d row)
{
    vec4di base = vec4d_to_index(vec4d_mul(row, vec4d_set_d(5.0)));
    int idx;
    UNROLL
    for (idx = 0; idx < N; idx++)
        res[idx] = vec4d_gather(table + idx, base);
}

/* The kernel on 4 lanes: the same computation as _exp_core, on
   canonical inputs xc (SoA); returns the mask of lanes where an
   exactness check failed (to be redone by the scalar code). */
#if N == 1
static vec4d
DF(_v_exp_core)(vec4d * res, const vec4d * xc)
{
    vec4d h, d, kd, iv, bad, p1, e1, s, p, T, y, q1, q2;
    vec4di base;

    kd = vec4d_round(vec4d_mul(xc[0], vec4d_set_d(1.4426950408889634)));
    h = vec4d_sub(xc[0], vec4d_mul(kd, vec4d_set_d(_dfloat_exp_ln2_42[0])));
    iv = vec4d_round(vec4d_mul(h, vec4d_set_d(256.0)));
    iv = vec4d_max(vec4d_min(iv, vec4d_set_d(89.0)), vec4d_set_d(-89.0));
    d = vec4d_sub(h, vec4d_mul(iv, vec4d_set_d(1.0 / 256.0)));
    bad = vec4d_zero();
    p1 = vec4d_mul(kd, vec4d_set_d(_dfloat_exp_ln2_42[1]));
    e1 = vec4d_fmsub_fused(kd, vec4d_set_d(_dfloat_exp_ln2_42[1]), p1);
    s = vec4d_sub(vec4d_sub(d, p1), e1);
    p = vec4d_fmadd_fused(s, vec4d_set_d(1.0 / 720.0), vec4d_set_d(1.0 / 120.0));
    p = vec4d_fmadd_fused(s, p, vec4d_set_d(1.0 / 24.0));
    p = vec4d_fmadd_fused(s, p, vec4d_set_d(1.0 / 6.0));
    p = vec4d_fmadd_fused(s, p, vec4d_set_d(0.5));
    p = vec4d_fmadd_fused(s, p, vec4d_set_d(1.0));
    p = vec4d_mul(s, p);
    base = vec4d_to_index(vec4d_mul(vec4d_add(iv, vec4d_set_d(89.0)), vec4d_set_d(5.0)));
    T = vec4d_gather(&_dfloat_exp_T1[0][0], base);
    y = vec4d_fmadd_fused(T, p, T);
    DF(_v_pow2_split)(&q1, &q2, kd);
    res[0] = vec4d_mul(vec4d_mul(y, q1), q2);
    return bad;
}
#else
static vec4d
DF(_v_exp_core)(vec4d * res, const vec4d * xc)
{
    vec4d s[N], w[N], q[N], r[N], u[N], c[N], T[N], Ts[N], Ts2[N], y[N], Kt[N], Y[N], A[N];
    vec4d h, d, s0, kd, iv, jv, bad, p1, p2;
    int kk, M, idx;

    V_EXP_REDUCE_HEAD(xc[0], kd, h, iv, d, jv, s0);
    bad = vec4d_zero();

    UNROLL
    for (idx = 0; idx < N; idx++)
        c[idx] = vec4d_set_d(_dfloat_exp_ln2_42[idx + 1]);
    DF(_v_mul_d_kernel)(Kt, NULL, c, kd, 0, 0);
    UNROLL
    for (idx = 0; idx + 1 < N; idx++)
        A[idx] = xc[idx + 1];
    A[N - 1] = vec4d_zero();
    DF(_v_addsub_kernel)(Y, NULL, A, Kt, 1, 0, 0);
    A[0] = s0;
    UNROLL
    for (idx = 1; idx < N; idx++)
        A[idx] = vec4d_zero();
    DF(_v_addsub_kernel)(s, NULL, A, Y, 0, 0, 0);

    /* T = exp(i/256) exp(j/65536), T s, T s^2, w = s^2 */
    DF(_v_gather)(c, &_dfloat_exp_T1[0][0], vec4d_add(iv, vec4d_set_d(89.0)));
    DF(_v_gather)(A, &_dfloat_exp_T2[0][0], vec4d_add(jv, vec4d_set_d(128.0)));
    DF(_v_mul_kernel)(T, NULL, c, A, 0, 0);
    DF(_v_sqr_kernel)(w, NULL, s, 0, 0);
    DF(_v_mul_kernel)(Ts, NULL, T, s, 0, 0);
    DF(_v_mul_kernel)(Ts2, NULL, T, w, 0, 0);

    /* the E and O chains */
    M = DF(_exp_level)[EXP_JE];
    for (idx = 0; idx < N; idx++)
    {
        q[idx] = vec4d_set_d((idx < M) ? EXP_CE(EXP_JE, idx) : 0.0);
        r[idx] = vec4d_set_d((idx < M) ? EXP_CO(EXP_JE, idx) : 0.0);
    }
    for (kk = EXP_JE - 1; kk >= 0; kk--)
    {
        M = DF(_exp_level)[kk];
        for (idx = 0; idx < M; idx++)
            c[idx] = vec4d_set_d(EXP_CE(kk, idx));
        DF(_vexp_mul_level)(u, w, q, M);
        DF(_vexp_add_level)(q, c, u, M);
        for (idx = 0; idx < M; idx++)
            c[idx] = vec4d_set_d(EXP_CO(kk, idx));
        DF(_vexp_mul_level)(u, w, r, M);
        DF(_vexp_add_level)(r, c, u, M);
    }

    /* y = T + T s E + T s^2 O */
    DF(_v_mul_kernel)(u, NULL, Ts, q, 0, 0);
    DF(_v_mul_kernel)(y, NULL, Ts2, r, 0, 0);
    DF(_v_add3_kernel)(res, NULL, T, u, y, 0, 0);
    DF(_v_pow2_split)(&p1, &p2, kd);
    UNROLL
    for (idx = 0; idx < N; idx++)
        res[idx] = vec4d_mul(vec4d_mul(res[idx], p1), p2);
    return bad;
}
#endif

#endif

#ifndef DFLOAT_HAVE_SIMD
/* plain (untracked) product and sum at level M <= N, on the first M
   components; the higher components of the output are zeroed */
FLINT_FORCE_INLINE void
DF(_exp_mul_level)(double * res, const double * x, const double * y, int M)
{
    DF(_level_mul)(res, x, y, M);
}

FLINT_FORCE_INLINE void
DF(_exp_add_level)(double * res, const double * x, const double * y, int M)
{
    DF(_level_addsub)(res, x, y, 0, M);
}

#endif

/* The kernel: res = exp(x) (1 + delta) with |delta| <= DFLOAT_EXP_EPS_N
   for a canonical (|x[i+1]| <= u |x[i]|, as checked by _is_canonical)
   x with -746 <= x[0] <= 709.9 (an absolute error of 2^-1070 is added
   to the guarantee when the result is subnormal), returning 0; returns
   nonzero (with res unspecified) in the rare cases the analysis does
   not cover, which the caller handles via arb. */
#if N == 1
/* the double kernel: the reduction with one table level, the Taylor
   polynomial of degree 6 and y = T (1 + p), in plain double arithmetic
   (see dev/dfloat_exp_bound.py, which bounds every rounding) */
FLINT_FORCE_INLINE int
DF(_exp_core_m)(double * res, const double * xc, int m1)
{
    double h, d, kd, p1, e1, s, p, T, y;
    slong k;
    int i;

    /* one table level: s = h - i/256 - k L1 (- k L2, dropped),
       |s| <= 2^-9 + 2^-31, and the Taylor polynomial of degree 6 */
    k = (slong) rint(xc[0] * 1.4426950408889634);
    kd = (double) k;
    h = xc[0] - kd * _dfloat_exp_ln2_42[0];
    i = (int) rint(h * 256.0);
    i = FLINT_MAX(FLINT_MIN(i, 89), -89);
    d = h - i / 256.0;
    DFLOAT_TWO_PROD(p1, e1, kd, _dfloat_exp_ln2_42[1]);
    s = (d - p1) - e1;
    p = fma(s, 1.0 / 720.0, 1.0 / 120.0);
    p = fma(s, p, 1.0 / 24.0);
    p = fma(s, p, 1.0 / 6.0);
    p = fma(s, p, 0.5);
    p = fma(s, p, 1.0);
    p = s * p;
    T = _dfloat_exp_T1[i + 89][0];
    if (m1)
    {
        /* expm1 (k = 0): (T - 1) + T p with T - 1 from the table */
        res[0] = fma(T, p, _dfloat_exp_E1[i + 89][0]);
        return 0;
    }
    y = fma(T, p, T);
    /* y in [0.7, 1.42]: 2^k y by an exponent addition when normal */
    if (FLINT_LIKELY(k >= -1021 && k <= 1022))
    {
        double_uint64_u u;
        u.f = y;
        u.i += ((uint64_t) k) << 52;
        res[0] = u.f;
    }
    else
        res[0] = ldexp(y, (int) k);
    return 0;
}

FLINT_FORCE_INLINE int
DF(_exp_core)(double * res, const double * xc)
{
    return DF(_exp_core_m)(res, xc, 0);
}
#else
/* m1: expm1 instead (for |x0| <= 0.34, so that k = 0): the constant
   term T of the sum is replaced by C = T - 1 = E1 + T1 E2 with E1 =
   exp(i/256) - 1 and E2 = exp(j/65536) - 1 tabulated to relative
   accuracy, which is what keeps expm1 relative for small x (C = 0
   exactly when i = j = 0, and then s = x exactly) */
FLINT_FORCE_INLINE int
DF(_exp_core_m)(double * res, const double * xc, int m1)
{
    double s[N], w[N], c[N], T[N], Ts[N], Ts2[N], y[N], Kt[N], Y[N], A[N];
    double E[N], O[N], v0[N], v1[N], C[N];
    double h, d, s0, kd;
    slong k;
    int i, j, kk, M, idx;

    EXP_REDUCE_HEAD(xc[0], k, kd, h, i, d, j, s0);

    /* T = exp(i/256) exp(j/65536) */
    for (idx = 0; idx < N; idx++)
    {
        c[idx] = _dfloat_exp_T1[i + 89][idx];
        A[idx] = _dfloat_exp_T2[j + 128][idx];
    }
    DF(_mul_impl)((DF(ptr)) T, (DF(srcptr)) c, (DF(srcptr)) A);
    if (m1)
    {
        /* C = E1 + T1 E2 */
        DF(_mul_impl)((DF(ptr)) C, (DF(srcptr)) c, (DF(srcptr)) _dfloat_exp_E2[j + 128]);
        DF(_add_impl)((DF(ptr)) C, (DF(srcptr)) _dfloat_exp_E1[i + 89], (DF(srcptr)) C);
    }

    /* Kt = k (L1 + ... + LN), the tail of k log 2, as N components */
    DF(_mul_d_kernel)(Kt, NULL, _dfloat_exp_ln2_42 + 1, kd, 0);

    /* Y = (x1, ..., x_{N-1}, 0) - Kt; s = s0 + Y */
    UNROLL
    for (idx = 0; idx + 1 < N; idx++)
        A[idx] = xc[idx + 1];
    A[N - 1] = 0.0;
    DF(_sub_kernel)(Y, NULL, A, Kt, 0);
    A[0] = s0;
    UNROLL
    for (idx = 1; idx < N; idx++)
        A[idx] = 0.0;
    DF(_add_kernel)(s, NULL, A, Y, 0);

    /* w = s^2, T s and T w, off the critical path of the chains */
    DF(_sqr_impl)((DF(ptr)) w, (DF(srcptr)) s);
    DF(_mul_impl)((DF(ptr)) Ts, (DF(srcptr)) T, (DF(srcptr)) s);
    DF(_mul_impl)((DF(ptr)) Ts2, (DF(srcptr)) T, (DF(srcptr)) w);

#ifdef DFLOAT_HAVE_SIMD
    {
        /* the four chains E0, E1, O0, O1 in v = w^2 in the four lanes,
           then (T s) E0, (T s w) E1, (T w) O0, (T w^2) O1 in one
           vector product */
        vec4d vv[N], vq[N], vu[N], va[N];
        const double * ct[4 * (EXP_J4 + 1)];
        double v[N], Tsw[N], Tw2[N], v2[N], v3[N];
        double t4[4];

        DF(_sqr_impl)((DF(ptr)) v, (DF(srcptr)) w);
        DF(_mul_impl)((DF(ptr)) Tsw, (DF(srcptr)) Ts, (DF(srcptr)) w);
        DF(_mul_impl)((DF(ptr)) Tw2, (DF(srcptr)) Ts2, (DF(srcptr)) w);
        for (kk = 0; kk <= EXP_J4; kk++)
        {
            /* lane 0: E0 = c_1, c_5, c_9, ...; lane 1: E1 = c_3, c_7, ...;
               lane 2: O0 = c_2, c_6, ...; lane 3: O1 = c_4, c_8, ... */
            ct[0 * (EXP_J4 + 1) + kk] = (4 * kk + 1 <= EXP_DEG) ? _dfloat_exp_invfac[4 * kk + 1] : NULL;
            ct[1 * (EXP_J4 + 1) + kk] = (4 * kk + 3 <= EXP_DEG) ? _dfloat_exp_invfac[4 * kk + 3] : NULL;
            ct[2 * (EXP_J4 + 1) + kk] = (4 * kk + 2 <= EXP_DEG) ? _dfloat_exp_invfac[4 * kk + 2] : NULL;
            ct[3 * (EXP_J4 + 1) + kk] = (4 * kk + 4 <= EXP_DEG) ? _dfloat_exp_invfac[4 * kk + 4] : NULL;
        }
        for (idx = 0; idx < N; idx++)
            vv[idx] = vec4d_set_d(v[idx]);
        DF(_horner4)(vq, vv, ct, EXP_J4, DF(_exp_level4));
        for (idx = 0; idx < N; idx++)
            va[idx] = vec4d_set_d4(Ts[idx], Tsw[idx], Ts2[idx], Tw2[idx]);
        DF(_v_mul_kernel)(vu, NULL, va, vq, 0, 0);
        for (idx = 0; idx < N; idx++)
        {
            vec4d_store_unaligned(t4, vu[idx]);
            v0[idx] = t4[0];
            v1[idx] = t4[1];
            v2[idx] = t4[2];
            v3[idx] = t4[3];
        }
        /* y = T + v0 + v1 + v2 + v3 in one sum */
        {
            const double * xs[5] = {m1 ? C : T, v0, v1, v2, v3};
            if (FLINT_UNLIKELY(DF(_s_addn_kernel)(y, NULL, xs, 5, 0, 0)))
                DF(_shift_zero_head)(y);
        }
        (void) E; (void) O; (void) M;
    }
#else
    {
        /* the chains E and O in w, sequentially (for expm1, whose
           result is s E + w O itself, the chains are not tapered:
           the exp schedule is tuned for errors that s scales down) */
        double u[N];
        M = m1 ? N : DF(_exp_level)[EXP_JE];
        for (idx = 0; idx < N; idx++)
        {
            E[idx] = (idx < M) ? EXP_CE(EXP_JE, idx) : 0.0;
            O[idx] = (idx < M) ? EXP_CO(EXP_JE, idx) : 0.0;
        }
        for (kk = EXP_JE - 1; kk >= 0; kk--)
        {
            M = m1 ? N : DF(_exp_level)[kk];
            for (idx = 0; idx < M; idx++)
                c[idx] = EXP_CE(kk, idx);
            DF(_exp_mul_level)(u, w, E, M);
            DF(_exp_add_level)(E, c, u, M);
            for (idx = 0; idx < M; idx++)
                c[idx] = EXP_CO(kk, idx);
            DF(_exp_mul_level)(u, w, O, M);
            DF(_exp_add_level)(O, c, u, M);
        }
        DF(_mul_impl)((DF(ptr)) v0, (DF(srcptr)) Ts, (DF(srcptr)) E);
        DF(_mul_impl)((DF(ptr)) v1, (DF(srcptr)) Ts2, (DF(srcptr)) O);
    }
#endif

#ifndef DFLOAT_HAVE_SIMD
    /* y = T + T s E + T s^2 O (one three-term sum) */
    if (FLINT_UNLIKELY(DF(_s_add3_kernel)(y, NULL, m1 ? C : T, v0, v1, 0, 0)))
        DF(_shift_zero_head)(y);
#endif
    if (m1)
    {
        for (idx = 0; idx < N; idx++)
            res[idx] = y[idx];
        return 0;
    }
    DF(_mul_2exp_si_impl)((DF(ptr)) res, (DF(srcptr)) y, k);
    return 0;
}
static int
DF(_exp_core)(double * res, const double * xc)
{
    return DF(_exp_core_m)(res, xc, 0);
}
#endif

/* ---- sine and cosine ----

   With k = round(x / (pi/2)), t = x - k pi/2 = A + s where
   A = i/64 + j/2^13 is a table point (|t| <= pi/4 + 2^-13, so
   |i| <= 51 -> clamped to 50, |j| <= 64) and |s| <= 2^-14 + 2^-52,
   and sin t = SA cos s + CA sin s, cos t = CA cos s - SA sin s with
   SA = sin A, CA = cos A from the two table levels (sin(a + b) = sin a
   cos b + cos a sin b, cos(a + b) = cos a cos b - sin a sin b).
   With w = s^2, sin s = s P(w) and cos s = 1 + w Q(w) (Taylor
   polynomials of degrees TRIG_DEG_SIN and TRIG_DEG_COS in s), so that

     sin t = SA + (SA w) Q(w) + (CA s) P(w),
     cos t = CA + (CA w) Q(w) - (SA s) P(w),

   and sin x, cos x are +-sin t, +-cos t by the quadrant k mod 4.  The
   reduction is exact up to the truncation of pi/2 and the dropped
   terms of the expansion arithmetic (an absolute error of about
   u^N |x|), which is what limits the accuracy near multiples of pi:
   the guarantee is |res - sin x| <= EPS (|x| + min(|x|, 1)) and
   |res - cos x| <= EPS (|x| + 1) for |x| <= TRIG_XMAX (dev/dfloat_exp_bound.py,
   which also chooses the schedules).  The scalar kernel splits the
   chains 4-way (P = P0(v) + w P1(v), Q = Q0(v) + w Q1(v), v = w^2),
   the vector kernel 2-way, as for exp. */
#define TRIG_XMAX 1647099.0   /* 2^20 pi/2 - 1: |k| <= 2^20 keeps k P0 exact */
/* below this, sin x = x and cos x = 1 (errors |x|^3/6 <= 2^-600 |x| and
   x^2/2 <= 2^-601), which keeps all intermediate quantities of the
   kernels normal (the bound derivation assumes no underflow) */
#define TRIG_XMIN 0x1p-300
#define TRIG_TINY(sn, cs, xc, want) \
    do { \
        int _i; \
        if ((want) & 1) \
            for (_i = 0; _i < N; _i++) \
                (sn)[_i] = (xc)[_i]; \
        if ((want) & 2) \
        { \
            (cs)[0] = 1.0; \
            for (_i = 1; _i < N; _i++) \
                (cs)[_i] = 0.0; \
        } \
    } while (0)
#if N == 1
#define TRIG_DEG_SIN 7
#define TRIG_DEG_COS 8
#define TRIG_J2 3
#define TRIG_J4 1
static const int DF(_trig_level)[TRIG_J2 + 1] = {1, 1, 1, 1};
static const int DF(_trig_level4)[TRIG_J4 + 1] = {1, 1};
#elif N == 2
#define TRIG_DEG_SIN 7
#define TRIG_DEG_COS 8
#define TRIG_J2 3
#define TRIG_J4 1
static const int DF(_trig_level)[TRIG_J2 + 1] = {2, 2, 1, 1};
static const int DF(_trig_level4)[TRIG_J4 + 1] = {2, 1};
#elif N == 3
#define TRIG_DEG_SIN 11
#define TRIG_DEG_COS 10
#define TRIG_J2 5
#define TRIG_J4 2
static const int DF(_trig_level)[TRIG_J2 + 1] = {3, 3, 2, 1, 1, 1};
static const int DF(_trig_level4)[TRIG_J4 + 1] = {3, 2, 1};
#else
#define TRIG_DEG_SIN 13
#define TRIG_DEG_COS 12
#define TRIG_J2 6
#define TRIG_J4 3
static const int DF(_trig_level)[TRIG_J2 + 1] = {4, 4, 3, 2, 2, 1, 1};
static const int DF(_trig_level4)[TRIG_J4 + 1] = {4, 3, 2, 1};
#endif
/* the coefficient of s^k (a 5-term row), or NULL beyond the degree */
#define TRIG_CS(k) (((k) <= TRIG_DEG_SIN) ? _dfloat_trig_coef[k] : NULL)
#define TRIG_CC(k) (((k) <= TRIG_DEG_COS) ? _dfloat_trig_coef[k] : NULL)

/* the reduction head for a double x0 with |x0| <= TRIG_XMAX: k =
   round(x0 2/pi) as a double kd and its residue q = k mod 4, and h =
   x0 - k P0 exactly (P0 has 33 bits and |k| <= 2^20; for |k| >= 2 the
   two are within a factor 2 of each other, and for |k| = 1 both are
   multiples of 2^-53 with a difference below 1).  The rounding is
   done by adding 1.5 2^52 (|x0 2/pi| < 2^51), which leaves k in the
   low bits of the significand (2^51 + k, so the low two bits are
   those of k also when k is negative). */
#define TRIG_REDUCE_HEAD(x0, q, kd, h) \
    do { \
        double_uint64_u _u; \
        _u.f = (x0) * 0.63661977236758138 + 0x1.8p52; \
        (kd) = _u.f - 0x1.8p52; \
        (q) = (int) (_u.i & 3); \
        (h) = (x0) - (kd) * _dfloat_trig_pi2_33[0]; \
    } while (0)
/* the quadrant: sin x = +-(sin t or cos t) by q, likewise cos x */
#define TRIG_SIN_SIGN(q) (((q) & 2) ? -1.0 : 1.0)
#define TRIG_COS_SIGN(q) ((((q) + 1) & 2) ? -1.0 : 1.0)

/* the table indices and the exact head s0 of t = t0 + tail:
   i = round(64 t0), d = t0 - i/64, j = round(2^13 d), s0 = d - j/2^13
   (exact as for exp: differences of multiples of ulp(t0), ulp(d) of at
   most 2^53 ulps) */
#define TRIG_TABLE_HEAD(t0, i, d, j, s0) \
    do { \
        (i) = (int) rint((t0) * 64.0); \
        (i) = FLINT_MAX(FLINT_MIN((i), 50), -50); \
        (d) = (t0) - (i) / 64.0; \
        (j) = (int) rint((d) * 8192.0); \
        (j) = FLINT_MAX(FLINT_MIN((j), 64), -64); \
        (s0) = (d) - (j) / 8192.0; \
    } while (0)

#if N == 1
/* the double kernel: one table level (i/64, |s| <= 2^-7 + 2^-14, from
   the signed table _dfloat_trig_A1; |t| <= pi/4 + 2^-30 gives |i| <=
   50 without a clamp), the tail k P1 by TwoProd, the polynomials by
   fma */
FLINT_FORCE_INLINE int
DF(_trig_core)(double * sn, double * cs, const double * xc, int want)
{
    double h, kd, p1, e1, t, d, id, s, w, P, Q, sa, ca, st, ct;
    int i, q;

    if (FLINT_UNLIKELY(fabs(xc[0]) < TRIG_XMIN))
    {
        TRIG_TINY(sn, cs, xc, want);
        return 0;
    }
    TRIG_REDUCE_HEAD(xc[0], q, kd, h);
    DFLOAT_TWO_PROD(p1, e1, kd, _dfloat_trig_pi2_33[1]);
    t = h - p1;
    id = rint(t * 64.0);
    i = (int) id;
    d = t - id * (1.0 / 64.0);
    s = d - e1;
    w = s * s;
    /* P = 1 - w/6 + w^2/120 - w^3/5040, Q = -1/2 + w/24 - w^2/720 + w^3/40320 */
    P = fma(w, _dfloat_trig_coef[7][0], _dfloat_trig_coef[5][0]);
    P = fma(w, P, _dfloat_trig_coef[3][0]);
    P = fma(w, P, 1.0);
    Q = fma(w, _dfloat_trig_coef[8][0], _dfloat_trig_coef[6][0]);
    Q = fma(w, Q, _dfloat_trig_coef[4][0]);
    Q = fma(w, Q, _dfloat_trig_coef[2][0]);
    sa = _dfloat_trig_A1[i + 50][0];
    ca = _dfloat_trig_A1[i + 50][1];
    /* sin t = sa + (sa w) Q + (ca s) P, cos t = ca + (ca w) Q - (sa s) P */
    st = fma(sa * w, Q, fma(ca * s, P, sa));
    ct = fma(ca * w, Q, fma(-sa * s, P, ca));
    if (want & 1)
        sn[0] = ((q & 1) ? ct : st) * TRIG_SIN_SIGN(q);
    if (want & 2)
        cs[0] = ((q & 1) ? st : ct) * TRIG_COS_SIGN(q);
    return 0;
}
#else
#ifndef DFLOAT_HAVE_SIMD
/* SA = sin A, CA = cos A for A = i/64 + j/2^13 (N components) */
FLINT_FORCE_INLINE void
DF(_trig_table)(double * SA, double * CA, int i, int j)
{
    double sa[N], ca[N], sb[N], cb[N], p1[N], p2[N];
    double si = (i < 0) ? -1.0 : 1.0, sj = (j < 0) ? -1.0 : 1.0;
    int idx, ai = FLINT_ABS(i), aj = FLINT_ABS(j);
    for (idx = 0; idx < N; idx++)
    {
        sa[idx] = _dfloat_trig_A[ai][0][idx] * si;
        ca[idx] = _dfloat_trig_A[ai][1][idx];
        sb[idx] = _dfloat_trig_B[aj][0][idx] * sj;
        cb[idx] = _dfloat_trig_B[aj][1][idx];
    }
    DF(_mul_impl)((DF(ptr)) p1, (DF(srcptr)) sa, (DF(srcptr)) cb);
    DF(_mul_impl)((DF(ptr)) p2, (DF(srcptr)) ca, (DF(srcptr)) sb);
    DF(_add_impl)((DF(ptr)) SA, (DF(srcptr)) p1, (DF(srcptr)) p2);
    DF(_mul_impl)((DF(ptr)) p1, (DF(srcptr)) ca, (DF(srcptr)) cb);
    DF(_mul_impl)((DF(ptr)) p2, (DF(srcptr)) sa, (DF(srcptr)) sb);
    DF(_sub_impl)((DF(ptr)) CA, (DF(srcptr)) p1, (DF(srcptr)) p2);
}
#endif

/* sn = sin x (want & 1), cs = cos x (want & 2) for a canonical x with
   |x0| <= TRIG_XMAX; returns nonzero if the reduced argument failed
   the canonical check (the caller then uses arb) */
static int
DF(_trig_core)(double * sn, double * cs, const double * xc, int want)
{
    double h, kd, d, s0;
    double Kt[N], A[N], T[N], s[N], w[N], st[N], ct[N];
    int i, j, idx, q, need_st, need_ct;

    if (FLINT_UNLIKELY(fabs(xc[0]) < TRIG_XMIN))
    {
        TRIG_TINY(sn, cs, xc, want);
        return 0;
    }
    TRIG_REDUCE_HEAD(xc[0], q, kd, h);
    DF(_mul_d_kernel)(Kt, NULL, _dfloat_trig_pi2_33 + 1, kd, 0);
    A[0] = h;
    for (idx = 1; idx < N; idx++)
        A[idx] = xc[idx];
    DF(_sub_kernel)(T, NULL, A, Kt, 0);
    /* canonical, so that the tail of t (which becomes the tail of s)
       is bounded by u^i |t| whatever the size of x: without this a
       large x_1 leaves a large second component behind */
    DF(_canonicalise_inplace)(T);
    if (FLINT_UNLIKELY(!DF(_is_canonical)(T)))
        return 1;
    TRIG_TABLE_HEAD(T[0], i, d, j, s0);
    s[0] = s0;
    for (idx = 1; idx < N; idx++)
        s[idx] = T[idx];

    need_st = ((want & 1) && !(q & 1)) || ((want & 2) && (q & 1));
    need_ct = ((want & 2) && !(q & 1)) || ((want & 1) && (q & 1));

    DF(_sqr_impl)((DF(ptr)) w, (DF(srcptr)) s);

#ifdef DFLOAT_HAVE_SIMD
    {
        /* The chains P0, P1, Q0, Q1 in v = w^2 in the four lanes, then
           sin t = SA + (SA w) Q0 + (SA w^2) Q1 + (CA s) P0 + (CA s w) P1
           and cos t = CA + (CA w) Q0 + (CA w^2) Q1 - (SA s) P0 - (SA s w) P1.
           The table products (sa cb, ca sb, ca cb, sa sb), the
           multipliers (SA w, CA s, CA w, SA s) and their products with
           w are each one vector multiplication, so that a pass costs
           one more (the products with the chains) and a five-term sum. */
        vec4d vv[N], vq[N], vu[N], va[N], vb[N], vt[N], vm[N], vmw[N];
        const double * ctab[4 * (TRIG_J4 + 1)];
        double v[N], X[N], l0[N], l1[N], l2[N], l3[N];
        double t4[4];
        int kk, pass;

        DF(_sqr_impl)((DF(ptr)) v, (DF(srcptr)) w);
        for (kk = 0; kk <= TRIG_J4; kk++)
        {
            /* P0: c1, c5, ...; P1: c3, c7, ...; Q0: c2, c6, ...; Q1: c4, c8, ... */
            ctab[0 * (TRIG_J4 + 1) + kk] = TRIG_CS(4 * kk + 1);
            ctab[1 * (TRIG_J4 + 1) + kk] = TRIG_CS(4 * kk + 3);
            ctab[2 * (TRIG_J4 + 1) + kk] = TRIG_CC(4 * kk + 2);
            ctab[3 * (TRIG_J4 + 1) + kk] = TRIG_CC(4 * kk + 4);
        }
        for (idx = 0; idx < N; idx++)
            vv[idx] = vec4d_set_d(v[idx]);
        DF(_horner4)(vq, vv, ctab, TRIG_J4, DF(_trig_level4));

        /* (SA, CA, SA, CA) = (p0 + p1, p2 - p3, ...) from the lanes
           (sa cb, ca sb, ca cb, sa sb) */
        {
            double sa[N], ca[N], sb[N], cb[N];
            double si = (i < 0) ? -1.0 : 1.0, sj = (j < 0) ? -1.0 : 1.0;
            int ai = FLINT_ABS(i), aj = FLINT_ABS(j);
            for (idx = 0; idx < N; idx++)
            {
                sa[idx] = _dfloat_trig_A[ai][0][idx] * si;
                ca[idx] = _dfloat_trig_A[ai][1][idx];
                sb[idx] = _dfloat_trig_B[aj][0][idx] * sj;
                cb[idx] = _dfloat_trig_B[aj][1][idx];
                va[idx] = vec4d_set_d4(sa[idx], ca[idx], ca[idx], sa[idx]);
                vb[idx] = vec4d_set_d4(cb[idx], sb[idx], cb[idx], sb[idx]);
            }
            DF(_v_mul_kernel)(vu, NULL, va, vb, 0, 0);
            for (idx = 0; idx < N; idx++)
            {
                va[idx] = vec4d_permute_0_2_0_2(vu[idx]);   /* p0, p2, p0, p2 */
                vb[idx] = vec4d_bit_xor(vec4d_permute_1_3_1_3(vu[idx]), vec4d_set_d4(0.0, -0.0, 0.0, -0.0));
            }
            DF(_v_addsub_kernel)(vt, NULL, va, vb, 0, 0, 0);
        }
        /* the multipliers m = (SA w, CA s, CA w, SA s) and m w */
        for (idx = 0; idx < N; idx++)
        {
            va[idx] = vec4d_permute_0_1_1_0(vt[idx]);   /* SA, CA, CA, SA */
            vb[idx] = vec4d_set_d4(w[idx], s[idx], w[idx], s[idx]);
        }
        DF(_v_mul_kernel)(vm, NULL, va, vb, 0, 0);
        for (idx = 0; idx < N; idx++)
            vb[idx] = vec4d_set_d(w[idx]);
        DF(_v_mul_kernel)(vmw, NULL, vm, vb, 0, 0);

        for (pass = 0; pass < 2; pass++)
        {
            /* pass 0: sin t (if needed), pass 1: cos t (if needed) */
            double * out = (pass == 0) ? st : ct;
            if (!((pass == 0) ? need_st : need_ct))
                continue;
            for (idx = 0; idx < N; idx++)
            {
                vec4d lo = vec4d_unpacklo(vm[idx], vmw[idx]);   /* m0, m0 w, m2, m2 w */
                vec4d hi = vec4d_unpackhi(vm[idx], vmw[idx]);   /* m1, m1 w, m3, m3 w */
                if (pass == 0)
                    /* CA s, CA s w, SA w, SA w^2 */
                    va[idx] = vec4d_permute2_0_2(hi, lo);
                else
                    /* -SA s, -SA s w, CA w, CA w^2 */
                    va[idx] = vec4d_bit_xor(vec4d_permute2_1_3(hi, lo), vec4d_set_d4(-0.0, -0.0, 0.0, 0.0));
            }
            DF(_v_mul_kernel)(vu, NULL, va, vq, 0, 0);
            for (idx = 0; idx < N; idx++)
            {
                vec4d_store_unaligned(t4, vu[idx]);
                l0[idx] = t4[0]; l1[idx] = t4[1]; l2[idx] = t4[2]; l3[idx] = t4[3];
                vec4d_store_unaligned(t4, vt[idx]);
                X[idx] = t4[pass];
            }
            {
                const double * xs[5] = {X, l0, l1, l2, l3};
                if (FLINT_UNLIKELY(DF(_s_addn_kernel)(out, NULL, xs, 5, 0, 0)))
                    DF(_shift_zero_head)(out);
            }
        }
    }
#else
    {
        /* the chains P and Q in w, sequentially */
        double P[N], Q[N], u[N], c[N], m0[N], m1[N], l0[N], l1[N], SA[N], CA[N];
        int kk, M, pass;
        DF(_trig_table)(SA, CA, i, j);
        M = DF(_trig_level)[TRIG_J2];
        for (idx = 0; idx < N; idx++)
        {
            P[idx] = (idx < M && TRIG_CS(2 * TRIG_J2 + 1) != NULL) ? TRIG_CS(2 * TRIG_J2 + 1)[idx] : 0.0;
            Q[idx] = (idx < M && TRIG_CC(2 * TRIG_J2 + 2) != NULL) ? TRIG_CC(2 * TRIG_J2 + 2)[idx] : 0.0;
        }
        for (kk = TRIG_J2 - 1; kk >= 0; kk--)
        {
            M = DF(_trig_level)[kk];
            for (idx = 0; idx < N; idx++)
                c[idx] = (idx < M && TRIG_CS(2 * kk + 1) != NULL) ? TRIG_CS(2 * kk + 1)[idx] : 0.0;
            DF(_level_mul)(u, w, P, M);
            DF(_level_addsub)(P, c, u, 0, M);
            for (idx = 0; idx < N; idx++)
                c[idx] = (idx < M && TRIG_CC(2 * kk + 2) != NULL) ? TRIG_CC(2 * kk + 2)[idx] : 0.0;
            DF(_level_mul)(u, w, Q, M);
            DF(_level_addsub)(Q, c, u, 0, M);
        }
        for (pass = 0; pass < 2; pass++)
        {
            const double * X = (pass == 0) ? SA : CA;
            const double * Y = (pass == 0) ? CA : SA;
            double * out = (pass == 0) ? st : ct;
            if (!((pass == 0) ? need_st : need_ct))
                continue;
            DF(_mul_impl)((DF(ptr)) m0, (DF(srcptr)) X, (DF(srcptr)) w);
            DF(_mul_impl)((DF(ptr)) m1, (DF(srcptr)) Y, (DF(srcptr)) s);
            if (pass == 1)
                DF(_neg_inline)((DF(ptr)) m1, (DF(srcptr)) m1);
            DF(_mul_impl)((DF(ptr)) l0, (DF(srcptr)) m0, (DF(srcptr)) Q);
            DF(_mul_impl)((DF(ptr)) l1, (DF(srcptr)) m1, (DF(srcptr)) P);
            if (FLINT_UNLIKELY(DF(_s_add3_kernel)(out, NULL, X, l0, l1, 0, 0)))
                DF(_shift_zero_head)(out);
        }
    }
#endif

    /* the quadrant */
    if (want & 1)
    {
        const double * src = (q & 1) ? ct : st;
        double sg = TRIG_SIN_SIGN(q);
        for (idx = 0; idx < N; idx++)
            sn[idx] = src[idx] * sg;
    }
    if (want & 2)
    {
        const double * src = (q & 1) ? st : ct;
        double sg = TRIG_COS_SIGN(q);
        for (idx = 0; idx < N; idx++)
            cs[idx] = src[idx] * sg;
    }
    return 0;
}
#endif

DFLOAT_FUNC_NOINLINE void DF(sin_cos)(DF(t) sn, DF(t) cs, DF(srcptr) x)
{
    double xc[N];
    int i;

    for (i = 0; i < N; i++)
        xc[i] = x->d[i];
    DF(_canonicalise_inplace)(xc);
    if (FLINT_UNLIKELY(!(fabs(xc[0]) <= TRIG_XMAX) || !DF(_is_canonical)(xc)))
    {
        _dfloat_trig_slow(sn->d, cs->d, N, xc, 3);
        return;
    }
    if (FLINT_UNLIKELY(DF(_trig_core)(sn->d, cs->d, xc, 3)))
        _dfloat_trig_slow(sn->d, cs->d, N, xc, 3);
}

DFLOAT_FUNC_NOINLINE void DF(sin)(DF(t) res, DF(srcptr) x)
{
    double xc[N], cs[N];
    int i;

    for (i = 0; i < N; i++)
        xc[i] = x->d[i];
    DF(_canonicalise_inplace)(xc);
    if (FLINT_UNLIKELY(!(fabs(xc[0]) <= TRIG_XMAX) || !DF(_is_canonical)(xc)))
    {
        _dfloat_trig_slow(res->d, cs, N, xc, 1);
        return;
    }
    if (FLINT_UNLIKELY(DF(_trig_core)(res->d, cs, xc, 1)))
        _dfloat_trig_slow(res->d, cs, N, xc, 1);
}

DFLOAT_FUNC_NOINLINE void DF(cos)(DF(t) res, DF(srcptr) x)
{
    double xc[N], sn[N];
    int i;

    for (i = 0; i < N; i++)
        xc[i] = x->d[i];
    DF(_canonicalise_inplace)(xc);
    if (FLINT_UNLIKELY(!(fabs(xc[0]) <= TRIG_XMAX) || !DF(_is_canonical)(xc)))
    {
        _dfloat_trig_slow(sn, res->d, N, xc, 2);
        return;
    }
    if (FLINT_UNLIKELY(DF(_trig_core)(sn, res->d, xc, 2)))
        _dfloat_trig_slow(sn, res->d, N, xc, 2);
}

/* the radius of sin / cos of [m +/- r] from the static bound:
   |res - f(m)| <= EPS (|m| + a) with a = min(|m|, 1) for sin and 1 for
   cos, and |f(x') - f(m)| <= r; a radius of 1 or more is replaced by
   the trivial enclosure [0 +/- 1] */
FLINT_FORCE_INLINE void
DB(_trig_finish)(double * res, double * rad, double am, double a, double r)
{
    double eps = DFLOAT_TRIG_EPS(N);
    double rr = _dfloat_rad_finish(eps * (am + a) + r);
    int i;
    if (FLINT_UNLIKELY(rr < DFLOAT_RAD_UNDERFLOW_THRESHOLD))
        rr += 0x1p-1069;
    if (!(rr < 1.0))
    {
        for (i = 0; i < N; i++)
            res[i] = 0.0;
        rr = 1.0;
    }
    *rad = rr;
}

/* (for N = 1, the whole line fails the range test below and goes to
   _dfloat_ball_trig_slow, which gives [0 +/- 1]; the explicit check is
   for the tails of N >= 2) */
DFLOAT_FUNC_NOINLINE void DB(sin_cos)(DB(t) sn, DB(t) cs, DB(srcptr) x)
{
    if (FLINT_UNLIKELY(N > 1 && DB(_is_W)(x)))
    {
        DB(_set_zero_pm)(sn, 1.0);
        DB(_set_zero_pm)(cs, 1.0);
        return;
    }
    double xc[N], m0, rx, am;
    int idx;

    rx = x->rad;
    for (idx = 0; idx < N; idx++)
        xc[idx] = x->d[idx];
    DF(_canonicalise_inplace)(xc);
    m0 = xc[0];
    if (FLINT_UNLIKELY(!(rx <= DBL_MAX) || !(fabs(m0) <= TRIG_XMAX) || !DF(_is_canonical)(xc)))
    {
        _dfloat_ball_trig_slow(sn->d, &sn->rad, cs->d, &cs->rad, N, xc, rx, 3);
        return;
    }
    if (FLINT_UNLIKELY(DF(_trig_core)(sn->d, cs->d, xc, 3)))
    {
        _dfloat_ball_trig_slow(sn->d, &sn->rad, cs->d, &cs->rad, N, xc, rx, 3);
        return;
    }
    am = _dfloat_abs_sum(xc, N);
    DB(_trig_finish)(sn->d, &sn->rad, am, FLINT_MIN(am, 1.0), rx);
    DB(_trig_finish)(cs->d, &cs->rad, am, 1.0, rx);
}

DFLOAT_FUNC_NOINLINE void DB(sin)(DB(t) res, DB(srcptr) x)
{
    if (FLINT_UNLIKELY(N > 1 && DB(_is_W)(x)))
    {
        DB(_set_zero_pm)(res, 1.0);
        return;
    }
    double xc[N], cs[N], crad, m0, rx, am;
    int idx;

    rx = x->rad;
    for (idx = 0; idx < N; idx++)
        xc[idx] = x->d[idx];
    DF(_canonicalise_inplace)(xc);
    m0 = xc[0];
    if (FLINT_UNLIKELY(!(rx <= DBL_MAX) || !(fabs(m0) <= TRIG_XMAX) || !DF(_is_canonical)(xc)))
    {
        _dfloat_ball_trig_slow(res->d, &res->rad, cs, &crad, N, xc, rx, 1);
        return;
    }
    if (FLINT_UNLIKELY(DF(_trig_core)(res->d, cs, xc, 1)))
    {
        _dfloat_ball_trig_slow(res->d, &res->rad, cs, &crad, N, xc, rx, 1);
        return;
    }
    am = _dfloat_abs_sum(xc, N);
    DB(_trig_finish)(res->d, &res->rad, am, FLINT_MIN(am, 1.0), rx);
}

DFLOAT_FUNC_NOINLINE void DB(cos)(DB(t) res, DB(srcptr) x)
{
    if (FLINT_UNLIKELY(N > 1 && DB(_is_W)(x)))
    {
        DB(_set_zero_pm)(res, 1.0);
        return;
    }
    double xc[N], sn[N], srad, m0, rx, am;
    int idx;

    rx = x->rad;
    for (idx = 0; idx < N; idx++)
        xc[idx] = x->d[idx];
    DF(_canonicalise_inplace)(xc);
    m0 = xc[0];
    if (FLINT_UNLIKELY(!(rx <= DBL_MAX) || !(fabs(m0) <= TRIG_XMAX) || !DF(_is_canonical)(xc)))
    {
        _dfloat_ball_trig_slow(sn, &srad, res->d, &res->rad, N, xc, rx, 2);
        return;
    }
    if (FLINT_UNLIKELY(DF(_trig_core)(sn, res->d, xc, 2)))
    {
        _dfloat_ball_trig_slow(sn, &srad, res->d, &res->rad, N, xc, rx, 2);
        return;
    }
    am = _dfloat_abs_sum(xc, N);
    DB(_trig_finish)(res->d, &res->rad, am, 1.0, rx);
}

#ifdef DFLOAT_HAVE_SIMD
/* the quadrant selection on 4 lanes: from sin t, cos t and k (as a
   double), sn = +-(sin t or cos t) and cs likewise */
FLINT_FORCE_INLINE void
DF(_v_trig_quadrant)(vec4d * sn, vec4d * cs, const vec4d * st, const vec4d * ct, vec4d kd)
{
    /* q = k mod 4 in [0, 3] */
    vec4d q = vec4d_sub(kd, vec4d_mul(vec4d_set_d(4.0), vec4d_floor(vec4d_mul(kd, vec4d_set_d(0.25)))));
    vec4d odd = vec4d_bit_or(vec4d_cmp_eq(q, vec4d_set_d(1.0)), vec4d_cmp_eq(q, vec4d_set_d(3.0)));
    vec4d sneg = vec4d_bit_and(vec4d_cmp_ge(q, vec4d_set_d(2.0)), vec4d_set_d(-0.0));
    vec4d cneg = vec4d_bit_and(vec4d_bit_or(vec4d_cmp_eq(q, vec4d_set_d(1.0)), vec4d_cmp_eq(q, vec4d_set_d(2.0))), vec4d_set_d(-0.0));
    int idx;
    UNROLL
    for (idx = 0; idx < N; idx++)
    {
        sn[idx] = vec4d_bit_xor(vec4d_blendv(st[idx], ct[idx], odd), sneg);
        cs[idx] = vec4d_bit_xor(vec4d_blendv(ct[idx], st[idx], odd), cneg);
    }
}

/* the table rows sin A, cos A for A = i/64 (iv) resp. j/2^13 (jv)
   from the tables of the form [.][2][5], with the sign of the sine
   from the sign of the index */
FLINT_FORCE_INLINE void
DF(_v_trig_gather)(vec4d * sa, vec4d * ca, const double * table, vec4d iv)
{
    vec4d ai = vec4d_abs(iv);
    vec4d sign = vec4d_bit_and(iv, vec4d_set_d(-0.0));
    vec4d row = vec4d_add(ai, ai);
    int idx;
    DF(_v_gather)(sa, table, row);
    DF(_v_gather)(ca, table + 5, row);
    UNROLL
    for (idx = 0; idx < N; idx++)
        sa[idx] = vec4d_bit_xor(sa[idx], sign);
}

/* The kernel on 4 lanes: sn = sin x, cs = cos x for canonical xc with
   TRIG_XMIN <= |x0| <= TRIG_XMAX; returns the mask of lanes where the
   reduced argument failed the canonical check (to be redone by the
   scalar code). */
#if N == 1
static vec4d
DF(_v_trig_core)(vec4d * sn, vec4d * cs, const vec4d * xc)
{
    vec4d kd, h, p1, e1, t, iv, d, s, w, P, Q, sa, ca, st, ct;

    kd = vec4d_round(vec4d_mul(xc[0], vec4d_set_d(0.63661977236758138)));
    h = vec4d_sub(xc[0], vec4d_mul(kd, vec4d_set_d(_dfloat_trig_pi2_33[0])));
    p1 = vec4d_mul(kd, vec4d_set_d(_dfloat_trig_pi2_33[1]));
    e1 = vec4d_fmsub_fused(kd, vec4d_set_d(_dfloat_trig_pi2_33[1]), p1);
    t = vec4d_sub(h, p1);
    iv = vec4d_round(vec4d_mul(t, vec4d_set_d(64.0)));
    /* the clamp keeps the table index in range on the lanes that go
       to the scalar code anyway (nan, out of range) */
    iv = vec4d_max(vec4d_min(iv, vec4d_set_d(50.0)), vec4d_set_d(-50.0));
    d = vec4d_sub(t, vec4d_mul(iv, vec4d_set_d(1.0 / 64.0)));
    s = vec4d_sub(d, e1);
    w = vec4d_mul(s, s);
    P = vec4d_fmadd_fused(w, vec4d_set_d(_dfloat_trig_coef[7][0]), vec4d_set_d(_dfloat_trig_coef[5][0]));
    P = vec4d_fmadd_fused(w, P, vec4d_set_d(_dfloat_trig_coef[3][0]));
    P = vec4d_fmadd_fused(w, P, vec4d_set_d(1.0));
    Q = vec4d_fmadd_fused(w, vec4d_set_d(_dfloat_trig_coef[8][0]), vec4d_set_d(_dfloat_trig_coef[6][0]));
    Q = vec4d_fmadd_fused(w, Q, vec4d_set_d(_dfloat_trig_coef[4][0]));
    Q = vec4d_fmadd_fused(w, Q, vec4d_set_d(_dfloat_trig_coef[2][0]));
    {
        vec4d t50 = vec4d_add(iv, vec4d_set_d(50.0));
        vec4di row = vec4d_to_index(vec4d_add(t50, t50));
        sa = vec4d_gather(&_dfloat_trig_A1[0][0], row);
        ca = vec4d_gather(&_dfloat_trig_A1[0][1], row);
    }
    st = vec4d_fmadd_fused(vec4d_mul(sa, w), Q, vec4d_fmadd_fused(vec4d_mul(ca, s), P, sa));
    ct = vec4d_fmadd_fused(vec4d_mul(ca, w), Q, vec4d_fmadd_fused(vec4d_neg(vec4d_mul(sa, s)), P, ca));
    DF(_v_trig_quadrant)(sn, cs, &st, &ct, kd);
    return vec4d_zero();
}
#else
static vec4d
DF(_v_trig_core)(vec4d * sn, vec4d * cs, const vec4d * xc)
{
    vec4d Kt[N], A[N], T[N], s[N], w[N], sa[N], ca[N], sb[N], cb[N], SA[N], CA[N];
    vec4d P[N], Q[N], u[N], c[N], m0[N], m1[N], l0[N], l1[N], st[N], ct[N];
    vec4d kd, h, iv, jv, d, s0, bad;
    int kk, M, idx;

    kd = vec4d_round(vec4d_mul(xc[0], vec4d_set_d(0.63661977236758138)));
    h = vec4d_sub(xc[0], vec4d_mul(kd, vec4d_set_d(_dfloat_trig_pi2_33[0])));
    UNROLL
    for (idx = 0; idx < N; idx++)
        c[idx] = vec4d_set_d(_dfloat_trig_pi2_33[idx + 1]);
    DF(_v_mul_d_kernel)(Kt, NULL, c, kd, 0, 0);
    A[0] = h;
    UNROLL
    for (idx = 1; idx < N; idx++)
        A[idx] = xc[idx];
    DF(_v_addsub_kernel)(T, NULL, A, Kt, 1, 0, 0);
    bad = DF(_v_canonicalise)(T);
    iv = vec4d_round(vec4d_mul(T[0], vec4d_set_d(64.0)));
    iv = vec4d_max(vec4d_min(iv, vec4d_set_d(50.0)), vec4d_set_d(-50.0));
    d = vec4d_sub(T[0], vec4d_mul(iv, vec4d_set_d(1.0 / 64.0)));
    jv = vec4d_round(vec4d_mul(d, vec4d_set_d(8192.0)));
    jv = vec4d_max(vec4d_min(jv, vec4d_set_d(64.0)), vec4d_set_d(-64.0));
    s0 = vec4d_sub(d, vec4d_mul(jv, vec4d_set_d(1.0 / 8192.0)));
    s[0] = s0;
    UNROLL
    for (idx = 1; idx < N; idx++)
        s[idx] = T[idx];

    /* SA = sa cb + ca sb, CA = ca cb - sa sb */
    DF(_v_trig_gather)(sa, ca, &_dfloat_trig_A[0][0][0], iv);
    DF(_v_trig_gather)(sb, cb, &_dfloat_trig_B[0][0][0], jv);
    DF(_v_mul_kernel)(m0, NULL, sa, cb, 0, 0);
    DF(_v_mul_kernel)(m1, NULL, ca, sb, 0, 0);
    DF(_v_addsub_kernel)(SA, NULL, m0, m1, 0, 0, 0);
    DF(_v_mul_kernel)(m0, NULL, ca, cb, 0, 0);
    DF(_v_mul_kernel)(m1, NULL, sa, sb, 0, 0);
    DF(_v_addsub_kernel)(CA, NULL, m0, m1, 1, 0, 0);
    DF(_v_sqr_kernel)(w, NULL, s, 0, 0);

    /* the chains P and Q in w */
    M = DF(_trig_level)[TRIG_J2];
    for (idx = 0; idx < N; idx++)
    {
        P[idx] = vec4d_set_d((idx < M && TRIG_CS(2 * TRIG_J2 + 1) != NULL) ? TRIG_CS(2 * TRIG_J2 + 1)[idx] : 0.0);
        Q[idx] = vec4d_set_d((idx < M && TRIG_CC(2 * TRIG_J2 + 2) != NULL) ? TRIG_CC(2 * TRIG_J2 + 2)[idx] : 0.0);
    }
    for (kk = TRIG_J2 - 1; kk >= 0; kk--)
    {
        M = DF(_trig_level)[kk];
        for (idx = 0; idx < M; idx++)
            c[idx] = vec4d_set_d((TRIG_CS(2 * kk + 1) != NULL) ? TRIG_CS(2 * kk + 1)[idx] : 0.0);
        DF(_vexp_mul_level)(u, w, P, M);
        DF(_vexp_add_level)(P, c, u, M);
        for (idx = 0; idx < M; idx++)
            c[idx] = vec4d_set_d((TRIG_CC(2 * kk + 2) != NULL) ? TRIG_CC(2 * kk + 2)[idx] : 0.0);
        DF(_vexp_mul_level)(u, w, Q, M);
        DF(_vexp_add_level)(Q, c, u, M);
    }

    /* sin t = SA + (SA w) Q + (CA s) P, cos t = CA + (CA w) Q - (SA s) P */
    DF(_v_mul_kernel)(m0, NULL, SA, w, 0, 0);
    DF(_v_mul_kernel)(m1, NULL, CA, s, 0, 0);
    DF(_v_mul_kernel)(l0, NULL, m0, Q, 0, 0);
    DF(_v_mul_kernel)(l1, NULL, m1, P, 0, 0);
    DF(_v_add3_kernel)(st, NULL, SA, l0, l1, 0, 0);
    DF(_v_mul_kernel)(m0, NULL, CA, w, 0, 0);
    DF(_v_mul_kernel)(m1, NULL, SA, s, 0, 0);
    UNROLL
    for (idx = 0; idx < N; idx++)
        m1[idx] = vec4d_neg(m1[idx]);
    DF(_v_mul_kernel)(l0, NULL, m0, Q, 0, 0);
    DF(_v_mul_kernel)(l1, NULL, m1, P, 0, 0);
    DF(_v_add3_kernel)(ct, NULL, CA, l0, l1, 0, 0);

    DF(_v_trig_quadrant)(sn, cs, st, ct, kd);
    return bad;
}
#endif

/* the lanes the vector kernel does not handle: non-canonical, tiny,
   large or non-finite heads */
FLINT_FORCE_INLINE vec4d
DF(_v_trig_bad)(vec4d * xc)
{
    vec4d bad = DF(_v_canonicalise)(xc);
    vec4d ax = vec4d_abs(xc[0]);
    bad = vec4d_bit_or(bad, vec4d_cmp_nle(ax, vec4d_set_d(TRIG_XMAX)));
    bad = vec4d_bit_or(bad, vec4d_cmp_nge(ax, vec4d_set_d(TRIG_XMIN)));
    return bad;
}

/* the radii of sin and cos of [m +/- r] as in DB(_trig_finish), on
   4 lanes: EPS (|m| + a) + r with a = min(|m|, 1) for the sine, 1 for
   the cosine, the trivial enclosure [0 +/- 1] beyond 1 */
FLINT_FORCE_INLINE void
DB(_v_trig_finish)(vec4d * y, vec4d * rad, vec4d am, vec4d rx, int is_sin)
{
    vec4d a = is_sin ? vec4d_min(am, vec4d_set_d(1.0)) : vec4d_set_d(1.0);
    vec4d rr = vec4d_add(vec4d_mul(vec4d_set_d(DFLOAT_TRIG_EPS(N)), vec4d_add(am, a)), rx);
    vec4d big;
    int idx;
    rr = vec4d_mul(rr, vec4d_set_d(DFLOAT_RAD_SLACK));
    rr = vec4d_blendv(rr, vec4d_add(rr, vec4d_set_d(0x1p-1069)),
        vec4d_cmp_lt(rr, vec4d_set_d(DFLOAT_RAD_UNDERFLOW_THRESHOLD)));
    big = vec4d_cmp_nlt(rr, vec4d_set_d(1.0));
    UNROLL
    for (idx = 0; idx < N; idx++)
        y[idx] = vec4d_blendv(y[idx], vec4d_zero(), big);
    *rad = vec4d_blendv(rr, vec4d_set_d(1.0), big);
}
#endif

void DFV(vec_sin)(DF(ptr) res, DF(srcptr) x, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DF(struct) xs[4];
        vec4d xc[N], sn[N], cs[N], bad;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N);
        bad = DF(_v_trig_bad)(xc);
        bad = vec4d_bit_or(bad, DF(_v_trig_core)(sn, cs, xc));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DF(_vstore4)(res[i].d, sn, N);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DF(sin)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DF(sin)(res + i, x + i);
}

void DFV(vec_cos)(DF(ptr) res, DF(srcptr) x, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DF(struct) xs[4];
        vec4d xc[N], sn[N], cs[N], bad;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N);
        bad = DF(_v_trig_bad)(xc);
        bad = vec4d_bit_or(bad, DF(_v_trig_core)(sn, cs, xc));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DF(_vstore4)(res[i].d, cs, N);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DF(cos)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DF(cos)(res + i, x + i);
}

void DBV(vec_sin)(DB(ptr) res, DB(srcptr) x, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DB(struct) xs[4];
        vec4d xc[N], sn[N], cs[N], bad, rx, am, rad;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N + 1);
        DB(_vload_rad)(&rx, x + i);
        bad = DF(_v_trig_bad)(xc);
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(rx, vec4d_set_d(DBL_MAX)));
        bad = vec4d_bit_or(bad, DF(_v_trig_core)(sn, cs, xc));
        am = DF(_vabs_sum)(xc);
        DB(_v_trig_finish)(sn, &rad, am, rx, 1);
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DB(_vstore4)(res + i, sn, rad);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DB(sin)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DB(sin)(res + i, x + i);
}

/* both from one kernel call (sn and cs may not alias x) */
void DFV(vec_sin_cos)(DF(ptr) sres, DF(ptr) cres, DF(srcptr) x, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DF(struct) xs[4];
        vec4d xc[N], sn[N], cs[N], bad;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N);
        bad = DF(_v_trig_bad)(xc);
        bad = vec4d_bit_or(bad, DF(_v_trig_core)(sn, cs, xc));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DF(_vstore4)(sres[i].d, sn, N);
        DF(_vstore4)(cres[i].d, cs, N);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DF(sin_cos)(sres + i + l, cres + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DF(sin_cos)(sres + i, cres + i, x + i);
}

void DBV(vec_sin_cos)(DB(ptr) sres, DB(ptr) cres, DB(srcptr) x, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DB(struct) xs[4];
        vec4d xc[N], sn[N], cs[N], bad, rx, am, rad;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N + 1);
        DB(_vload_rad)(&rx, x + i);
        bad = DF(_v_trig_bad)(xc);
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(rx, vec4d_set_d(DBL_MAX)));
        bad = vec4d_bit_or(bad, DF(_v_trig_core)(sn, cs, xc));
        am = DF(_vabs_sum)(xc);
        DB(_v_trig_finish)(sn, &rad, am, rx, 1);
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DB(_vstore4)(sres + i, sn, rad);
        DB(_v_trig_finish)(cs, &rad, am, rx, 0);
        DB(_vstore4)(cres + i, cs, rad);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DB(sin_cos)(sres + i + l, cres + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DB(sin_cos)(sres + i, cres + i, x + i);
}

void DBV(vec_cos)(DB(ptr) res, DB(srcptr) x, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DB(struct) xs[4];
        vec4d xc[N], sn[N], cs[N], bad, rx, am, rad;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N + 1);
        DB(_vload_rad)(&rx, x + i);
        bad = DF(_v_trig_bad)(xc);
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(rx, vec4d_set_d(DBL_MAX)));
        bad = vec4d_bit_or(bad, DF(_v_trig_core)(sn, cs, xc));
        am = DF(_vabs_sum)(xc);
        DB(_v_trig_finish)(cs, &rad, am, rx, 0);
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DB(_vstore4)(res + i, cs, rad);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DB(cos)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DB(cos)(res + i, x + i);
}

/* ---- logarithm ----

   With x = 2^e m, m in [M0, 2 M0), M0 = 0.70703125, and two reciprocal
   tables: r = 1/m_mid rounded for the 128 intervals of the bit
   pattern of m (index ix = the 7 bits below the exponent of bits(m)
   - bits(M0); width 2^-8 in m below 1, 2^-7 above; r = 1 on the two
   intervals adjacent to 1), so that z = m r - 1 has |z| <= 2^-7;
   r'_j = 1/(1 + j/2^13) rounded (j = round(2^13 z) from a double
   approximation of z), and z2 = m (r r') - 1 with |z2| <= 2^-14 (1 +
   2^-7), computed in one go from the exact double-double product
   r r' (no division, and both tables compensate their rounded
   entries exactly: log x = e log 2 + L + L'_j + log1p(z2) with L =
   -log r).  The reduction is exact where r = r' = 1, which is what
   preserves the relative accuracy near x = 1; its absolute error is
   otherwise about u^N |z2| against a result of at least 2^-8 (r != 1)
   or about |z| >= 2^-14 (r' != 1).  With w = z2^2, log1p(z2) = z2 +
   w E(w) + z2 w O(w),
   E(w) = -1/2 - w/4 - ..., O(w) = 1/3 + w/5 + ... (degree LOG_DEG
   in z2), evaluated as four chains in v = w^2 (scalar) or two in w
   (vector), and the result is the five-term sum of D = e log 2 + L +
   L'_j + z2 (assembled while the chains run) and the four products.
   The guarantee is |res - log x| <= EPS |log x| (dev/dfloat_exp_bound.py). */
#if N == 1
#define LOG_DEG 8
#elif N == 2
#define LOG_DEG 8
#define LOG_J2 3
#define LOG_J4 1
static const int DF(_log_level)[LOG_J2 + 1] = {2, 2, 1, 1};
static const int DF(_log_level4)[LOG_J4 + 1] = {2, 1};
#elif N == 3
#define LOG_DEG 12
#define LOG_J2 5
#define LOG_J4 2
static const int DF(_log_level)[LOG_J2 + 1] = {3, 3, 2, 2, 1, 1};
static const int DF(_log_level4)[LOG_J4 + 1] = {3, 2, 1};
#else
#define LOG_DEG 15
#define LOG_J2 7
#define LOG_J4 3
static const int DF(_log_level)[LOG_J2 + 1] = {4, 4, 3, 3, 2, 2, 1, 1};
static const int DF(_log_level4)[LOG_J4 + 1] = {4, 3, 2, 1};
#endif
/* the coefficient of z^k (a 5-term row), or NULL beyond the degree */
#define LOG_C(k) (((k) <= LOG_DEG) ? _dfloat_log_coef[k] : NULL)
/* below this |z2|, log1p(z2) = z2 (a relative error of at most
   |z2|/2 <= 2^-251; above it the powers of z2 in the chains stay
   normal, as the bound derivation assumes) */
#define LOG_ZMIN 0x1p-250

/* the exponent e, the scaling 2^-e (as a double) and the level-1
   index of a double x0 in [2^-1022, 2^1022): x0 = 2^e m with m in
   [M0, 2 M0), from the bits of x0 less those of M0 (the exponent
   field of the difference is e, the 7 bits below it the index) */
#define LOG_M0_BITS UINT64_C(0x3fe6a00000000000)
#define LOG_EXPONENT(x0, e, scale, ix) \
    do { \
        double_uint64_u _u; \
        uint64_t _t; \
        _u.f = (x0); \
        _t = _u.i - LOG_M0_BITS; \
        (e) = (int) ((int64_t) _t >> 52); \
        (ix) = (int) ((_t >> 45) & 127); \
        _u.i = ((uint64_t) (1023 - (e))) << 52; \
        (scale) = _u.f; \
    } while (0)

#if N == 1
/* the double kernel: one table level, z = fma(m, r, -1) rounded once
   (a relative error; exact where r = 1), the polynomial of degree 8
   by fma, and the sum (e log 2 + L) + e L1 + z + w P in that order */
FLINT_FORCE_INLINE int
DF(_log_kernel)(double * res, double m, double t, int e, int ix)
{
    double r, z, w, P, d;

    r = _dfloat_log_R1[ix];
    /* z = (m + t) r - 1 with t the low part of a double-double
       argument (log1p: 0 for log) */
    z = fma(m, r, -1.0) + t * r;
    w = z * z;
    /* log1p(z) = z + w (c2 + c3 z + w (c4 + c5 z + w (c6 + c7 z + w c8))) */
    P = fma(w, _dfloat_log_coef[8][0], fma(z, _dfloat_log_coef[7][0], _dfloat_log_coef[6][0]));
    P = fma(w, P, fma(z, _dfloat_log_coef[5][0], _dfloat_log_coef[4][0]));
    P = fma(w, P, fma(z, _dfloat_log_coef[3][0], _dfloat_log_coef[2][0]));
    d = (double) e * _dfloat_exp_ln2_42[0] + _dfloat_log_L1[ix][0];
    d = fma((double) e, _dfloat_exp_ln2_42[1], d) + z;
    res[0] = fma(w, P, d);
    return 0;
}

FLINT_FORCE_INLINE int
DF(_log_core)(double * res, const double * xc)
{
    double scale;
    int e, ix;
    LOG_EXPONENT(xc[0], e, scale, ix);
    return DF(_log_kernel)(res, xc[0] * scale, 0.0, e, ix);
}
#else
/* res = log(2^e m) for m = m0 + tail in [M0, 2 M0), with tail (N
   terms, of size u m0 or below) and mfull = the first N terms of
   (m0, tail) the two forms of the argument used by the reduction;
   e and the table index i from LOG_EXPONENT of m0.  For log the tail
   is (m1, ..., m_{N-1}, 0) and mfull = m; for log1p the argument
   1 + x is carried as m0 = fl(1 + x0) and the tail (1 + x0 - m0, x1,
   ..., x_{N-1}), so that it is exact (the dropped x_{N-1} Rl is of
   size u^(N+1)).  Returns nonzero if the reduced argument failed the
   canonical check (the caller then uses arb). */
static int
DF(_log_kernel)(double * res, double m0, const double * tail, const double * mfull, int e, int i)
{
    double r, ed, za;
    double P[N], z[N], z2[N], w[N], c[N], eL[N], D[N], l0[N], l1[N], l2[N], l3[N];
    int j, idx;

    r = _dfloat_log_R1[i];
    /* the level-2 index from the double approximation of z = m r - 1
       (any j with |z - j/2^13| <= 2^-14 + 2u serves: L'_j compensates
       the rounded r'_j exactly) */
    za = fma(m0, r, -1.0);
    j = (int) rint(za * 8192.0);
    j = FLINT_MAX(FLINT_MIN(j, 64), -64);
    /* z2 = m (r r') - 1 with r r' = Rh + Rl exactly: z2 = (h0 - 1) +
       l0 + (m1 + ...) Rh + m Rl with m0 Rh = h0 + l0 (TwoProd) and
       h0 - 1 exact (h0 in [1 - 2^-13, 1 + 2^-13]): the head of the
       product is cancelled exactly and the other terms are of size
       u, computed with N components (absolute errors u^N u), so that
       the N components of z2 carry the accuracy u^N |z2| (a product
       rounded to N components would only give u^N); both products
       are exact where r = r' = 1 */
    {
        double Rh, Rl, h0, l0;
        DFLOAT_TWO_PROD(Rh, Rl, r, _dfloat_log_R2[j + 64]);
        DFLOAT_TWO_PROD(h0, l0, m0, Rh);
        c[0] = h0 - 1.0;
        D[0] = l0;
        for (idx = 1; idx < N; idx++)
            c[idx] = D[idx] = 0.0;
        DF(_mul_d_kernel)(P, NULL, tail, Rh, 0);
        DF(_mul_d_kernel)(z, NULL, mfull, Rl, 0);
        {
            const double * xs[4] = {c, D, P, z};
            if (FLINT_UNLIKELY(DF(_s_addn_kernel)(z2, NULL, xs, 4, 0, 0)))
                DF(_shift_zero_head)(z2);
        }
    }
    DF(_canonicalise_inplace)(z2);
    if (FLINT_UNLIKELY(!DF(_is_canonical)(z2)))
        return 1;

    /* D = e L0 + e (L1 + ... + LN) + L_i + L'_j + z2 (e L0 is exact:
       L0 has 42 bits and |e| < 2^11; the tail L1 + ... is a
       canonical expansion by itself, unlike (L0, L1, ...)) */
    ed = (double) e;
    c[0] = ed * _dfloat_exp_ln2_42[0];
    for (idx = 1; idx < N; idx++)
        c[idx] = 0.0;
    DF(_mul_d_kernel)(eL, NULL, _dfloat_exp_ln2_42 + 1, ed, 0);
    {
        const double * xs[5] = {c, eL, _dfloat_log_L1[i], _dfloat_log_L2[j + 64], z2};
        if (FLINT_UNLIKELY(DF(_s_addn_kernel)(D, NULL, xs, 5, 0, 0)))
            DF(_shift_zero_head)(D);
    }
    if (FLINT_UNLIKELY(fabs(z2[0]) < LOG_ZMIN))
    {
        for (idx = 0; idx < N; idx++)
            res[idx] = D[idx];
        return 0;
    }

    DF(_sqr_kernel)(w, NULL, z2, 0);
#ifdef DFLOAT_HAVE_SIMD
    {
        /* the chains E0, E1, O0, O1 in v = w^2, the multipliers
           (w, v, z2 w, z2 v), and the five-term sum */
        vec4d vv[N], vq[N], va[N], vu[N];
        const double * ctab[4 * (LOG_J4 + 1)];
        double v[N], zw[N], zv[N];
        double t4[4];
        int kk;

        DF(_sqr_kernel)(v, NULL, w, 0);
        for (kk = 0; kk <= LOG_J4; kk++)
        {
            /* E0: c2, c6, ...; E1: c4, c8, ...; O0: c3, c7, ...; O1: c5, c9, ... */
            ctab[0 * (LOG_J4 + 1) + kk] = LOG_C(4 * kk + 2);
            ctab[1 * (LOG_J4 + 1) + kk] = LOG_C(4 * kk + 4);
            ctab[2 * (LOG_J4 + 1) + kk] = LOG_C(4 * kk + 3);
            ctab[3 * (LOG_J4 + 1) + kk] = LOG_C(4 * kk + 5);
        }
        for (idx = 0; idx < N; idx++)
            vv[idx] = vec4d_set_d(v[idx]);
        DF(_horner4)(vq, vv, ctab, LOG_J4, DF(_log_level4));
        DF(_mul_kernel)(zw, NULL, z2, w, 0);
        DF(_mul_kernel)(zv, NULL, z2, v, 0);
        for (idx = 0; idx < N; idx++)
            va[idx] = vec4d_set_d4(w[idx], v[idx], zw[idx], zv[idx]);
        DF(_v_mul_kernel)(vu, NULL, va, vq, 0, 0);
        for (idx = 0; idx < N; idx++)
        {
            vec4d_store_unaligned(t4, vu[idx]);
            l0[idx] = t4[0]; l1[idx] = t4[1]; l2[idx] = t4[2]; l3[idx] = t4[3];
        }
        {
            const double * xs[5] = {D, l0, l1, l2, l3};
            if (FLINT_UNLIKELY(DF(_s_addn_kernel)(res, NULL, xs, 5, 0, 0)))
                DF(_shift_zero_head)(res);
        }
    }
#else
    {
        /* the chains E and O in w, sequentially */
        double E[N], O[N], u[N], zw[N];
        int kk, M;
        M = DF(_log_level)[LOG_J2];
        for (idx = 0; idx < N; idx++)
        {
            E[idx] = (idx < M && LOG_C(2 * LOG_J2 + 2) != NULL) ? LOG_C(2 * LOG_J2 + 2)[idx] : 0.0;
            O[idx] = (idx < M && LOG_C(2 * LOG_J2 + 3) != NULL) ? LOG_C(2 * LOG_J2 + 3)[idx] : 0.0;
        }
        for (kk = LOG_J2 - 1; kk >= 0; kk--)
        {
            M = DF(_log_level)[kk];
            for (idx = 0; idx < N; idx++)
                c[idx] = (idx < M && LOG_C(2 * kk + 2) != NULL) ? LOG_C(2 * kk + 2)[idx] : 0.0;
            DF(_level_mul)(u, w, E, M);
            DF(_level_addsub)(E, c, u, 0, M);
            for (idx = 0; idx < N; idx++)
                c[idx] = (idx < M && LOG_C(2 * kk + 3) != NULL) ? LOG_C(2 * kk + 3)[idx] : 0.0;
            DF(_level_mul)(u, w, O, M);
            DF(_level_addsub)(O, c, u, 0, M);
        }
        DF(_mul_kernel)(zw, NULL, z2, w, 0);
        DF(_mul_kernel)(l0, NULL, w, E, 0);
        DF(_mul_kernel)(l1, NULL, zw, O, 0);
        if (FLINT_UNLIKELY(DF(_s_add3_kernel)(res, NULL, D, l0, l1, 0, 0)))
            DF(_shift_zero_head)(res);
        (void) l2; (void) l3;
    }
#endif
    return 0;
}

/* the argument (2^e m with m in [M0, 2 M0)) of a canonical x with
   2^-1022 <= x0 < 2^1022 */
static int
DF(_log_core)(double * res, const double * xc)
{
    double scale, m[N], tail[N];
    int e, i, idx;
    LOG_EXPONENT(xc[0], e, scale, i);
    for (idx = 0; idx < N; idx++)
        m[idx] = xc[idx] * scale;
    for (idx = 0; idx + 1 < N; idx++)
        tail[idx] = m[idx + 1];
    tail[N - 1] = 0.0;
    return DF(_log_kernel)(res, m[0], tail, m, e, i);
}
#endif

DFLOAT_FUNC_NOINLINE void DF(log)(DF(t) res, DF(srcptr) x)
{
    double xc[N];
    int i;

    for (i = 0; i < N; i++)
        xc[i] = x->d[i];
    DF(_canonicalise_inplace)(xc);
    if (FLINT_UNLIKELY(!(xc[0] >= 0x1p-1022) || !(xc[0] < 0x1p1022) || !DF(_is_canonical)(xc)))
    {
        _dfloat_log_slow(res->d, N, xc);
        return;
    }
    if (FLINT_UNLIKELY(DF(_log_core)(res->d, xc)))
        _dfloat_log_slow(res->d, N, xc);
}

/* log1p: the argument 1 + x carried exactly as (s, t, x1, ...) with
   s + t = 1 + x0 (TwoSum); returns nonzero if 1 + x is out of the
   kernel's range or the kernel failed */
/* shared by the plain and ball functions; inlined for N <= 2, where the
   call would cost 5-15% */
#if N >= 3
static DFLOAT_NOINLINE int
#else
FLINT_FORCE_INLINE int
#endif
DF(_log1p_core)(double * res, const double * xc)
{
    double s, t, scale;
    int e, i;
#if N > 1
    if (fabs(xc[0]) > 0.5)
    {
        /* |log(1 + x)| >= 0.4: the argument 1 + x rounded to N
           components (a relative error u^N, an absolute one in the
           result) and the plain kernel */
        double y[N], one[N];
        int idx;
        one[0] = 1.0;
        for (idx = 1; idx < N; idx++)
            one[idx] = 0.0;
        DF(_add_kernel)(y, NULL, one, xc, 0);
        DF(_canonicalise_inplace)(y);
        if (FLINT_UNLIKELY(!(y[0] >= 0x1p-1022) || !(y[0] < 0x1p1022) || !DF(_is_canonical)(y)))
            return 1;
        return DF(_log_core)(res, y);
    }
#endif
    /* |x0| <= 1/2: 1 + x = s + t + x1 + ... exactly, with s in [1/2,
       3/2] and the tail of the size u s or below */
    DFLOAT_TWO_SUM(s, t, 1.0, xc[0]);
    if (FLINT_UNLIKELY(!(s >= 0x1p-1022) || !(s < 0x1p1022)))
        return 1;
    LOG_EXPONENT(s, e, scale, i);
#if N == 1
    return DF(_log_kernel)(res, s * scale, t * scale, e, i);
#else
    {
        double tail[N], mfull[N];
        int idx;
        tail[0] = t * scale;
        for (idx = 1; idx < N; idx++)
            tail[idx] = xc[idx] * scale;
        mfull[0] = s * scale;
        for (idx = 1; idx < N; idx++)
            mfull[idx] = tail[idx - 1];
        return DF(_log_kernel)(res, s * scale, tail, mfull, e, i);
    }
#endif
}

DFLOAT_FUNC_NOINLINE void DF(log1p)(DF(t) res, DF(srcptr) x)
{
    double xc[N];
    int i;

    for (i = 0; i < N; i++)
        xc[i] = x->d[i];
    DF(_canonicalise_inplace)(xc);
    if (FLINT_UNLIKELY(!(fabs(xc[0]) <= DBL_MAX) || !DF(_is_canonical)(xc) || DF(_log1p_core)(res->d, xc)))
        _dfloat_log1p_slow(res->d, N, xc);
}

/* the logarithm of [m +/- r]: the kernel on m and the radius
   EPS |y| + r / (m - r) (the derivative 1/t is largest at the left
   end); GR_UNABLE unless m - r > 0 */
DFLOAT_FUNC_NOINLINE int DB(log)(DB(t) res, DB(srcptr) x)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x)))
        return GR_UNABLE;
    double xc[N], m0, rx, mlow, rr;
    int idx;

    rx = x->rad;
    for (idx = 0; idx < N; idx++)
        xc[idx] = x->d[idx];
    DF(_canonicalise_inplace)(xc);
    m0 = xc[0];
    /* a lower bound on |x| over the ball: the head less 2u of it
       (the tail) less the radius, rounded down */
    mlow = _dfloat_lower((m0 - 0x1p-52 * m0) - rx);
    if (FLINT_UNLIKELY(!(rx <= DBL_MAX) || !(mlow > 0.0) || !(m0 >= 0x1p-1022) || !(m0 < 0x1p1022) || !DF(_is_canonical)(xc)))
    {
        if (DB_ALL_LE(x, 0.0))
            return GR_DOMAIN;
        return _dfloat_ball_log_slow(res->d, &res->rad, N, xc, rx);
    }
    if (FLINT_UNLIKELY(DF(_log_core)(res->d, xc)))
        return _dfloat_ball_log_slow(res->d, &res->rad, N, xc, rx);
    rr = DFLOAT_LOG_EPS(N) * _dfloat_abs_sum(res->d, N) + rx / mlow;
    rr = _dfloat_rad_finish(rr);
    if (FLINT_UNLIKELY(rr < DFLOAT_RAD_UNDERFLOW_THRESHOLD))
    {
        /* log 1 = 0 exactly; otherwise a tiny result may have lost
           accuracy to underflow */
        if (m0 == 1.0 && rx == 0.0 && DF(_tail_abs_sum)(xc) == 0.0)
        {
            DB(_zero_inline)(res);
            rr = 0.0;
        }
        else
            rr += 0x1p-1069;
    }
    res->rad = rr;
    DB(_fix_nonfinite)(res->d, &res->rad);
    return GR_SUCCESS;
}

/* log1p of [m +/- r]: the radius EPS |y| + r / (1 + m - r) */
DFLOAT_FUNC_NOINLINE int DB(log1p)(DB(t) res, DB(srcptr) x)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x)))
        return GR_UNABLE;
    double xc[N], m0, rx, mlow, rr;
    int idx;

    rx = x->rad;
    for (idx = 0; idx < N; idx++)
        xc[idx] = x->d[idx];
    DF(_canonicalise_inplace)(xc);
    m0 = xc[0];
    /* a lower bound on 1 + x over the ball */
    mlow = _dfloat_lower(((1.0 + m0) - 0x1p-52 * fabs(m0)) - rx);
    if (FLINT_UNLIKELY(!(rx <= DBL_MAX) || !(mlow > 0.0) || !(fabs(m0) <= DBL_MAX) || !DF(_is_canonical)(xc)))
    {
        if (DB_ALL_LE(x, -1.0))
            return GR_DOMAIN;
        return _dfloat_ball_log1p_slow(res->d, &res->rad, N, xc, rx);
    }
    if (m0 == 0.0 && _dfloat_abs_sum(xc, N) == 0.0)
    {
        /* |log1p(x')| <= r / (1 - r); tiny r needs the bump, since
           r (1 + r) rounds to r */
        DB(_zero_inline)(res);
        rr = (rx == 0.0) ? 0.0 : _dfloat_rad_finish(rx / mlow);
        if (FLINT_UNLIKELY(rr != 0.0 && rr < DFLOAT_RAD_UNDERFLOW_THRESHOLD))
            rr += 0x1p-1069;
        res->rad = rr;
        return GR_SUCCESS;
    }
    if (FLINT_UNLIKELY(DF(_log1p_core)(res->d, xc)))
        return _dfloat_ball_log1p_slow(res->d, &res->rad, N, xc, rx);
    rr = DFLOAT_LOG1P_EPS(N) * _dfloat_abs_sum(res->d, N) + rx / mlow;
    rr = _dfloat_rad_finish(rr);
    if (FLINT_UNLIKELY(rr < DFLOAT_RAD_UNDERFLOW_THRESHOLD))
        rr += 0x1p-1069;
    res->rad = rr;
    DB(_fix_nonfinite)(res->d, &res->rad);
    return GR_SUCCESS;
}

/* ---- arctangent ----

   For x >= 0 (odd symmetry): with c = i/64, i = round(64 x) for x <= 1
   and i = round(64/x) for x > 1, atan x = atan c + atan z with z = (x
   - c)/(1 + x c) (x <= 1), or atan x = (pi/2 - atan c) - atan z with
   z = (1 - c x)/(x + c) (x > 1), |z| <= 2^-7.  One division; the
   numerator x - c is exact and the head cancellation of 1 - c x
   exact, so that the absolute error of z is about u^N |z| for x <=
   1 (relative, since the result is at least atan c - |z| >= 2^-7
   when i > 0, and z = x exactly when i = 0, where no division is
   done) and about u^N for x > 1 (against a result of at least
   pi/4).  With w = z^2, atan z = z + z w E(w), E(w) = -1/3 + w/5 -
   ..., evaluated as four chains in u = w^4 (scalar; the multipliers
   z w, z w^2, z w^3, z w^4 are one vector product) or two in w
   (vector).  A second reduction level would shorten the chains by
   two tapered steps at the cost of a second division, which is not
   worth it (see the documentation).  The guarantee is |res - atan x|
   <= EPS |atan x|; for |x| < 2^-200 the result is x and for |x| >=
   2^900 it is pi/2 - 1/x0. */
#if N == 1
#define ATAN_K 2   /* coefficients of E: k = 0..2 (z^3, z^5, z^7) */
#elif N == 2
#define ATAN_K 6
#define ATAN_J2 3
#define ATAN_J8 1
static const int DF(_atan_level)[ATAN_J2 + 1] = {2, 2, 1, 1};
static const int DF(_atan_level8)[ATAN_J8 + 1] = {2, 1};
#elif N == 3
#define ATAN_K 10
#define ATAN_J2 5
#define ATAN_J8 2
static const int DF(_atan_level)[ATAN_J2 + 1] = {3, 3, 2, 2, 1, 1};
static const int DF(_atan_level8)[ATAN_J8 + 1] = {3, 2, 1};
#else
#define ATAN_K 13
#define ATAN_J2 6
#define ATAN_J8 3
static const int DF(_atan_level)[ATAN_J2 + 1] = {4, 4, 3, 3, 2, 2, 1};
static const int DF(_atan_level8)[ATAN_J8 + 1] = {4, 3, 2, 1};
#endif
/* the coefficient of z^(2k+3) in atan (the coefficient k + 1 of the
   table), or NULL beyond the degree */
#define ATAN_C(k) (((k) <= ATAN_K) ? _dfloat_atan_coef[(k) + 1] : NULL)
#define ATAN_XMIN 0x1p-200
#define ATAN_XMAX 0x1p900

#if N == 1
/* the double kernel: the reduction with a double division, the
   polynomial by fma */
FLINT_FORCE_INLINE int
DF(_atan_core)(double * res, const double * xc)
{
    double x, c, z, w, E, C;
    int i, big;

    x = fabs(xc[0]);
    big = (x > 1.0);
    i = (int) rint(big ? 64.0 / x : 64.0 * x);
    c = i * (1.0 / 64.0);
    if (big)
    {
        z = fma(-c, x, 1.0) / (x + c);
        C = _dfloat_pi_2[0][0] - _dfloat_atan_C[i][0];
        z = -z;
    }
    else
    {
        z = (i == 0) ? x : (x - c) / fma(x, c, 1.0);
        C = _dfloat_atan_C[i][0];
    }
    w = z * z;
    E = fma(w, _dfloat_atan_coef[3][0], _dfloat_atan_coef[2][0]);
    E = fma(w, E, _dfloat_atan_coef[1][0]);
    res[0] = (C + fma(z * w, E, z)) * ((xc[0] < 0.0) ? -1.0 : 1.0);
    return 0;
}
#else
/* res = atan x for a canonical x with 2^-200 <= |x0| < 2^900; returns
   nonzero if the reduced argument failed the canonical check */
static int
DF(_atan_core)(double * res, const double * xc)
{
    double c, sg;
    double x[N], n[N], d[N], z[N], w[N], t[N], C[N], D[N], l0[N], l1[N], l2[N], l3[N];
    int i, idx, big;

    sg = (xc[0] < 0.0) ? -1.0 : 1.0;
    for (idx = 0; idx < N; idx++)
        x[idx] = xc[idx] * sg;
    big = (x[0] > 1.0);
    i = (int) rint(big ? 64.0 / x[0] : 64.0 * x[0]);
    i = FLINT_MAX(FLINT_MIN(i, 64), 0);
    c = i * (1.0 / 64.0);
    if (big)
    {
        /* z = (1 - c x)/(x + c), atan x = (pi/2 - atan c) + atan(-z) */
        DF(_mul_d_kernel)(t, NULL, x, c, 0);
        n[0] = 1.0 - t[0];
        for (idx = 1; idx < N; idx++)
            n[idx] = -t[idx];
        /* the head has cancelled: renormalise (for N > 2 the
           canonicalisation's TwoSum chain does it) */
        if (N == 2)
            DFLOAT_TWO_SUM(n[0], n[1], n[0], n[1]);
        DF(_canonicalise_inplace)(n);
        t[0] = c;
        for (idx = 1; idx < N; idx++)
            t[idx] = 0.0;
        DF(_add_kernel)(d, NULL, x, t, 0);
        DF(_div_kernel)(z, n, d);
        for (idx = 0; idx < N; idx++)
            z[idx] = -z[idx];
        DF(_sub_kernel)(C, NULL, _dfloat_pi_2[0], _dfloat_atan_C[i], 0);
    }
    else if (i == 0)
    {
        for (idx = 0; idx < N; idx++)
        {
            z[idx] = x[idx];
            C[idx] = 0.0;
        }
    }
    else
    {
        /* z = (x - c)/(1 + x c): x0 - c exact (Sterbenz) */
        n[0] = x[0] - c;
        for (idx = 1; idx < N; idx++)
            n[idx] = x[idx];
        if (N == 2)
            DFLOAT_TWO_SUM(n[0], n[1], n[0], n[1]);
        DF(_canonicalise_inplace)(n);
        DF(_mul_d_kernel)(t, NULL, x, c, 0);
        C[0] = 1.0;
        for (idx = 1; idx < N; idx++)
            C[idx] = 0.0;
        DF(_add_kernel)(d, NULL, C, t, 0);
        DF(_div_kernel)(z, n, d);
        for (idx = 0; idx < N; idx++)
            C[idx] = _dfloat_atan_C[i][idx];
    }
    DF(_canonicalise_inplace)(z);
    if (FLINT_UNLIKELY(!DF(_is_canonical)(z)))
        return 1;
    /* D = C + z */
    DF(_add_kernel)(D, NULL, C, z, 0);
    DF(_sqr_kernel)(w, NULL, z, 0);

#ifdef DFLOAT_HAVE_SIMD
    {
        /* the chains E0..E3 in u = w^4, the multipliers p (1, w, w^2,
           w^3) with p = z w, and the five-term sum */
        vec4d vv[N], vq[N], va[N], vb[N], vu[N];
        const double * ctab[4 * (ATAN_J8 + 1)];
        double w2[N], w3[N], u[N], p[N];
        double t4[4];
        int kk;

        DF(_sqr_kernel)(w2, NULL, w, 0);
        DF(_sqr_kernel)(u, NULL, w2, 0);
        for (kk = 0; kk <= ATAN_J8; kk++)
        {
            ctab[0 * (ATAN_J8 + 1) + kk] = ATAN_C(4 * kk);
            ctab[1 * (ATAN_J8 + 1) + kk] = ATAN_C(4 * kk + 1);
            ctab[2 * (ATAN_J8 + 1) + kk] = ATAN_C(4 * kk + 2);
            ctab[3 * (ATAN_J8 + 1) + kk] = ATAN_C(4 * kk + 3);
        }
        for (idx = 0; idx < N; idx++)
            vv[idx] = vec4d_set_d(u[idx]);
        DF(_horner4)(vq, vv, ctab, ATAN_J8, DF(_atan_level8));
        DF(_mul_kernel)(p, NULL, z, w, 0);
        DF(_mul_kernel)(w3, NULL, w, w2, 0);
        for (idx = 0; idx < N; idx++)
        {
            va[idx] = vec4d_set_d(p[idx]);
            vb[idx] = vec4d_set_d4((idx == 0) ? 1.0 : 0.0, w[idx], w2[idx], w3[idx]);
        }
        DF(_v_mul_kernel)(vu, NULL, va, vb, 0, 0);
        DF(_v_mul_kernel)(va, NULL, vu, vq, 0, 0);
        for (idx = 0; idx < N; idx++)
        {
            vec4d_store_unaligned(t4, va[idx]);
            l0[idx] = t4[0]; l1[idx] = t4[1]; l2[idx] = t4[2]; l3[idx] = t4[3];
        }
        {
            const double * xs[5] = {D, l0, l1, l2, l3};
            if (FLINT_UNLIKELY(DF(_s_addn_kernel)(res, NULL, xs, 5, 0, 0)))
                DF(_shift_zero_head)(res);
        }
    }
#else
    {
        /* the chain E in w, sequentially */
        double E[N], u[N], cc[N], p[N];
        int kk, M;
        M = DF(_atan_level)[ATAN_J2];
        for (idx = 0; idx < N; idx++)
        {
            E[idx] = (idx < M && ATAN_C(2 * ATAN_J2) != NULL) ? ATAN_C(2 * ATAN_J2)[idx] : 0.0;
            l1[idx] = (idx < M && ATAN_C(2 * ATAN_J2 + 1) != NULL) ? ATAN_C(2 * ATAN_J2 + 1)[idx] : 0.0;
        }
        /* two chains E0(v), E1(v) in v = w^2 (E = E0 + w E1) */
        DF(_sqr_kernel)(u, NULL, w, 0);
        for (kk = ATAN_J2 - 1; kk >= 0; kk--)
        {
            M = DF(_atan_level)[kk];
            for (idx = 0; idx < N; idx++)
                cc[idx] = (idx < M && ATAN_C(2 * kk) != NULL) ? ATAN_C(2 * kk)[idx] : 0.0;
            DF(_level_mul)(p, u, E, M);
            DF(_level_addsub)(E, cc, p, 0, M);
            for (idx = 0; idx < N; idx++)
                cc[idx] = (idx < M && ATAN_C(2 * kk + 1) != NULL) ? ATAN_C(2 * kk + 1)[idx] : 0.0;
            DF(_level_mul)(p, u, l1, M);
            DF(_level_addsub)(l1, cc, p, 0, M);
        }
        DF(_mul_kernel)(p, NULL, z, w, 0);
        DF(_mul_kernel)(l0, NULL, p, E, 0);
        DF(_mul_kernel)(p, NULL, p, w, 0);
        DF(_mul_kernel)(l2, NULL, p, l1, 0);
        if (FLINT_UNLIKELY(DF(_s_add3_kernel)(res, NULL, D, l0, l2, 0, 0)))
            DF(_shift_zero_head)(res);
        (void) l3;
    }
#endif
    for (idx = 0; idx < N; idx++)
        res[idx] *= sg;
    return 0;
}
#endif

/* the tiny and huge cases, and the dispatch to the kernel; xc is
   canonical with a finite head; returns nonzero if the kernel failed */
FLINT_FORCE_INLINE int
DF(_atan_eval)(double * res, const double * xc)
{
    double ax = fabs(xc[0]);
    int idx;
    if (FLINT_UNLIKELY(ax < ATAN_XMIN))
    {
        for (idx = 0; idx < N; idx++)
            res[idx] = xc[idx];
        return 0;
    }
    if (FLINT_UNLIKELY(!(ax < ATAN_XMAX)))
    {
        /* pi/2 - 1/x0 (the tail of x is below 2^847 and 1/x0 below
           2^-900, so the two roundings and the neglected terms are far
           below u^N pi/2) */
        double t[N];
        t[0] = -1.0 / xc[0];
        for (idx = 1; idx < N; idx++)
            t[idx] = 0.0;
        DF(_add_kernel)(res, NULL, _dfloat_pi_2[0], t, 0);
        if (xc[0] < 0.0)
            for (idx = 0; idx < N; idx++)
                res[idx] = -res[idx];
        return 0;
    }
    return DF(_atan_core)(res, xc);
}

DFLOAT_FUNC_NOINLINE void DF(atan)(DF(t) res, DF(srcptr) x)
{
    double xc[N];
    int i;

    for (i = 0; i < N; i++)
        xc[i] = x->d[i];
    DF(_canonicalise_inplace)(xc);
    if (FLINT_UNLIKELY(!(fabs(xc[0]) <= DBL_MAX) || !DF(_is_canonical)(xc)))
    {
        _dfloat_atan_slow(res->d, N, xc);
        return;
    }
    if (FLINT_UNLIKELY(DF(_atan_eval)(res->d, xc)))
        _dfloat_atan_slow(res->d, N, xc);
}

/* the arctangent of [m +/- r]: the kernel on m and the radius
   EPS |y| + r (the derivative is at most 1) */
DFLOAT_FUNC_NOINLINE void DB(atan)(DB(t) res, DB(srcptr) x)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x)))
    {
        DB(_set_zero_pm)(res, 0x1.921fb54442d19p+0);   /* pi/2 rounded up */
        return;
    }
    double xc[N], rx, rr;
    int idx;

    rx = x->rad;
    for (idx = 0; idx < N; idx++)
        xc[idx] = x->d[idx];
    DF(_canonicalise_inplace)(xc);
    if (FLINT_UNLIKELY(!(rx <= DBL_MAX) || !(fabs(xc[0]) <= DBL_MAX) || !DF(_is_canonical)(xc)))
    {
        _dfloat_ball_atan_slow(res->d, &res->rad, N, xc, rx);
        return;
    }
    if (FLINT_UNLIKELY(DF(_atan_eval)(res->d, xc)))
    {
        _dfloat_ball_atan_slow(res->d, &res->rad, N, xc, rx);
        return;
    }
    rr = DFLOAT_ATAN_EPS(N) * _dfloat_abs_sum(res->d, N) + rx;
    rr = _dfloat_rad_finish(rr);
    if (FLINT_UNLIKELY(rr < DFLOAT_RAD_UNDERFLOW_THRESHOLD))
        rr += 0x1p-1069;
    res->rad = rr;
}

#ifdef DFLOAT_HAVE_SIMD
/* an integer-valued vec4d from 64-bit integer lanes in [-2^51, 2^51)
   (the bits of 1.5 2^52 + k represent that number exactly) */
FLINT_FORCE_INLINE vec4d
DF(_v_epi64_to_pd)(vec4d k)
{
    return vec4d_sub(vec4d_bits_add(k, vec4d_bits_set_u64(UINT64_C(0x4338000000000000))), vec4d_set_d(0x1.8p52));
}

/* the logarithm kernel on 4 lanes (N >= 2 only: the N = 1 vector
   function uses the scalar code, which is as fast); returns the
   mask of lanes to redo in scalar (a failed canonical check, a tiny
   z2) */
#if N >= 2
static vec4d
DF(_v_log_core)(vec4d * res, const vec4d * xc)
{
    vec4d m[N], P[N], z[N], z2[N], w[N], c[N], D[N], eL[N], q[N], r[N], u[N], zw[N], l0[N], l1[N];
    vec4d scale, rr, Rh, Rl, h0, l0d, jv, ed, bad;
    vec4d t, eb, ixd;
    int idx, kk, M;

    /* the exponent and the level-1 index from the bits */
    /* (the lanes of t, eb and the integer ix are uint64_t bit patterns) */
    t = vec4d_bits_sub(xc[0], vec4d_bits_set_u64(LOG_M0_BITS));
    eb = vec4d_bits_shr(vec4d_bits_add(t, vec4d_bits_set_u64(UINT64_C(1023) << 52)), 52);   /* e + 1023 */
    ixd = DF(_v_epi64_to_pd)(vec4d_bit_and(vec4d_bits_shr(t, 45), vec4d_bits_set_u64(127)));
    scale = vec4d_bits_shl(vec4d_bits_sub(vec4d_bits_set_u64(2046), eb), 52);   /* 2^-e */
    ed = vec4d_sub(DF(_v_epi64_to_pd)(eb), vec4d_set_d(1023.0));
    UNROLL
    for (idx = 0; idx < N; idx++)
        m[idx] = vec4d_mul(xc[idx], scale);
    {
        rr = vec4d_gather(_dfloat_log_R1, vec4d_to_index(ixd));
        DF(_v_gather)(c, &_dfloat_log_L1[0][0], ixd);    /* L (reused below as a term) */
    }
    /* the level-2 index from the double approximation of z, then
       z2 = m (r r') - 1 = (h0 - 1) + l0 + (m1 + ...) Rh + m Rl as in
       the scalar kernel */
    jv = vec4d_round(vec4d_mul(vec4d_fmsub_fused(m[0], rr, vec4d_set_d(1.0)), vec4d_set_d(8192.0)));
    jv = vec4d_max(vec4d_min(jv, vec4d_set_d(64.0)), vec4d_set_d(-64.0));
    {
        vec4d jd = vec4d_add(jv, vec4d_set_d(64.0));
        vec4d rp = vec4d_gather(_dfloat_log_R2, vec4d_to_index(jd));
        DF(_v_gather)(q, &_dfloat_log_L2[0][0], jd);     /* L' */
        Rh = vec4d_mul(rr, rp);
        Rl = vec4d_fmsub_fused(rr, rp, Rh);
    }
    h0 = vec4d_mul(m[0], Rh);
    l0d = vec4d_fmsub_fused(m[0], Rh, h0);
    UNROLL
    for (idx = 0; idx + 1 < N; idx++)
        P[idx] = m[idx + 1];
    P[N - 1] = vec4d_zero();
    DF(_v_mul_d_kernel)(u, NULL, P, Rh, 0, 0);
    DF(_v_mul_d_kernel)(z, NULL, m, Rl, 0, 0);
    D[0] = vec4d_sub(h0, vec4d_set_d(1.0));
    eL[0] = l0d;
    UNROLL
    for (idx = 1; idx < N; idx++)
        D[idx] = eL[idx] = vec4d_zero();
    {
        const vec4d * xs[4] = {D, eL, u, z};
        DF(_v_addn_kernel)(z2, NULL, xs, 4, 0, 0);
    }
    bad = DF(_v_canonicalise)(z2);
    bad = vec4d_bit_or(bad, vec4d_cmp_nge(vec4d_abs(z2[0]), vec4d_set_d(LOG_ZMIN)));
    /* D = e L0 + e (L1 + ...) + L + L' + z2 */
    UNROLL
    for (idx = 0; idx < N; idx++)
        r[idx] = vec4d_set_d(_dfloat_exp_ln2_42[idx + 1]);
    DF(_v_mul_d_kernel)(eL, NULL, r, ed, 0, 0);
    u[0] = vec4d_mul(ed, vec4d_set_d(_dfloat_exp_ln2_42[0]));
    UNROLL
    for (idx = 1; idx < N; idx++)
        u[idx] = vec4d_zero();
    {
        const vec4d * xs[5] = {u, eL, c, q, z2};
        DF(_v_addn_kernel)(D, NULL, xs, 5, 0, 0);
    }
    /* the chains E and O in w */
    DF(_v_sqr_kernel)(w, NULL, z2, 0, 0);
    M = DF(_log_level)[LOG_J2];
    for (idx = 0; idx < N; idx++)
    {
        q[idx] = vec4d_set_d((idx < M && LOG_C(2 * LOG_J2 + 2) != NULL) ? LOG_C(2 * LOG_J2 + 2)[idx] : 0.0);
        r[idx] = vec4d_set_d((idx < M && LOG_C(2 * LOG_J2 + 3) != NULL) ? LOG_C(2 * LOG_J2 + 3)[idx] : 0.0);
    }
    for (kk = LOG_J2 - 1; kk >= 0; kk--)
    {
        M = DF(_log_level)[kk];
        for (idx = 0; idx < M; idx++)
            c[idx] = vec4d_set_d((LOG_C(2 * kk + 2) != NULL) ? LOG_C(2 * kk + 2)[idx] : 0.0);
        DF(_vexp_mul_level)(u, w, q, M);
        DF(_vexp_add_level)(q, c, u, M);
        for (idx = 0; idx < M; idx++)
            c[idx] = vec4d_set_d((LOG_C(2 * kk + 3) != NULL) ? LOG_C(2 * kk + 3)[idx] : 0.0);
        DF(_vexp_mul_level)(u, w, r, M);
        DF(_vexp_add_level)(r, c, u, M);
    }
    DF(_v_mul_kernel)(zw, NULL, z2, w, 0, 0);
    DF(_v_mul_kernel)(l0, NULL, w, q, 0, 0);
    DF(_v_mul_kernel)(l1, NULL, zw, r, 0, 0);
    DF(_v_add3_kernel)(res, NULL, D, l0, l1, 0, 0);
    return bad;
}

/* the arctangent kernel on 4 lanes (N >= 2), for lanes with
   ATAN_XMIN <= |x0| < ATAN_XMAX; returns the mask of lanes to redo */
static vec4d
DF(_v_atan_core)(vec4d * res, const vec4d * xc)
{
    vec4d x[N], n[N], d[N], t[N], z[N], w[N], v[N], C[N], D[N], E0[N], E1[N], u[N], c[N], p[N], l0[N], l1[N];
    vec4d sg, big, iv, cv, i0, bad;
    int idx, kk, M;

    sg = vec4d_bit_and(xc[0], vec4d_set_d(-0.0));
    UNROLL
    for (idx = 0; idx < N; idx++)
        x[idx] = vec4d_bit_xor(xc[idx], sg);
    big = vec4d_cmp_gt(x[0], vec4d_set_d(1.0));
    iv = vec4d_blendv(vec4d_mul(x[0], vec4d_set_d(64.0)), vec4d_div(vec4d_set_d(64.0), x[0]), big);
    iv = vec4d_round(iv);
    iv = vec4d_max(vec4d_min(iv, vec4d_set_d(64.0)), vec4d_zero());
    cv = vec4d_mul(iv, vec4d_set_d(1.0 / 64.0));
    i0 = vec4d_bit_andnot(big, vec4d_cmp_eq(iv, vec4d_zero()));   /* x <= 1/128: z = x */
    /* t = c x; n = x - c or 1 - c x; d = 1 + c x or x + c */
    DF(_v_mul_d_kernel)(t, NULL, x, cv, 0, 0);
    n[0] = vec4d_blendv(vec4d_sub(x[0], cv), vec4d_sub(vec4d_set_d(1.0), t[0]), big);
    UNROLL
    for (idx = 1; idx < N; idx++)
        n[idx] = vec4d_blendv(x[idx], vec4d_neg(t[idx]), big);
    /* the head has cancelled: renormalise (for N > 2 the
       canonicalisation's TwoSum chain does it) */
    if (N == 2)
        V_TWO_SUM(n[0], n[1], n[0], n[1]);
    bad = DF(_v_canonicalise)(n);
    u[0] = vec4d_set_d(1.0);
    UNROLL
    for (idx = 1; idx < N; idx++)
        u[idx] = vec4d_zero();
    DF(_v_addsub_kernel)(d, NULL, u, t, 0, 0, 0);
    c[0] = cv;
    UNROLL
    for (idx = 1; idx < N; idx++)
        c[idx] = vec4d_zero();
    DF(_v_addsub_kernel)(p, NULL, x, c, 0, 0, 0);
    UNROLL
    for (idx = 0; idx < N; idx++)
        d[idx] = vec4d_blendv(d[idx], p[idx], big);
    DF(_v_div_level)(z, n, d, N, N);
    UNROLL
    for (idx = 0; idx < N; idx++)
    {
        z[idx] = vec4d_blendv(z[idx], x[idx], i0);
        z[idx] = vec4d_bit_xor(z[idx], vec4d_bit_and(big, vec4d_set_d(-0.0)));
    }
    bad = vec4d_bit_or(bad, DF(_v_canonicalise)(z));
    /* C = atan c, or pi/2 - atan c */
    DF(_v_gather)(t, &_dfloat_atan_C[0][0], iv);
    UNROLL
    for (idx = 0; idx < N; idx++)
        c[idx] = vec4d_set_d(_dfloat_pi_2[0][idx]);
    DF(_v_addsub_kernel)(p, NULL, c, t, 1, 0, 0);
    UNROLL
    for (idx = 0; idx < N; idx++)
        C[idx] = vec4d_blendv(t[idx], p[idx], big);
    DF(_v_addsub_kernel)(D, NULL, C, z, 0, 0, 0);
    /* the chains E0, E1 in v = w^2 */
    DF(_v_sqr_kernel)(w, NULL, z, 0, 0);
    DF(_v_sqr_kernel)(v, NULL, w, 0, 0);
    M = DF(_atan_level)[ATAN_J2];
    for (idx = 0; idx < N; idx++)
    {
        E0[idx] = vec4d_set_d((idx < M && ATAN_C(2 * ATAN_J2) != NULL) ? ATAN_C(2 * ATAN_J2)[idx] : 0.0);
        E1[idx] = vec4d_set_d((idx < M && ATAN_C(2 * ATAN_J2 + 1) != NULL) ? ATAN_C(2 * ATAN_J2 + 1)[idx] : 0.0);
    }
    for (kk = ATAN_J2 - 1; kk >= 0; kk--)
    {
        M = DF(_atan_level)[kk];
        for (idx = 0; idx < M; idx++)
            c[idx] = vec4d_set_d((ATAN_C(2 * kk) != NULL) ? ATAN_C(2 * kk)[idx] : 0.0);
        DF(_vexp_mul_level)(u, v, E0, M);
        DF(_vexp_add_level)(E0, c, u, M);
        for (idx = 0; idx < M; idx++)
            c[idx] = vec4d_set_d((ATAN_C(2 * kk + 1) != NULL) ? ATAN_C(2 * kk + 1)[idx] : 0.0);
        DF(_vexp_mul_level)(u, v, E1, M);
        DF(_vexp_add_level)(E1, c, u, M);
    }
    DF(_v_mul_kernel)(p, NULL, z, w, 0, 0);
    DF(_v_mul_kernel)(l0, NULL, p, E0, 0, 0);
    DF(_v_mul_kernel)(t, NULL, p, w, 0, 0);
    DF(_v_mul_kernel)(l1, NULL, t, E1, 0, 0);
    DF(_v_add3_kernel)(res, NULL, D, l0, l1, 0, 0);
    UNROLL
    for (idx = 0; idx < N; idx++)
        res[idx] = vec4d_bit_xor(res[idx], sg);
    return bad;
}

/* the lanes the vector kernels do not handle */
FLINT_FORCE_INLINE vec4d
DF(_v_log_bad)(vec4d * xc)
{
    vec4d bad = DF(_v_canonicalise)(xc);
    bad = vec4d_bit_or(bad, vec4d_cmp_nge(xc[0], vec4d_set_d(0x1p-1022)));
    bad = vec4d_bit_or(bad, vec4d_cmp_nlt(xc[0], vec4d_set_d(0x1p1022)));
    return bad;
}

FLINT_FORCE_INLINE vec4d
DF(_v_atan_bad)(vec4d * xc)
{
    vec4d bad = DF(_v_canonicalise)(xc);
    vec4d ax = vec4d_abs(xc[0]);
    bad = vec4d_bit_or(bad, vec4d_cmp_nge(ax, vec4d_set_d(ATAN_XMIN)));
    bad = vec4d_bit_or(bad, vec4d_cmp_nlt(ax, vec4d_set_d(ATAN_XMAX)));
    return bad;
}
#endif
#endif

void DFV(vec_log)(DF(ptr) res, DF(srcptr) x, slong len)
{
    slong i = 0;
#if defined(DFLOAT_HAVE_SIMD) && N >= 2
    for (; i + 4 <= len; i += 4)
    {
        DF(struct) xs[4];
        vec4d xc[N], y[N], bad;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N);
        bad = DF(_v_log_bad)(xc);
        bad = vec4d_bit_or(bad, DF(_v_log_core)(y, xc));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DF(_vstore4)(res[i].d, y, N);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DF(log)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DF(log)(res + i, x + i);
}

void DFV(vec_atan)(DF(ptr) res, DF(srcptr) x, slong len)
{
    slong i = 0;
#if defined(DFLOAT_HAVE_SIMD) && N >= 2
    for (; i + 4 <= len; i += 4)
    {
        DF(struct) xs[4];
        vec4d xc[N], y[N], bad;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N);
        bad = DF(_v_atan_bad)(xc);
        bad = vec4d_bit_or(bad, DF(_v_atan_core)(y, xc));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DF(_vstore4)(res[i].d, y, N);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DF(atan)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DF(atan)(res + i, x + i);
}

int DBV(vec_log)(DB(ptr) res, DB(srcptr) x, slong len)
{
    slong i = 0;
    int status = GR_SUCCESS;
#if defined(DFLOAT_HAVE_SIMD) && N >= 2
    for (; i + 4 <= len; i += 4)
    {
        DB(struct) xs[4];
        vec4d xc[N], y[N], bad, rx, mlow, rad, ay;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N + 1);
        DB(_vload_rad)(&rx, x + i);
        bad = DF(_v_log_bad)(xc);
        /* a lower bound on the ball: the head less its tail (2u of
           it) less the radius, with slack for the roundings */
        mlow = vec4d_sub(vec4d_sub(xc[0], vec4d_mul(xc[0], vec4d_set_d(0x1p-52))), rx);
        mlow = vec4d_mul(mlow, vec4d_set_d(1.0 - 0x1p-50));
        bad = vec4d_bit_or(bad, vec4d_cmp_ngt(mlow, vec4d_set_d(0x1p-1000)));
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(rx, vec4d_set_d(DBL_MAX)));
        bad = vec4d_bit_or(bad, DF(_v_log_core)(y, xc));
        ay = DF(_vabs_sum)(y);
        rad = vec4d_add(vec4d_mul(ay, vec4d_set_d(DFLOAT_LOG_EPS(N))), vec4d_div(rx, mlow));
        rad = vec4d_mul(rad, vec4d_set_d(DFLOAT_RAD_SLACK));
        rad = vec4d_blendv(rad, vec4d_add(rad, vec4d_set_d(0x1p-1069)),
            vec4d_cmp_lt(rad, vec4d_set_d(DFLOAT_RAD_UNDERFLOW_THRESHOLD)));
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(rad, vec4d_set_d(DBL_MAX)));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DB(_vstore4)(res + i, y, rad);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    status |= DB(log)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        status |= DB(log)(res + i, x + i);
    return status;
}

void DBV(vec_atan)(DB(ptr) res, DB(srcptr) x, slong len)
{
    slong i = 0;
#if defined(DFLOAT_HAVE_SIMD) && N >= 2
    for (; i + 4 <= len; i += 4)
    {
        DB(struct) xs[4];
        vec4d xc[N], y[N], bad, rx, rad, ay;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N + 1);
        DB(_vload_rad)(&rx, x + i);
        bad = DF(_v_atan_bad)(xc);
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(rx, vec4d_set_d(DBL_MAX)));
        bad = vec4d_bit_or(bad, DF(_v_atan_core)(y, xc));
        ay = DF(_vabs_sum)(y);
        rad = vec4d_add(vec4d_mul(ay, vec4d_set_d(DFLOAT_ATAN_EPS(N))), rx);
        rad = vec4d_mul(rad, vec4d_set_d(DFLOAT_RAD_SLACK));
        rad = vec4d_blendv(rad, vec4d_add(rad, vec4d_set_d(0x1p-1069)),
            vec4d_cmp_lt(rad, vec4d_set_d(DFLOAT_RAD_UNDERFLOW_THRESHOLD)));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DB(_vstore4)(res + i, y, rad);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DB(atan)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DB(atan)(res + i, x + i);
}

/* a ball from a 5-term table row (the first N components, the rest
   and the table residual in the radius) */
FLINT_FORCE_INLINE void
DB(_set_table_row)(DB(t) res, const double * row)
{
    double r = _dfloat_exp_tab_err;
    int i;
    for (i = 0; i < N; i++)
        res->d[i] = row[i];
    for (i = N; i < 5; i++)
        r += fabs(row[i]);
    res->rad = _dfloat_rad_finish(r);
}

/* ---- compositions: tan, sinh, cosh, tanh, asin, acos, atan2, pow ----

   Built from the kernels above with the plain arithmetic; the error
   bounds follow from those of the parts (see the documentation).  The
   ball versions are the same compositions on balls. */

DFLOAT_FUNC_NOINLINE void DF(tan)(DF(t) res, DF(srcptr) x_)
{
    DF(t) x;
    DF(_canonicalise_inline)(x, x_);
    DF(t) s, c;
    DF(sin_cos)(s, c, x);
    DF(_div_inline)(res, s, c);
}

DFLOAT_FUNC_NOINLINE int DB(tan)(DB(t) res, DB(srcptr) x_)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x_)))
        return GR_UNABLE;
    DB(t) x;
    DB(_canonicalise_inline)(x, x_);
    DB(t) s, c;
    DB(sin_cos)(s, c, x);
    if (FLINT_UNLIKELY(DB(div)(res, s, c) != GR_SUCCESS))
    {
        /* near a pole the cosine's absolute error covers zero: arb
           decides with the relative accuracy it needs (the poles are
           irrational, so a ball never lies outside the domain) */
        arb_t t;
        int status;
        arb_init(t);
        DB(get_arb)(t, x);
        arb_tan(t, t, 53 * N + 60);
        status = arb_is_finite(t) ? GR_SUCCESS : GR_UNABLE;
        if (status == GR_SUCCESS)
            DB(set_arb)(res, t);
        arb_clear(t);
        return status;
    }
    return GR_SUCCESS;
}

/* sinh |x| = expm1 (expm1 + 2) / (2 (expm1 + 1)) for |x| < 1 (relative
   near 0) and (E - 1/E)/2 beyond (E >= e, no cancellation); the
   operands of the arithmetic never alias */
DFLOAT_FUNC_NOINLINE void DF(sinh)(DF(t) res, DF(srcptr) x_)
{
    DF(t) x;
    DF(_canonicalise_inline)(x, x_);
    DF(t) a, e, t, u, v, w;
    DF(_abs_inline)(a, x);
    if (a->d[0] < 1.0)
    {
        DF(expm1)(e, a);
        DF(_set_d_inline)(t, 2.0);
        DF(_add_impl)(u, e, t);
        DF(_mul_impl)(v, e, u);
        DF(_set_d_inline)(t, 1.0);
        DF(_add_impl)(u, e, t);
        DF(_mul_2exp_si_impl)(w, u, 1);
        DF(_div_inline)(res, v, w);
    }
    else
    {
        DF(exp)(e, a);
        DF(_inv_inline)(t, e);
        DF(_sub_impl)(u, e, t);
        DF(_mul_2exp_si_impl)(res, u, -1);
    }
    if (x->d[0] < 0.0)
        DF(_neg_inline)(res, res);
}

DFLOAT_FUNC_NOINLINE void DB(sinh)(DB(t) res, DB(srcptr) x_)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x_)))
    {
        DB(_indeterminate_inline)(res);
        return;
    }
    DB(t) x;
    DB(_canonicalise_inline)(x, x_);
    DB(t) a, e, t, u, v, w;
    int neg = (x->d[0] < 0.0);
    DB(_abs_inline)(a, x);
    if (a->d[0] < 1.0)
    {
        DB(expm1)(e, a);
        DB(_set_d_inline)(t, 2.0);
        DB(_add_impl)(u, e, t);
        DB(_mul_impl)(v, e, u);
        DB(_set_d_inline)(t, 1.0);
        DB(_add_impl)(u, e, t);
        DB(_mul_2exp_si_impl)(w, u, 1);
        if (DB(div)(res, v, w) != GR_SUCCESS)
            DB(_indeterminate_inline)(res);
    }
    else
    {
        DB(exp)(e, a);
        if (DB(inv)(t, e) != GR_SUCCESS)
            DB(_indeterminate_inline)(res);
        else
        {
            DB(_sub_impl)(u, e, t);
            DB(_mul_2exp_si_impl)(res, u, -1);
        }
    }
    if (neg)
        DB(_neg_inline)(res, res);
}

/* cosh x = (E + 1/E)/2 with E = exp |x| */
DFLOAT_FUNC_NOINLINE void DF(cosh)(DF(t) res, DF(srcptr) x_)
{
    DF(t) x;
    DF(_canonicalise_inline)(x, x_);
    DF(t) a, e, t, u;
    DF(_abs_inline)(a, x);
    DF(exp)(e, a);
    DF(_inv_inline)(t, e);
    DF(_add_impl)(u, e, t);
    DF(_mul_2exp_si_impl)(res, u, -1);
}

DFLOAT_FUNC_NOINLINE void DB(cosh)(DB(t) res, DB(srcptr) x_)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x_)))
    {
        DB(_indeterminate_inline)(res);
        return;
    }
    DB(t) x;
    DB(_canonicalise_inline)(x, x_);
    DB(t) a, e, t, u;
    DB(_abs_inline)(a, x);
    DB(exp)(e, a);
    if (DB(inv)(t, e) != GR_SUCCESS)
    {
        DB(_indeterminate_inline)(res);
        return;
    }
    DB(_add_impl)(u, e, t);
    DB(_mul_2exp_si_impl)(res, u, -1);
}

/* tanh |x| = expm1(2|x|) / (expm1(2|x|) + 2) for |x| < 1 (relative
   near 0) and 1 - 2/(E + 1), E = exp(2|x|), beyond (tanh >= 0.76) */
DFLOAT_FUNC_NOINLINE void DF(tanh)(DF(t) res, DF(srcptr) x_)
{
    DF(t) x;
    DF(_canonicalise_inline)(x, x_);
    DF(t) a, b, e, t, u, v;
    DF(_abs_inline)(a, x);
    DF(_mul_2exp_si_impl)(b, a, 1);
    if (b->d[0] < 2.0)
    {
        DF(expm1)(e, b);
        DF(_set_d_inline)(t, 2.0);
        DF(_add_impl)(u, e, t);
        DF(_div_inline)(res, e, u);
    }
    else
    {
        DF(exp)(e, b);
        DF(_set_d_inline)(t, 1.0);
        DF(_add_impl)(u, e, t);
        DF(_set_d_inline)(t, 2.0);
        DF(_div_inline)(v, t, u);
        DF(_set_d_inline)(t, 1.0);
        DF(_sub_impl)(res, t, v);
    }
    if (x->d[0] < 0.0)
        DF(_neg_inline)(res, res);
}

DFLOAT_FUNC_NOINLINE void DB(tanh)(DB(t) res, DB(srcptr) x_)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x_)))
    {
        DB(_set_zero_pm)(res, 1.0);
        return;
    }
    DB(t) x;
    DB(_canonicalise_inline)(x, x_);
    DB(t) a, b, e, t, u, v;
    int neg = (x->d[0] < 0.0), status;
    DB(_abs_inline)(a, x);
    DB(_mul_2exp_si_impl)(b, a, 1);
    if (b->d[0] < 2.0)
    {
        DB(expm1)(e, b);
        DB(_set_d_inline)(t, 2.0);
        DB(_add_impl)(u, e, t);
        status = DB(div)(res, e, u);
    }
    else
    {
        DB(exp)(e, b);
        DB(_set_d_inline)(t, 1.0);
        DB(_add_impl)(u, e, t);
        DB(_set_d_inline)(t, 2.0);
        status = DB(div)(v, t, u);
        DB(_set_d_inline)(t, 1.0);
        DB(_sub_impl)(res, t, v);
    }
    if (status != GR_SUCCESS)
        DB(_indeterminate_inline)(res);
    if (neg)
        DB(_neg_inline)(res, res);
}

/* asin x = atan(x / sqrt((1 - x)(1 + x))) (1 - x and 1 + x are exact
   where they cancel), acos x = 2 atan(sqrt((1 - x) / (1 + x))) */
DFLOAT_FUNC_NOINLINE void DF(asin)(DF(t) res, DF(srcptr) x_)
{
    DF(t) x;
    DF(_canonicalise_inline)(x, x_);
    DF(t) one, a, b, t, u, v;
    DF(_set_d_inline)(one, 1.0);
    DF(_sub_impl)(a, one, x);
    DF(_add_impl)(b, one, x);
    DF(_mul_impl)(t, a, b);
    DF(_sqrt_inline)(u, t);
    DF(_div_inline)(v, x, u);
    DF(atan)(res, v);
}

/* whether y is an integer: GR_SUCCESS with e if y is an exact integer,
   GR_DOMAIN if the ball certainly contains no integer, else GR_UNABLE
   (the whole line). A quick test in doubles decides most balls; the
   rest (an integer near an end, or huge values) are decided exactly. */
static int
DB(_get_fmpz)(fmpz_t e, DB(srcptr) y)
{
    arf_t t, u, r;
    int status;
    if (DB(_is_W)(y))
        return GR_UNABLE;
    if (y->rad != 0.0)
    {
        /* the ends are in [a, a2] and [b2, b] up to the margin */
        double T = _dfloat_rad_finish(DF(_tail_abs_sum)(y->d));
        double s = _dfloat_rad_finish(T + y->rad);
        double a = y->d[0] - s, b = y->d[0] + s, a2 = y->d[0] + T - y->rad, b2 = y->d[0] - T + y->rad, k, margin;
        margin = FLINT_MAX(fabs(a), fabs(b)) * 0x1p-52 + DBL_MIN;
        if (fabs(a) < 0x1p51 && fabs(b) < 0x1p51)
        {
            k = ceil(a - margin);       /* the least integer that may be in */
            if (k > b + margin)
                return GR_DOMAIN;
            if (k >= a2 + margin && k <= b2 - margin)
                return GR_UNABLE;       /* k is certainly inside */
        }
        /* exact: ceil(mid - rad) > mid + rad */
        arf_init(t); arf_init(u); arf_init(r);
        _dfloat_get_arf(t, y->d, N);
        arf_set_d(r, y->rad);
        arf_add(u, t, r, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_sub(t, t, r, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_ceil(t, t);
        status = (arf_cmp(t, u) > 0) ? GR_DOMAIN : GR_UNABLE;
        arf_clear(t); arf_clear(u); arf_clear(r);
        return status;
    }
    arf_init(t);
    _dfloat_get_arf(t, y->d, N);
    if (arf_is_int(t))
    {
        arf_get_fmpz(e, t, ARF_RND_DOWN);
        status = GR_SUCCESS;
    }
    else
        status = GR_DOMAIN;
    arf_clear(t);
    return status;
}

/* asin and acos where the formulas fail (at and near +-1, where 1 - x^2
   contains zero): for a ball inside [-1, 1], the values at the exact
   endpoints (the functions are monotone); GR_DOMAIN for a ball outside
   [-1, 1], GR_UNABLE otherwise */
static int
DB(_asin_acos_fallback)(DB(t) res, DB(srcptr) x, int is_acos)
{
    arf_t m, r, lo, hi;
    arb_t a, b;
    slong prec = 53 * N + 30;
    if (!(DB_ALL_GE(x, -1.0) && DB_ALL_LE(x, 1.0)))
        return (DB_ALL_LT(x, -1.0) || DB_ALL_GT(x, 1.0)) ? GR_DOMAIN : GR_UNABLE;
    arf_init(m); arf_init(r); arf_init(lo); arf_init(hi);
    arb_init(a); arb_init(b);
    _dfloat_get_arf(m, x->d, N);
    arf_set_d(r, x->rad);
    arf_sub(lo, m, r, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_add(hi, m, r, ARF_PREC_EXACT, ARF_RND_DOWN);
    /* (the comparisons above have slack: clamp the exact endpoints) */
    if (arf_cmp_si(lo, -1) < 0) arf_set_si(lo, -1);
    if (arf_cmp_si(hi, 1) > 0) arf_set_si(hi, 1);
    arb_set_arf(a, lo);
    arb_set_arf(b, hi);
    if (is_acos)
    {
        arb_acos(a, a, prec);
        arb_acos(b, b, prec);
    }
    else
    {
        arb_asin(a, a, prec);
        arb_asin(b, b, prec);
    }
    arb_union(a, a, b, prec);
    DB(set_arb)(res, a);
    arf_clear(m); arf_clear(r); arf_clear(lo); arf_clear(hi);
    arb_clear(a); arb_clear(b);
    return GR_SUCCESS;
}

DFLOAT_FUNC_NOINLINE int DB(asin)(DB(t) res, DB(srcptr) x_)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x_)))
        return GR_UNABLE;
    DB(t) x;
    DB(_canonicalise_inline)(x, x_);
    DB(t) one, a, b, t, u, v;
    int status;
    DB(_set_d_inline)(one, 1.0);
    DB(_sub_impl)(a, one, x);
    DB(_add_impl)(b, one, x);
    DB(_mul_impl)(t, a, b);
    status = DB(sqrt)(u, t);
    if (status == GR_SUCCESS)
        status = DB(div)(v, x, u);
    if (status != GR_SUCCESS)
        return DB(_asin_acos_fallback)(res, x, 0);
    DB(atan)(res, v);
    return GR_SUCCESS;
}

DFLOAT_FUNC_NOINLINE void DF(acos)(DF(t) res, DF(srcptr) x_)
{
    DF(t) x;
    DF(_canonicalise_inline)(x, x_);
    DF(t) one, a, b, t, u, v;
    DF(_set_d_inline)(one, 1.0);
    DF(_sub_impl)(a, one, x);
    DF(_add_impl)(b, one, x);
    DF(_div_inline)(t, a, b);
    DF(_sqrt_inline)(u, t);
    DF(atan)(v, u);
    DF(_mul_2exp_si_impl)(res, v, 1);
}

DFLOAT_FUNC_NOINLINE int DB(acos)(DB(t) res, DB(srcptr) x_)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x_)))
        return GR_UNABLE;
    DB(t) x;
    DB(_canonicalise_inline)(x, x_);
    DB(t) one, a, b, t, u, v;
    int status;
    DB(_set_d_inline)(one, 1.0);
    DB(_sub_impl)(a, one, x);
    DB(_add_impl)(b, one, x);
    status = DB(div)(t, a, b);
    if (status == GR_SUCCESS)
        status = DB(sqrt)(u, t);
    if (status != GR_SUCCESS)
        return DB(_asin_acos_fallback)(res, x, 1);
    DB(atan)(v, u);
    DB(_mul_2exp_si_impl)(res, v, 1);
    return GR_SUCCESS;
}

/* atan2(y, x) = atan(y/x) for x > 0, atan(y/x) +- pi for x < 0, +- pi/2
   for x = 0 (0 for y = 0); nonfinite arguments follow the double
   atan2 (the multiples of pi/4 it returns are made exact) */
DFLOAT_FUNC_NOINLINE void DF(atan2)(DF(t) res, DF(srcptr) y_, DF(srcptr) x_)
{
    DF(t) x, y;
    DF(_canonicalise_inline)(x, x_);
    DF(_canonicalise_inline)(y, y_);
    DF(t) q, p, r;
    double x0 = x->d[0], y0 = y->d[0];
    if (!(fabs(x0) <= DBL_MAX) || !(fabs(y0) <= DBL_MAX) || (x0 == 0.0 && DF(_is_zero_inline)(x)))
    {
        double v = atan2(y0, x0);
        if (v != v)
            DF(_set_d_inline)(res, v);
        else
        {
            DF(_set_inline)((DF(ptr)) p, (DF(srcptr)) _dfloat_pi_2[0]);
            DF(_mul_d_impl)(res, p, v / 1.5707963267948966);   /* an exact multiple of 1/2 */
        }
        return;
    }
    DF(_div_inline)(q, y, x);
    DF(atan)(r, q);
    if (x0 < 0.0)
    {
        DF(_set_inline)((DF(ptr)) q, (DF(srcptr)) _dfloat_pi_2[0]);
        DF(_mul_2exp_si_impl)(p, q, 1);
        if (y0 < 0.0 || (y0 == 0.0 && DF(sgn)(y) < 0))
            DF(_sub_impl)(res, r, p);
        else
            DF(_add_impl)(res, r, p);
    }
    else
        DF(_set_inline)(res, r);
}

/* the ball version always returns a finite enclosure: for x away from
   zero, atan(y/x) (+- pi for x < 0, or [0 +/- pi] if y also contains
   zero: the branch cut); for x containing zero and y away from it,
   +-pi/2 - atan(x/y) (an identity valid for all x); for both
   containing zero, [0 +/- pi] (0 exactly for the exact zeros) */
DFLOAT_FUNC_NOINLINE void DB(atan2)(DB(t) res, DB(srcptr) y_, DB(srcptr) x_)
{
    if (FLINT_UNLIKELY(DB(_is_W)(x_) || DB(_is_W)(y_)))
    {
        DB(_set_zero_pm)(res, 0x1.921fb54442d19p+1);   /* pi rounded up */
        return;
    }
    DB(t) x, y, q, p, r;
    int xz, yz;
    DB(_canonicalise_inline)(x, x_);
    DB(_canonicalise_inline)(y, y_);
    xz = DB(_contains_zero_inline)(x);
    yz = DB(_contains_zero_inline)(y);
    if (!xz)
    {
        if (DB(div)(q, y, x) != GR_SUCCESS)
        {
            DB(_indeterminate_inline)(res);
            return;
        }
        DB(atan)(r, q);
        if (x->d[0] > 0.0)
            DB(_set_inline)(res, r);
        else if (DB(_is_zero_inline)(y))
        {
            /* atan2(0, x) = pi for x < 0 */
            DB(_set_table_row)(q, _dfloat_pi_2[0]);
            DB(_mul_2exp_si_impl)(res, q, 1);
        }
        else if (yz)
        {
            DB(_zero_inline)(res);
            res->rad = 0x1.921fb54442d19p+1;   /* pi rounded up */
        }
        else
        {
            DB(_set_table_row)(q, _dfloat_pi_2[0]);
            DB(_mul_2exp_si_impl)(p, q, 1);
            if (y->d[0] < 0.0)
                DB(_sub_impl)(res, r, p);
            else
                DB(_add_impl)(res, r, p);
        }
    }
    else if (!yz)
    {
        if (DB(div)(q, x, y) != GR_SUCCESS)
        {
            DB(_indeterminate_inline)(res);
            return;
        }
        DB(atan)(r, q);
        DB(_set_table_row)(p, _dfloat_pi_2[0]);
        if (y->d[0] < 0.0)
            DB(_neg_inline)(p, p);
        DB(_sub_impl)(res, p, r);
    }
    else if (DB(_is_zero_inline)(x) && DB(_is_zero_inline)(y))
        DB(_zero_inline)(res);
    else
    {
        DB(_zero_inline)(res);
        res->rad = 0x1.921fb54442d19p+1;
    }
}

/* x^e for an integer e with |e| < 2^40 by binary powering (a relative
   error of about 2 log2 |e| u^N) */
static void
DF(_pow_si)(DF(t) res, DF(srcptr) x, slong e)
{
    DF(t) b, r, t;
    ulong a = (e < 0) ? (ulong) (-e) : (ulong) e;
    /* (1/x)^|e| for e < 0: a tiny x is the more common case */
    if (e < 0)
        DF(_inv_inline)(b, x);
    else
        DF(_set_inline)(b, x);
    DF(_one_inline)(r);
    while (a)
    {
        if (a & 1)
        {
            DF(_mul_impl)(t, r, b);
            DF(_set_inline)(r, t);
        }
        a >>= 1;
        if (a)
        {
            DF(_sqr_impl)(t, b);
            DF(_set_inline)(b, t);
        }
    }
    DF(_set_inline)(res, r);
}

/* pow(x, y): an integer y (|y| < 2^40) by binary powering, otherwise
   exp(y log x) for x > 0 (nan for x < 0; 0^y = 0, 1, inf by the sign
   of y; nonfinite arguments follow the double pow) */
DFLOAT_FUNC_NOINLINE void DF(pow)(DF(t) res, DF(srcptr) x_, DF(srcptr) y_)
{
    DF(t) x, y, l, m;
    double x0, y0;
    DF(_canonicalise_inline)(x, x_);
    DF(_canonicalise_inline)(y, y_);
    x0 = x->d[0];
    y0 = y->d[0];
    if (fabs(y0) < 0x1p40 && y0 == rint(y0) && DF(_tail_abs_sum)(y->d) == 0.0 && fabs(x0) <= DBL_MAX)
    {
        DF(_pow_si)(res, x, (slong) y0);
        return;
    }
    if (!(x0 > 0.0) || !(x0 <= DBL_MAX) || !(fabs(y0) <= DBL_MAX))
    {
        if (x0 == 0.0)
            DF(_set_d_inline)(res, (y0 > 0.0) ? 0.0 : (y0 == 0.0) ? 1.0 : D_INF);
        else
            DF(_set_d_inline)(res, pow(x0, y0));
        return;
    }
    DF(log)(l, x);
    DF(_mul_impl)(m, l, y);
    DF(exp)(res, m);
}

/* the ball version: an exact integer y by binary powering (GR_UNABLE
   for a negative y and a base containing zero), otherwise x must be
   positive */
/* pow on balls, with the real domain of x^y: x > 0; x = 0 and y >= 0;
   x < 0 and y an integer (GR_SUCCESS for balls inside it, GR_DOMAIN for
   balls outside, GR_UNABLE otherwise). An exact integer y (any base) by
   binary powering for |y| < 2^40 and through arb beyond; otherwise
   exp(y log x) for x > 0. */
DFLOAT_FUNC_NOINLINE int DB(pow)(DB(t) res, DB(srcptr) x_, DB(srcptr) y_)
{
    DB(t) x, y, l, m;
    fmpz_t e;
    int status, ystat;
    DB(_canonicalise_inline)(x, x_);     /* (the whole line becomes [0 +/- inf]) */
    DB(_canonicalise_inline)(y, y_);
    fmpz_init(e);
    ystat = DB(_get_fmpz)(e, y);
    if (ystat == GR_SUCCESS)
    {
        if (fmpz_sgn(e) < 0 && DB(_contains_zero_inline)(x))
            status = DB(_is_zero_inline)(x) ? GR_DOMAIN : GR_UNABLE;
        else if (fmpz_bits(e) <= 40)
        {
            DB(t) b, r, t;
            ulong a = FLINT_UABS(fmpz_get_si(e));
            if (fmpz_sgn(e) < 0)
                DB(inv)(b, x);           /* succeeds: x does not contain 0 */
            else
                DB(_set_inline)(b, x);
            DB(_one_inline)(r);                  /* (x^0 = 1 for every x) */
            while (a)
            {
                if (a & 1)
                {
                    DB(_mul_impl)(t, r, b);
                    DB(_set_inline)(r, t);
                }
                a >>= 1;
                if (a)
                {
                    DB(_sqr_impl)(t, b);
                    DB(_set_inline)(b, t);
                }
            }
            DB(_set_inline)(res, r);
            status = GR_SUCCESS;
        }
        else
        {
            arb_t t;
            arb_init(t);
            DB(get_arb)(t, x);
            arb_pow_fmpz(t, t, e, 53 * N + 30);
            DB(set_arb)(res, t);
            arb_clear(t);
            status = GR_SUCCESS;
        }
    }
    else if (DB(_is_W)(x))
        status = GR_UNABLE;          /* contains negative bases and 0 */
    else if (DB_ALL_GT(x, 0.0))
    {
        /* (y may be the whole line: then so is the result, unless
           log x = 0 exactly) */
        status = DB(log)(l, x);
        if (status == GR_SUCCESS)
        {
            DB(_mul_impl)(m, l, y);
            DB(exp)(res, m);
        }
    }
    else if (DB(_is_zero_inline)(x))
    {
        /* 0^y = 0 for y > 0, undefined for y < 0 */
        if (DB_ALL_GT(y, 0.0))
        {
            DB(_zero_inline)(res);
            status = GR_SUCCESS;
        }
        else
            status = DB_ALL_LT(y, 0.0) ? GR_DOMAIN : GR_UNABLE;
    }
    else if (DB_ALL_LT(x, 0.0))
        status = (ystat == GR_DOMAIN) ? GR_DOMAIN : GR_UNABLE;
    else if (DB_ALL_LE(x, 0.0))
        /* x in [lo, 0], lo < 0: out only if y has no integer and y < 0 */
        status = (ystat == GR_DOMAIN && DB_ALL_LT(y, 0.0)) ? GR_DOMAIN : GR_UNABLE;
    else
        status = GR_UNABLE;
    fmpz_clear(e);
    return status;
}

void DFV(vec_exp)(DF(ptr) res, DF(srcptr) x, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    for (; i + 4 <= len; i += 4)
    {
        DF(struct) xs[4];
        vec4d xc[N], y[N], bad, x0;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N);
        bad = DF(_v_canonicalise)(xc);
        x0 = xc[0];
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(x0, vec4d_set_d(709.9)));
        bad = vec4d_bit_or(bad, vec4d_cmp_nge(x0, vec4d_set_d(-746.0)));
        bad = vec4d_bit_or(bad, DF(_v_exp_core)(y, xc));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DF(_vstore4)(res[i].d, y, N);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DF(exp)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DF(exp)(res + i, x + i);
}

void DBV(vec_exp)(DB(ptr) res, DB(srcptr) x, slong len)
{
    slong i = 0;
#ifdef DFLOAT_HAVE_SIMD
    const double eps = DFLOAT_EXP_EPS(N);
    for (; i + 4 <= len; i += 4)
    {
        DB(struct) xs[4];
        vec4d xc[N], y[N], bad, x0, rx, E, ay, rad;
        int mask, l;
        DF(_vload4)(xc, x[i].d, N + 1);
        DB(_vload_rad)(&rx, x + i);
        bad = DF(_v_canonicalise)(xc);
        x0 = xc[0];
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(x0, vec4d_set_d(709.9)));
        bad = vec4d_bit_or(bad, vec4d_cmp_nge(x0, vec4d_set_d(-746.0)));
        bad = vec4d_bit_or(bad, vec4d_cmp_eq(x0, vec4d_zero()));
        /* the input radius: r <= 1 here, the rest goes to scalar */
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(rx, vec4d_set_d(1.0)));
        bad = vec4d_bit_or(bad, DF(_v_exp_core)(y, xc));
        E = vec4d_mul(vec4d_add(rx, vec4d_mul(rx, rx)), vec4d_set_d(DFLOAT_RAD_SLACK));
        ay = vec4d_mul(DF(_vabs_sum)(y), vec4d_set_d(1.0 + 2.0 * eps));
        rad = vec4d_add(vec4d_mul(ay, E), vec4d_mul(ay, vec4d_set_d(2.0 * eps)));
        rad = vec4d_mul(rad, vec4d_set_d(DFLOAT_RAD_SLACK));
        /* the underflow bump of the radius, as in the scalar code;
           nonfinite results go to scalar */
        rad = vec4d_blendv(rad, vec4d_add(rad, vec4d_set_d(0x1p-1069)),
            vec4d_cmp_lt(rad, vec4d_set_d(DFLOAT_RAD_UNDERFLOW_THRESHOLD)));
        bad = vec4d_bit_or(bad, vec4d_cmp_nle(vec4d_add(ay, rad), vec4d_set_d(DBL_MAX)));
        /* the inputs of the lanes redone in scalar, saved before the
           store (the output may alias them) */
        mask = vec4d_movemask(bad);
        if (FLINT_UNLIKELY(mask))
        {
            memcpy(xs, x + i, 4 * sizeof(xs[0]));
        }
        DB(_vstore4)(res + i, y, rad);
        if (FLINT_UNLIKELY(mask))
            for (l = 0; l < 4; l++)
                if (mask & (1 << l))
                    DB(exp)(res + i + l, xs + l);
    }
#endif
    for (; i < len; i++)
        DB(exp)(res + i, x + i);
}

/* ------------------------------------------------------------------ */
/* Generic ring wrappers                                               */
/* ------------------------------------------------------------------ */

#define GRF(name) DFV(gr_##name)
#define GRB(name) DBV(gr_##name)
#define PREC_VIA_ARB (53 * N + 30)

/* ---- constants: table lookups (the first N components of a 5-term
   expansion; the ball adds the rest and the table residual) ---- */

#define DEF_CONST(name, idx) \
    static int GRF(name)(DF(ptr) res, gr_ctx_t ctx) \
    { int i; for (i = 0; i < N; i++) res->d[i] = _dfloat_const[idx][i]; return GR_SUCCESS; } \
    static int GRB(name)(DB(ptr) res, gr_ctx_t ctx) \
    { DB(_set_table_row)(res, _dfloat_const[idx]); return GR_SUCCESS; }

DEF_CONST(pi, DFLOAT_CONST_PI)
DEF_CONST(euler, DFLOAT_CONST_EULER)
DEF_CONST(catalan, DFLOAT_CONST_CATALAN)
static int GRF(expm1)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(expm1)(res, x); return GR_SUCCESS; }
static int GRB(expm1)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(expm1)(res, x); return GR_SUCCESS; }
static int GRF(log1p)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx)
{
    if (!(x->d[0] > -1.0))
    {
        if (x->d[0] == -1.0 && DF(_tail_abs_sum)(x->d) == 0.0) { DF(_set_d_inline)(res, -D_INF); return GR_SUCCESS; }
        if (DF(_is_nan_inline)(x)) { DF(_set_d_inline)(res, D_NAN); return GR_SUCCESS; }
        if (x->d[0] < -1.0) return GR_DOMAIN;
    }
    DF(log1p)(res, x);
    return GR_SUCCESS;
}
static int GRB(log1p)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(log1p)(res, x); }
static int GRF(sin)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(sin)(res, x); return GR_SUCCESS; }
static int GRF(cos)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(cos)(res, x); return GR_SUCCESS; }
static int GRF(sin_cos)(DF(ptr) sn, DF(ptr) cs, DF(srcptr) x, gr_ctx_t ctx) { DF(sin_cos)(sn, cs, x); return GR_SUCCESS; }
static int GRB(sin)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(sin)(res, x); return GR_SUCCESS; }
static int GRB(cos)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(cos)(res, x); return GR_SUCCESS; }
static int GRB(sin_cos)(DB(ptr) sn, DB(ptr) cs, DB(srcptr) x, gr_ctx_t ctx) { DB(sin_cos)(sn, cs, x); return GR_SUCCESS; }
static int GRF(tan)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(tan)(res, x); return GR_SUCCESS; }
static int GRB(tan)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(tan)(res, x); }
static int GRF(sinh)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(sinh)(res, x); return GR_SUCCESS; }
static int GRB(sinh)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(sinh)(res, x); return GR_SUCCESS; }
static int GRF(cosh)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(cosh)(res, x); return GR_SUCCESS; }
static int GRB(cosh)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(cosh)(res, x); return GR_SUCCESS; }
static int GRF(tanh)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(tanh)(res, x); return GR_SUCCESS; }
static int GRB(tanh)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(tanh)(res, x); return GR_SUCCESS; }
/* the domain [-1, 1] of asin and acos: the value of x is |x0| +- the
   tail; a head beyond 1 by more than the tail is out */
static int GRF(asin)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx)
{
    if (fabs(x->d[0]) > 1.0 + DF(_tail_abs_sum)(x->d))
        return GR_DOMAIN;
    DF(asin)(res, x);
    return GR_SUCCESS;
}
static int GRB(asin)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(asin)(res, x); }
static int GRF(acos)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx)
{
    if (fabs(x->d[0]) > 1.0 + DF(_tail_abs_sum)(x->d))
        return GR_DOMAIN;
    DF(acos)(res, x);
    return GR_SUCCESS;
}
static int GRB(acos)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(acos)(res, x); }
static int GRF(atan2)(DF(ptr) res, DF(srcptr) y, DF(srcptr) x, gr_ctx_t ctx) { DF(atan2)(res, y, x); return GR_SUCCESS; }
static int GRB(atan2)(DB(ptr) res, DB(srcptr) y, DB(srcptr) x, gr_ctx_t ctx) { DB(atan2)(res, y, x); return GR_SUCCESS; }
/* pow: an integer exponent by binary powering (any base), otherwise
   exp(y log x) for x > 0 */
static int GRF(get_fmpz)(fmpz_t res, DF(srcptr) x, gr_ctx_t ctx);
static int GRF(pow)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, gr_ctx_t ctx)
{
    fmpz_t e;
    int status;
    fmpz_init(e);
    if (GRF(get_fmpz)(e, y, ctx) == GR_SUCCESS && fmpz_bits(e) <= 40)
        status = gr_pow_fmpz(res, x, e, ctx);
    else if (x->d[0] > 0.0 || (x->d[0] == 0.0 && DF(sgn)(x) >= 0))
    {
        DF(pow)(res, x, y);
        status = GR_SUCCESS;
    }
    else if (DF(_is_nan_inline)(x) || DF(_is_nan_inline)(y))
    {
        DF(_set_d_inline)(res, D_NAN);
        status = GR_SUCCESS;
    }
    else
        status = GR_DOMAIN;
    fmpz_clear(e);
    return status;
}
static int GRB(pow)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx) { return DB(pow)(res, x, y); }

/* ---- plain ring ---- */

static int GRF(ctx_write)(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Floating-point numbers as ");
    status |= gr_stream_write_si(out, N);
    status |= gr_stream_write(out, "-term double expansions (d");
    status |= gr_stream_write_si(out, N);
    status |= gr_stream_write(out, DFLOAT_CTX_STRONG(ctx) ? ", strong)" : ")");
    return status;
}

static int GRF(ctx_get_real_prec)(slong * res, gr_ctx_t ctx) { *res = 53 * N; return GR_SUCCESS; }
static int GRF(ctx_set_real_prec)(gr_ctx_t ctx, slong prec) { return GR_UNABLE; }

static void GRF(init)(DF(ptr) x, gr_ctx_t ctx) { DF(_zero_inline)(x); }
static void GRF(clear)(DF(ptr) x, gr_ctx_t ctx) { }
static void GRF(swap)(DF(ptr) x, DF(ptr) y, gr_ctx_t ctx) { DF(struct) t = *x; *x = *y; *y = t; }

static int GRF(randtest)(DF(ptr) res, flint_rand_t state, gr_ctx_t ctx)
{
    if (n_randint(state, 16) == 0)
        DF(randtest_special)(res, state);
    else
        DF(randtest)(res, state);
    return GR_SUCCESS;
}

static int GRF(write)(gr_stream_t out, DF(srcptr) x, gr_ctx_t ctx)
{
    return gr_stream_write_free(out, DF(get_str)(x, 0));
}

static int GRF(write_n)(gr_stream_t out, DF(srcptr) x, slong n, gr_ctx_t ctx)
{
    return gr_stream_write_free(out, DF(get_str)(x, FLINT_MAX(n, 1)));
}

static int GRF(zero)(DF(ptr) res, gr_ctx_t ctx) { DF(_zero_inline)(res); return GR_SUCCESS; }
static int GRF(one)(DF(ptr) res, gr_ctx_t ctx) { DF(_one_inline)(res); return GR_SUCCESS; }
static int GRF(neg_one)(DF(ptr) res, gr_ctx_t ctx) { DF(_set_d_inline)(res, -1.0); return GR_SUCCESS; }
static int GRF(set)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(_set_inline)(res, x); return GR_SUCCESS; }
static int GRF(set_si)(DF(ptr) res, slong x, gr_ctx_t ctx)
{
    fmpz_t t;
    if (DFLOAT_ULONG_IS_EXACT(FLINT_UABS(x)))
        DF(_set_d_inline)(res, (double) x);
    else
    {
        fmpz_init_set_si(t, x);
        DF(set_fmpz)(res, t);
        fmpz_clear(t);
    }
    return GR_SUCCESS;
}
static int GRF(set_ui)(DF(ptr) res, ulong x, gr_ctx_t ctx)
{
    fmpz_t t;
    if (DFLOAT_ULONG_IS_EXACT(x))
        DF(_set_d_inline)(res, (double) x);
    else
    {
        fmpz_init_set_ui(t, x);
        DF(set_fmpz)(res, t);
        fmpz_clear(t);
    }
    return GR_SUCCESS;
}
static int GRF(set_fmpz)(DF(ptr) res, const fmpz_t x, gr_ctx_t ctx) { DF(set_fmpz)(res, x); return GR_SUCCESS; }
static int GRF(set_fmpq)(DF(ptr) res, const fmpq_t x, gr_ctx_t ctx) { DF(set_fmpq)(res, x); return GR_SUCCESS; }
static int GRF(set_d)(DF(ptr) res, double x, gr_ctx_t ctx) { DF(_set_d_inline)(res, x); return GR_SUCCESS; }
static int GRF(set_str)(DF(ptr) res, const char * s, gr_ctx_t ctx) { return DF(set_str)(res, s); }

static int GRF(set_other)(DF(ptr) res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            DF(set_fmpz)(res, x); return GR_SUCCESS;
        case GR_CTX_FMPQ:
            DF(set_fmpq)(res, x); return GR_SUCCESS;
        case GR_CTX_REAL_FLOAT_ARF:
            DF(set_arf)(res, x); return GR_SUCCESS;
        case GR_CTX_RR_ARB:
            DF(set_arb)(res, x); return GR_SUCCESS;
        case GR_CTX_DFLOAT:
        case GR_CTX_DFLOAT_BALL:
            dfloat_set_dfloat(res->d, N, (const double *) x, DFLOAT_CTX_N(x_ctx));
            return GR_SUCCESS;
        case GR_CTX_DFLOAT_COMPLEX:
        case GR_CTX_DFLOAT_COMPLEX_BALL:
        {
            /* the real part, if the imaginary part (midpoint) is zero */
            int nx = DFLOAT_CTX_N(x_ctx);
            if (!(_dfloat_abs_sum((const double *) x + nx + DFLOAT_CTX_BALL(x_ctx), nx) == 0.0))
                return GR_DOMAIN;
            dfloat_set_dfloat(res->d, N, (const double *) x, nx);
            return GR_SUCCESS;
        }
        default:
        {
            gr_ctx_t actx;
            arb_t t;
            int status;
            gr_ctx_init_real_arb(actx, PREC_VIA_ARB);
            arb_init(t);
            status = gr_set_other(t, x, x_ctx, actx);
            if (status == GR_SUCCESS)
                DF(set_arb)(res, t);
            arb_clear(t);
            gr_ctx_clear(actx);
            return status;
        }
    }
}

static int GRF(get_fmpz)(fmpz_t res, DF(srcptr) x, gr_ctx_t ctx)
{
    arf_t t;
    int status;
    if (!DF(_is_finite_inline)(x))
        return GR_DOMAIN;
    arf_init(t);
    DF(get_arf)(t, x);
    if (arf_is_int(t))
    {
        arf_get_fmpz(res, t, ARF_RND_DOWN);
        status = GR_SUCCESS;
    }
    else
        status = GR_DOMAIN;
    arf_clear(t);
    return status;
}

static int GRF(get_si)(slong * res, DF(srcptr) x, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    status = GRF(get_fmpz)(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_fits_si(t))
            *res = fmpz_get_si(t);
        else
            status = GR_DOMAIN;
    }
    fmpz_clear(t);
    return status;
}

static int GRF(get_ui)(ulong * res, DF(srcptr) x, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    status = GRF(get_fmpz)(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_sgn(t) >= 0 && fmpz_abs_fits_ui(t))
            *res = fmpz_get_ui(t);
        else
            status = GR_DOMAIN;
    }
    fmpz_clear(t);
    return status;
}

static int GRF(get_d)(double * res, DF(srcptr) x, gr_ctx_t ctx) { *res = x->d[0]; return GR_SUCCESS; }

static truth_t GRF(is_zero)(DF(srcptr) x, gr_ctx_t ctx)
{
    if (DF(_is_nan_inline)(x)) return T_UNKNOWN;
    return DF(_is_zero_inline)(x) ? T_TRUE : T_FALSE;
}
static truth_t GRF(is_one)(DF(srcptr) x, gr_ctx_t ctx)
{
    if (DF(_is_nan_inline)(x)) return T_UNKNOWN;
    return DF(_is_one_inline)(x) ? T_TRUE : T_FALSE;
}
static truth_t GRF(is_neg_one)(DF(srcptr) x, gr_ctx_t ctx)
{
    DF(t) t;
    if (DF(_is_nan_inline)(x)) return T_UNKNOWN;
    DF(_neg_inline)(t, x);
    return DF(_is_one_inline)(t) ? T_TRUE : T_FALSE;
}
static truth_t GRF(equal)(DF(srcptr) x, DF(srcptr) y, gr_ctx_t ctx)
{
    if (DF(_is_nan_inline)(x) || DF(_is_nan_inline)(y)) return T_UNKNOWN;
    return DF(_equal_inline)(x, y) ? T_TRUE : T_FALSE;
}

static int GRF(neg)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(_neg_inline)(res, x); return GR_SUCCESS; }
static int GRF(add)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, gr_ctx_t ctx) { DF(_add_impl)(res, x, y); return GR_SUCCESS; }
static int GRF(sub)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, gr_ctx_t ctx) { DF(_sub_impl)(res, x, y); return GR_SUCCESS; }
static int GRF(mul)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, gr_ctx_t ctx) { DF(_mul_impl)(res, x, y); return GR_SUCCESS; }
static int GRF(sqr)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(_sqr_impl)(res, x); return GR_SUCCESS; }
static int GRF(mul_two)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(_mul_2exp_si_impl)(res, x, 1); return GR_SUCCESS; }
static int GRF(mul_2exp_si)(DF(ptr) res, DF(srcptr) x, slong e, gr_ctx_t ctx) { DF(_mul_2exp_si_impl)(res, x, e); return GR_SUCCESS; }
static int GRF(div)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, gr_ctx_t ctx)
{
    if (DF(_is_zero_inline)(y))
        return GR_DOMAIN;
    DF(_div_inline)(res, x, y);
    return GR_SUCCESS;
}
static int GRF(inv)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx)
{
    if (DF(_is_zero_inline)(x))
        return GR_DOMAIN;
    DF(_inv_inline)(res, x);
    return GR_SUCCESS;
}
static int GRF(sqrt)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx)
{
    if (x->d[0] < 0.0)
        return GR_DOMAIN;
    DF(_sqrt_inline)(res, x);
    return GR_SUCCESS;
}
static int GRF(rsqrt)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx)
{
    if (x->d[0] < 0.0 || DF(_is_zero_inline)(x))
        return GR_DOMAIN;
    DF(_rsqrt_inline)(res, x);
    return GR_SUCCESS;
}
static int GRF(abs)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(_abs_inline)(res, x); return GR_SUCCESS; }
static int GRF(sgn)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx)
{
    if (DF(_is_nan_inline)(x)) return GR_UNABLE;
    DF(_set_d_inline)(res, (double) DF(sgn)(x));
    return GR_SUCCESS;
}
static int GRF(cmp)(int * res, DF(srcptr) x, DF(srcptr) y, gr_ctx_t ctx)
{
    if (DF(_is_nan_inline)(x) || DF(_is_nan_inline)(y)) { *res = 0; return GR_UNABLE; }
    *res = DF(_cmp_inline)(x, y);
    return GR_SUCCESS;
}
static int GRF(cmpabs)(int * res, DF(srcptr) x, DF(srcptr) y, gr_ctx_t ctx)
{
    DF(t) a, b;
    DF(_abs_inline)(a, x);
    DF(_abs_inline)(b, y);
    return GRF(cmp)(res, a, b, ctx);
}
static int GRF(pos_inf)(DF(ptr) res, gr_ctx_t ctx) { DF(_set_d_inline)(res, D_INF); return GR_SUCCESS; }
static int GRF(neg_inf)(DF(ptr) res, gr_ctx_t ctx) { DF(_set_d_inline)(res, -D_INF); return GR_SUCCESS; }
static int GRF(nan)(DF(ptr) res, gr_ctx_t ctx) { DF(_set_d_inline)(res, D_NAN); return GR_SUCCESS; }
static int GRF(exp)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(exp)(res, x); return GR_SUCCESS; }
static int GRF(log)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx)
{
    if (!(x->d[0] > 0.0))
    {
        if (DF(_is_zero_inline)(x)) { DF(_set_d_inline)(res, -D_INF); return GR_SUCCESS; }
        if (DF(_is_nan_inline)(x)) { DF(_set_d_inline)(res, D_NAN); return GR_SUCCESS; }
        return GR_DOMAIN;
    }
    DF(log)(res, x);
    return GR_SUCCESS;
}
static int GRF(atan)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(atan)(res, x); return GR_SUCCESS; }
static int GRF(im)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) { DF(_zero_inline)(res); return GR_SUCCESS; }

static void GRF(vec_init)(DF(ptr) vec, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        DF(_zero_inline)(vec + i);
}
static void GRF(vec_clear)(DF(ptr) vec, slong len, gr_ctx_t ctx) { }
static int GRF(vec_set)(DF(ptr) res, DF(srcptr) x, slong len, gr_ctx_t ctx)
{
    if (res != x)
        memcpy(res, x, len * sizeof(DF(struct)));
    return GR_SUCCESS;
}
static int GRF(vec_zero)(DF(ptr) vec, slong len, gr_ctx_t ctx) { GRF(vec_init)(vec, len, ctx); return GR_SUCCESS; }
static int GRF(vec_add)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, slong len, gr_ctx_t ctx) { DFV(vec_add)(res, x, y, len); return GR_SUCCESS; }
static int GRF(vec_sub)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, slong len, gr_ctx_t ctx) { DFV(vec_sub)(res, x, y, len); return GR_SUCCESS; }
static int GRF(vec_mul)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, slong len, gr_ctx_t ctx) { DFV(vec_mul)(res, x, y, len); return GR_SUCCESS; }
static int GRF(vec_mul_scalar)(DF(ptr) res, DF(srcptr) x, slong len, DF(srcptr) c, gr_ctx_t ctx) { DFV(vec_mul_scalar)(res, x, len, c); return GR_SUCCESS; }
static int GRF(vec_addmul_scalar)(DF(ptr) res, DF(srcptr) x, slong len, DF(srcptr) c, gr_ctx_t ctx) { DFV(vec_addmul_scalar)(res, x, len, c); return GR_SUCCESS; }
static int GRF(vec_submul_scalar)(DF(ptr) res, DF(srcptr) x, slong len, DF(srcptr) c, gr_ctx_t ctx) { DFV(vec_submul_scalar)(res, x, len, c); return GR_SUCCESS; }
static int GRF(vec_div)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        if (DF(_is_zero_inline)(y + i))
            return GR_DOMAIN;
    DFV(vec_div)(res, x, y, len);
    return GR_SUCCESS;
}
static int GRF(vec_dot)(DF(ptr) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, slong len, gr_ctx_t ctx) { return DFV(_vec_dot_dispatch)(res, initial, subtract, x, y, len, 0); }
static int GRF(vec_dot_rev)(DF(ptr) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, slong len, gr_ctx_t ctx) { return DFV(_vec_dot_dispatch)(res, initial, subtract, x, y, len, 1); }

/* elementwise elementary functions, conversions and data movement */
static int GRF(vec_exp)(DF(ptr) res, DF(srcptr) x, slong len, gr_ctx_t ctx) { DFV(vec_exp)(res, x, len); return GR_SUCCESS; }
static int GRF(vec_sin)(DF(ptr) res, DF(srcptr) x, slong len, gr_ctx_t ctx) { DFV(vec_sin)(res, x, len); return GR_SUCCESS; }
static int GRF(vec_cos)(DF(ptr) res, DF(srcptr) x, slong len, gr_ctx_t ctx) { DFV(vec_cos)(res, x, len); return GR_SUCCESS; }
static int GRF(vec_sin_cos)(DF(ptr) sn, DF(ptr) cs, DF(srcptr) x, slong len, gr_ctx_t ctx) { DFV(vec_sin_cos)(sn, cs, x, len); return GR_SUCCESS; }
static int GRF(vec_sqrt)(DF(ptr) res, DF(srcptr) x, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        if (x[i].d[0] < 0.0)
            return GR_DOMAIN;
    DFV(vec_sqrt)(res, x, len);
    return GR_SUCCESS;
}
static int GRF(vec_rsqrt)(DF(ptr) res, DF(srcptr) x, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        if (x[i].d[0] < 0.0 || DF(_is_zero_inline)(x + i))
            return GR_DOMAIN;
    DFV(vec_rsqrt)(res, x, len);
    return GR_SUCCESS;
}
static int GRF(vec_log)(DF(ptr) res, DF(srcptr) x, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!(x[i].d[0] > 0.0) && !DF(_is_zero_inline)(x + i) && !DF(_is_nan_inline)(x + i))
            return GR_DOMAIN;
    DFV(vec_log)(res, x, len);
    return GR_SUCCESS;
}
/* from another dfloat ring: exact when it has at most N components;
   otherwise the leading N components when the source is canonical
   (a strong context), and a renormalisation in general */
static int GRF(vec_set_other)(DF(ptr) res, gr_srcptr x, gr_ctx_t x_ctx, slong len, gr_ctx_t ctx)
{
    slong i;
    int k;
    if (x_ctx->which_ring == GR_CTX_DFLOAT || x_ctx->which_ring == GR_CTX_DFLOAT_BALL)
    {
        int nx = DFLOAT_CTX_N(x_ctx);
        slong sz = x_ctx->sizeof_elem;
        if (nx <= N || DFLOAT_CTX_STRONG(x_ctx))
        {
            int nc = FLINT_MIN(nx, N);
            for (i = 0; i < len; i++)
            {
                const double * xi = (const double *) GR_ENTRY(x, i, sz);
                for (k = 0; k < nc; k++)
                    res[i].d[k] = xi[k];
                for (; k < N; k++)
                    res[i].d[k] = 0.0;
            }
        }
        else
        {
            double t[DFLOAT_MAX_N];
            for (i = 0; i < len; i++)
            {
                const double * xi = (const double *) GR_ENTRY(x, i, sz);
                for (k = 0; k < nx; k++)
                    t[k] = xi[k];
                _dfloat_renorm(res[i].d, N, NULL, t, nx);
            }
        }
        return GR_SUCCESS;
    }
    return gr_generic_vec_set_other(res, x, x_ctx, len, ctx);
}
static int GRF(vec_gather)(DF(ptr) res, DF(srcptr) vec, const slong * idx, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        if (i + 16 < len)
            DFLOAT_PREFETCH(vec + idx[i + 16]);
        res[i] = vec[idx[i]];
    }
    return GR_SUCCESS;
}
static int GRF(vec_scatter)(DF(ptr) vec, const slong * idx, DF(srcptr) src, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        vec[idx[i]] = src[i];
    return GR_SUCCESS;
}

/* ---- strong (canonical) variants: the same, followed by a
   canonicalisation of the result ---- */

#define GRF_STRONG_BINARY(name) \
static int GRF(name##_strong)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, gr_ctx_t ctx) \
{ int status = GRF(name)(res, x, y, ctx); DF(_canonicalise_inplace)(res->d); return status; }
#define GRF_STRONG_UNARY(name) \
static int GRF(name##_strong)(DF(ptr) res, DF(srcptr) x, gr_ctx_t ctx) \
{ int status = GRF(name)(res, x, ctx); DF(_canonicalise_inplace)(res->d); return status; }
#define GRF_STRONG_VEC_BINARY(name) \
static int GRF(name##_strong)(DF(ptr) res, DF(srcptr) x, DF(srcptr) y, slong len, gr_ctx_t ctx) \
{ slong i; int status = GRF(name)(res, x, y, len, ctx); \
  for (i = 0; i < len; i++) { DF(_canonicalise_inplace)(res[i].d); } return status; }
#define GRF_STRONG_VEC_SCALAR(name) \
static int GRF(name##_strong)(DF(ptr) res, DF(srcptr) x, slong len, DF(srcptr) c, gr_ctx_t ctx) \
{ slong i; int status = GRF(name)(res, x, len, c, ctx); \
  for (i = 0; i < len; i++) { DF(_canonicalise_inplace)(res[i].d); } return status; }
#define GRF_STRONG_DOT(name) \
static int GRF(name##_strong)(DF(ptr) res, DF(srcptr) initial, int subtract, DF(srcptr) x, DF(srcptr) y, slong len, gr_ctx_t ctx) \
{ int status = GRF(name)(res, initial, subtract, x, y, len, ctx); DF(_canonicalise_inplace)(res->d); return status; }

GRF_STRONG_BINARY(add)
GRF_STRONG_BINARY(sub)
GRF_STRONG_BINARY(mul)
GRF_STRONG_UNARY(sqr)
GRF_STRONG_BINARY(div)
GRF_STRONG_UNARY(inv)
GRF_STRONG_UNARY(sqrt)
GRF_STRONG_UNARY(exp)
GRF_STRONG_UNARY(sin)
GRF_STRONG_UNARY(cos)
GRF_STRONG_UNARY(log)
GRF_STRONG_UNARY(atan)
static int GRF(sin_cos_strong)(DF(ptr) sn, DF(ptr) cs, DF(srcptr) x, gr_ctx_t ctx)
{ int status = GRF(sin_cos)(sn, cs, x, ctx); DF(_canonicalise_inplace)(sn->d); DF(_canonicalise_inplace)(cs->d); return status; }
#define GRF_STRONG_VEC_UNARY(name) \
static int GRF(name##_strong)(DF(ptr) res, DF(srcptr) x, slong len, gr_ctx_t ctx) \
{ slong i; int status = GRF(name)(res, x, len, ctx); \
  for (i = 0; i < len; i++) { DF(_canonicalise_inplace)(res[i].d); } return status; }
GRF_STRONG_UNARY(rsqrt)
GRF_STRONG_VEC_UNARY(vec_sqrt)
GRF_STRONG_VEC_UNARY(vec_rsqrt)
GRF_STRONG_VEC_UNARY(vec_exp)
GRF_STRONG_VEC_UNARY(vec_log)
GRF_STRONG_VEC_UNARY(vec_sin)
GRF_STRONG_VEC_UNARY(vec_cos)
static int GRF(vec_sin_cos_strong)(DF(ptr) sn, DF(ptr) cs, DF(srcptr) x, slong len, gr_ctx_t ctx)
{ slong i; int status = GRF(vec_sin_cos)(sn, cs, x, len, ctx);
  for (i = 0; i < len; i++) { DF(_canonicalise_inplace)(sn[i].d); DF(_canonicalise_inplace)(cs[i].d); } return status; }
GRF_STRONG_VEC_BINARY(vec_add)
GRF_STRONG_VEC_BINARY(vec_sub)
GRF_STRONG_VEC_BINARY(vec_mul)
GRF_STRONG_VEC_BINARY(vec_div)
GRF_STRONG_VEC_SCALAR(vec_mul_scalar)
GRF_STRONG_VEC_SCALAR(vec_addmul_scalar)
GRF_STRONG_VEC_SCALAR(vec_submul_scalar)
GRF_STRONG_DOT(vec_dot)
GRF_STRONG_DOT(vec_dot_rev)

/* overrides applied on top of the weak table */
gr_method_tab_input DFV(gr_strong_methods_input)[] =
{
    {GR_METHOD_ADD,             (gr_funcptr) GRF(add_strong)},
    {GR_METHOD_SUB,             (gr_funcptr) GRF(sub_strong)},
    {GR_METHOD_MUL,             (gr_funcptr) GRF(mul_strong)},
    {GR_METHOD_SQR,             (gr_funcptr) GRF(sqr_strong)},
    {GR_METHOD_DIV,             (gr_funcptr) GRF(div_strong)},
    {GR_METHOD_INV,             (gr_funcptr) GRF(inv_strong)},
    {GR_METHOD_SQRT,            (gr_funcptr) GRF(sqrt_strong)},
    {GR_METHOD_EXP,             (gr_funcptr) GRF(exp_strong)},
    {GR_METHOD_SIN,             (gr_funcptr) GRF(sin_strong)},
    {GR_METHOD_COS,             (gr_funcptr) GRF(cos_strong)},
    {GR_METHOD_SIN_COS,         (gr_funcptr) GRF(sin_cos_strong)},
    {GR_METHOD_LOG,             (gr_funcptr) GRF(log_strong)},
    {GR_METHOD_ATAN,            (gr_funcptr) GRF(atan_strong)},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) GRF(vec_add_strong)},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) GRF(vec_sub_strong)},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) GRF(vec_mul_strong)},
    {GR_METHOD_VEC_DIV,         (gr_funcptr) GRF(vec_div_strong)},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) GRF(vec_mul_scalar_strong)},
    {GR_METHOD_VEC_ADDMUL_SCALAR, (gr_funcptr) GRF(vec_addmul_scalar_strong)},
    {GR_METHOD_VEC_SUBMUL_SCALAR, (gr_funcptr) GRF(vec_submul_scalar_strong)},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) GRF(vec_dot_strong)},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) GRF(vec_dot_rev_strong)},
    {GR_METHOD_RSQRT,           (gr_funcptr) GRF(rsqrt_strong)},
    {GR_METHOD_VEC_SQRT,        (gr_funcptr) GRF(vec_sqrt_strong)},
    {GR_METHOD_VEC_RSQRT,       (gr_funcptr) GRF(vec_rsqrt_strong)},
    {GR_METHOD_VEC_EXP,         (gr_funcptr) GRF(vec_exp_strong)},
    {GR_METHOD_VEC_LOG,         (gr_funcptr) GRF(vec_log_strong)},
    {GR_METHOD_VEC_SIN,         (gr_funcptr) GRF(vec_sin_strong)},
    {GR_METHOD_VEC_COS,         (gr_funcptr) GRF(vec_cos_strong)},
    {GR_METHOD_VEC_SIN_COS,     (gr_funcptr) GRF(vec_sin_cos_strong)},
    {0,                         (gr_funcptr) NULL},
};

gr_method_tab_input DFV(gr_methods_input)[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) GRF(ctx_write)},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,     (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_APPROX_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,        (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_HAS_REAL_PREC,       (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC,       (gr_funcptr) GRF(ctx_set_real_prec)},
    {GR_METHOD_CTX_GET_REAL_PREC,       (gr_funcptr) GRF(ctx_get_real_prec)},
    {GR_METHOD_INIT,            (gr_funcptr) GRF(init)},
    {GR_METHOD_CLEAR,           (gr_funcptr) GRF(clear)},
    {GR_METHOD_SWAP,            (gr_funcptr) GRF(swap)},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) GRF(set)},
    {GR_METHOD_RANDTEST,        (gr_funcptr) GRF(randtest)},
    {GR_METHOD_WRITE,           (gr_funcptr) GRF(write)},
    {GR_METHOD_WRITE_N,         (gr_funcptr) GRF(write_n)},
    {GR_METHOD_ZERO,            (gr_funcptr) GRF(zero)},
    {GR_METHOD_ONE,             (gr_funcptr) GRF(one)},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) GRF(neg_one)},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) GRF(is_zero)},
    {GR_METHOD_IS_ONE,          (gr_funcptr) GRF(is_one)},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) GRF(is_neg_one)},
    {GR_METHOD_EQUAL,           (gr_funcptr) GRF(equal)},
    {GR_METHOD_SET,             (gr_funcptr) GRF(set)},
    {GR_METHOD_SET_SI,          (gr_funcptr) GRF(set_si)},
    {GR_METHOD_SET_UI,          (gr_funcptr) GRF(set_ui)},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) GRF(set_fmpz)},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) GRF(set_fmpq)},
    {GR_METHOD_SET_D,           (gr_funcptr) GRF(set_d)},
    {GR_METHOD_SET_STR,         (gr_funcptr) GRF(set_str)},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) GRF(set_other)},
    {GR_METHOD_GET_SI,          (gr_funcptr) GRF(get_si)},
    {GR_METHOD_GET_UI,          (gr_funcptr) GRF(get_ui)},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) GRF(get_fmpz)},
    {GR_METHOD_GET_D,           (gr_funcptr) GRF(get_d)},
    {GR_METHOD_NEG,             (gr_funcptr) GRF(neg)},
    {GR_METHOD_ADD,             (gr_funcptr) GRF(add)},
    {GR_METHOD_SUB,             (gr_funcptr) GRF(sub)},
    {GR_METHOD_MUL,             (gr_funcptr) GRF(mul)},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) GRF(mul_two)},
    {GR_METHOD_SQR,             (gr_funcptr) GRF(sqr)},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) GRF(mul_2exp_si)},
    {GR_METHOD_DIV,             (gr_funcptr) GRF(div)},
    {GR_METHOD_INV,             (gr_funcptr) GRF(inv)},
    {GR_METHOD_SQRT,            (gr_funcptr) GRF(sqrt)},
    {GR_METHOD_RSQRT,           (gr_funcptr) GRF(rsqrt)},
    {GR_METHOD_POS_INF,         (gr_funcptr) GRF(pos_inf)},
    {GR_METHOD_NEG_INF,         (gr_funcptr) GRF(neg_inf)},
    {GR_METHOD_UINF,            (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_UNDEFINED,       (gr_funcptr) GRF(nan)},
    {GR_METHOD_UNKNOWN,         (gr_funcptr) GRF(nan)},
    {GR_METHOD_ABS,             (gr_funcptr) GRF(abs)},
    {GR_METHOD_CONJ,            (gr_funcptr) GRF(set)},
    {GR_METHOD_RE,              (gr_funcptr) GRF(set)},
    {GR_METHOD_IM,              (gr_funcptr) GRF(im)},
    {GR_METHOD_SGN,             (gr_funcptr) GRF(sgn)},
    {GR_METHOD_CSGN,            (gr_funcptr) GRF(sgn)},
    {GR_METHOD_CMP,             (gr_funcptr) GRF(cmp)},
    {GR_METHOD_CMPABS,          (gr_funcptr) GRF(cmpabs)},
    {GR_METHOD_I,               (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_PI,              (gr_funcptr) GRF(pi)},
    {GR_METHOD_EULER,           (gr_funcptr) GRF(euler)},
    {GR_METHOD_CATALAN,         (gr_funcptr) GRF(catalan)},
    {GR_METHOD_EXP,             (gr_funcptr) GRF(exp)},
    {GR_METHOD_EXPM1,           (gr_funcptr) GRF(expm1)},
    {GR_METHOD_LOG,             (gr_funcptr) GRF(log)},
    {GR_METHOD_LOG1P,           (gr_funcptr) GRF(log1p)},
    {GR_METHOD_SIN,             (gr_funcptr) GRF(sin)},
    {GR_METHOD_COS,             (gr_funcptr) GRF(cos)},
    {GR_METHOD_SIN_COS,         (gr_funcptr) GRF(sin_cos)},
    {GR_METHOD_TAN,             (gr_funcptr) GRF(tan)},
    {GR_METHOD_SINH,            (gr_funcptr) GRF(sinh)},
    {GR_METHOD_COSH,            (gr_funcptr) GRF(cosh)},
    {GR_METHOD_TANH,            (gr_funcptr) GRF(tanh)},
    {GR_METHOD_ATAN,            (gr_funcptr) GRF(atan)},
    {GR_METHOD_ATAN2,           (gr_funcptr) GRF(atan2)},
    {GR_METHOD_POW,             (gr_funcptr) GRF(pow)},
    {GR_METHOD_ASIN,            (gr_funcptr) GRF(asin)},
    {GR_METHOD_ACOS,            (gr_funcptr) GRF(acos)},
    {GR_METHOD_VEC_INIT,        (gr_funcptr) GRF(vec_init)},
    {GR_METHOD_VEC_CLEAR,       (gr_funcptr) GRF(vec_clear)},
    {GR_METHOD_VEC_SET,         (gr_funcptr) GRF(vec_set)},
    {GR_METHOD_VEC_ZERO,        (gr_funcptr) GRF(vec_zero)},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) GRF(vec_add)},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) GRF(vec_sub)},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) GRF(vec_mul)},
    {GR_METHOD_VEC_DIV,         (gr_funcptr) GRF(vec_div)},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) GRF(vec_mul_scalar)},
    {GR_METHOD_VEC_ADDMUL_SCALAR, (gr_funcptr) GRF(vec_addmul_scalar)},
    {GR_METHOD_VEC_SUBMUL_SCALAR, (gr_funcptr) GRF(vec_submul_scalar)},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) GRF(vec_dot)},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) GRF(vec_dot_rev)},
    {GR_METHOD_VEC_SQRT,        (gr_funcptr) GRF(vec_sqrt)},
    {GR_METHOD_VEC_RSQRT,       (gr_funcptr) GRF(vec_rsqrt)},
    {GR_METHOD_VEC_EXP,         (gr_funcptr) GRF(vec_exp)},
    {GR_METHOD_VEC_LOG,         (gr_funcptr) GRF(vec_log)},
    {GR_METHOD_VEC_SIN,         (gr_funcptr) GRF(vec_sin)},
    {GR_METHOD_VEC_COS,         (gr_funcptr) GRF(vec_cos)},
    {GR_METHOD_VEC_SIN_COS,     (gr_funcptr) GRF(vec_sin_cos)},
    {GR_METHOD_VEC_SET_OTHER,   (gr_funcptr) GRF(vec_set_other)},
    {GR_METHOD_VEC_GATHER,      (gr_funcptr) GRF(vec_gather)},
    {GR_METHOD_VEC_SCATTER,     (gr_funcptr) GRF(vec_scatter)},
    {0,                         (gr_funcptr) NULL},
};

/* ---- ball ring ---- */

static int GRB(ctx_write)(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Real numbers as balls with ");
    status |= gr_stream_write_si(out, N);
    status |= gr_stream_write(out, "-term double expansion midpoints (d");
    status |= gr_stream_write_si(out, N);
    status |= gr_stream_write(out, DFLOAT_CTX_STRONG(ctx) ? "b, strong)" : "b)");
    return status;
}

static void GRB(init)(DB(ptr) x, gr_ctx_t ctx) { DB(_zero_inline)(x); }
static void GRB(clear)(DB(ptr) x, gr_ctx_t ctx) { }
static void GRB(swap)(DB(ptr) x, DB(ptr) y, gr_ctx_t ctx) { DB(struct) t = *x; *x = *y; *y = t; }

static int GRB(randtest)(DB(ptr) res, flint_rand_t state, gr_ctx_t ctx)
{
    if (n_randint(state, 16) == 0)
        DB(randtest_special)(res, state);
    else
        DB(randtest)(res, state);
    return GR_SUCCESS;
}

static int GRB(write)(gr_stream_t out, DB(srcptr) x, gr_ctx_t ctx)
{
    if (x->rad == 0.0)
    {
        if (DB(_is_zero_inline)(x)) return gr_stream_write(out, "0");
        if (x->d[0] == 1.0 && DF(_tail_abs_sum)(x->d) == 0.0) return gr_stream_write(out, "1");
        if (x->d[0] == -1.0 && DF(_tail_abs_sum)(x->d) == 0.0) return gr_stream_write(out, "-1");
    }
    return gr_stream_write_free(out, DB(get_str)(x, 0));
}

static int GRB(write_n)(gr_stream_t out, DB(srcptr) x, slong n, gr_ctx_t ctx)
{
    if (DB(_is_W)(x))
        return gr_stream_write_free(out, DB(get_str)(x, 0));
    return gr_stream_write_free(out, _dfloat_get_str(x->d, N, -1.0, FLINT_MAX(n, 1)));
}

static int GRB(zero)(DB(ptr) res, gr_ctx_t ctx) { DB(_zero_inline)(res); return GR_SUCCESS; }
static int GRB(one)(DB(ptr) res, gr_ctx_t ctx) { DB(_one_inline)(res); return GR_SUCCESS; }
static int GRB(neg_one)(DB(ptr) res, gr_ctx_t ctx) { DB(_set_d_inline)(res, -1.0); return GR_SUCCESS; }
static int GRB(set)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx)
{
    DB(_set_inline)(res, x);
    if (FLINT_UNLIKELY(DB(_is_W)(res)))
        DB(_indeterminate_inline)(res);
    return GR_SUCCESS;
}
static int GRB(set_si)(DB(ptr) res, slong x, gr_ctx_t ctx)
{
    fmpz_t t;
    if (DFLOAT_ULONG_IS_EXACT(FLINT_UABS(x)))
        DB(_set_d_inline)(res, (double) x);
    else
    {
        fmpz_init_set_si(t, x);
        DB(set_fmpz)(res, t);
        fmpz_clear(t);
    }
    return GR_SUCCESS;
}
static int GRB(set_ui)(DB(ptr) res, ulong x, gr_ctx_t ctx)
{
    fmpz_t t;
    if (DFLOAT_ULONG_IS_EXACT(x))
        DB(_set_d_inline)(res, (double) x);
    else
    {
        fmpz_init_set_ui(t, x);
        DB(set_fmpz)(res, t);
        fmpz_clear(t);
    }
    return GR_SUCCESS;
}
static int GRB(set_fmpz)(DB(ptr) res, const fmpz_t x, gr_ctx_t ctx) { DB(set_fmpz)(res, x); return GR_SUCCESS; }
static int GRB(set_fmpq)(DB(ptr) res, const fmpq_t x, gr_ctx_t ctx) { DB(set_fmpq)(res, x); return GR_SUCCESS; }
static int GRB(set_d)(DB(ptr) res, double x, gr_ctx_t ctx)
{
    if (!(fabs(x) <= DBL_MAX))
        return GR_DOMAIN;
    DB(_set_d_inline)(res, x);
    return GR_SUCCESS;
}
static int GRB(set_str)(DB(ptr) res, const char * s, gr_ctx_t ctx) { return DB(set_str)(res, s); }

static int GRB(set_other)(DB(ptr) res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            DB(set_fmpz)(res, x); return GR_SUCCESS;
        case GR_CTX_FMPQ:
            DB(set_fmpq)(res, x); return GR_SUCCESS;
        case GR_CTX_REAL_FLOAT_ARF:
            if (!arf_is_finite(x)) return GR_DOMAIN;
            DB(set_arf)(res, x); return GR_SUCCESS;
        case GR_CTX_RR_ARB:
            /* a nonfinite arb ball is the whole line */
            DB(set_arb)(res, x); return GR_SUCCESS;
        case GR_CTX_DFLOAT:
        {
            double err = 0.0;
            double t[DFLOAT_MAX_N];
            int i, nx = DFLOAT_CTX_N(x_ctx);
            if (!(_dfloat_abs_sum(x, nx) <= DBL_MAX)) return GR_DOMAIN;
            for (i = 0; i < nx; i++) t[i] = ((const double *) x)[i];
            for (; i < N; i++) t[i] = 0.0;
            _dfloat_renorm(res->d, N, &err, t, FLINT_MAX(nx, N));
            res->rad = _dfloat_rad_finish(err);
            return GR_SUCCESS;
        }
        case GR_CTX_DFLOAT_BALL:
            dfloat_ball_set_ball(res->d, N, (const double *) x, DFLOAT_CTX_N(x_ctx));
            return GR_SUCCESS;
        case GR_CTX_DFLOAT_COMPLEX:
        case GR_CTX_DFLOAT_COMPLEX_BALL:
        {
            /* the real part, if the imaginary part is exactly zero
               (GR_UNABLE if it only contains zero) */
            gr_ctx_t rctx;
            arb_t t;
            int nx = DFLOAT_CTX_N(x_ctx), ball = DFLOAT_CTX_BALL(x_ctx), status;
            arb_init(t);
            if (!_dfloat_part_get_arb(t, (const double *) x + nx + ball, nx, ball))
                status = GR_DOMAIN;
            else if (!arb_is_zero(t))
                status = arb_contains_zero(t) ? GR_UNABLE : GR_DOMAIN;
            else
            {
                GR_MUST_SUCCEED(gr_ctx_init_dfloat(rctx, nx, ball ? DFLOAT_BALL : 0));
                status = GRB(set_other)(res, x, rctx, ctx);
                gr_ctx_clear(rctx);
            }
            arb_clear(t);
            return status;
        }
        default:
        {
            gr_ctx_t actx;
            arb_t t;
            int status;
            gr_ctx_init_real_arb(actx, PREC_VIA_ARB);
            arb_init(t);
            status = gr_set_other(t, x, x_ctx, actx);
            if (status == GR_SUCCESS)
                DB(set_arb)(res, t);     /* nonfinite: the whole line */
            arb_clear(t);
            gr_ctx_clear(actx);
            return status;
        }
    }
}

static int GRB(set_interval_mid_rad)(DB(ptr) res, DB(srcptr) m, DB(srcptr) r, gr_ctx_t ctx)
{
    double rad = fabs(r->d[0]) + _dfloat_rad_finish(DF(_tail_abs_sum)(r->d) + r->rad);
    DB(_set_inline)(res, m);
    res->rad = _dfloat_rad_finish(res->rad + rad);
    DB(_fix_nonfinite)(res->d, &res->rad);
    return GR_SUCCESS;
}

static int GRB(get_fmpz)(fmpz_t res, DB(srcptr) x, gr_ctx_t ctx) { return DB(_get_fmpz)(res, x); }

static int GRB(get_si)(slong * res, DB(srcptr) x, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    status = GRB(get_fmpz)(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_fits_si(t))
            *res = fmpz_get_si(t);
        else
            status = GR_DOMAIN;
    }
    fmpz_clear(t);
    return status;
}

static int GRB(get_ui)(ulong * res, DB(srcptr) x, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;
    fmpz_init(t);
    status = GRB(get_fmpz)(t, x, ctx);
    if (status == GR_SUCCESS)
    {
        if (fmpz_sgn(t) >= 0 && fmpz_abs_fits_ui(t))
            *res = fmpz_get_ui(t);
        else
            status = GR_DOMAIN;
    }
    fmpz_clear(t);
    return status;
}

static int GRB(get_d)(double * res, DB(srcptr) x, gr_ctx_t ctx) { *res = DB(get_d)(x); return GR_SUCCESS; }

static truth_t GRB(is_zero)(DB(srcptr) x, gr_ctx_t ctx)
{
    if (DB(_is_zero_inline)(x)) return T_TRUE;
    if (DB(_contains_zero_inline)(x)) return T_UNKNOWN;
    return T_FALSE;
}
static truth_t GRB(_is_d)(DB(srcptr) x, double c)
{
    if (DB(_is_W)(x))
        return T_UNKNOWN;
    if (x->rad == 0.0)
    {
        if (x->d[0] == c && DF(_tail_abs_sum)(x->d) == 0.0) return T_TRUE;
        /* exact and not equal to c, unless a non-canonical
           representation hides it: decide through a subtraction */
        return DB(contains_d)(x, c) ? T_TRUE : T_FALSE;
    }
    {
        double lo, hi, d0;
        DF(_diff_bounds_d)(&lo, &hi, &d0, x->d, c);
        if (_dfloat_rad_finish(hi) <= x->rad) return T_UNKNOWN;     /* contains c */
        /* certainly not: |mid - c| > rad */
        if (lo > _dfloat_rad_finish(x->rad)) return T_FALSE;
        return T_UNKNOWN;
    }
}
static truth_t GRB(is_one)(DB(srcptr) x, gr_ctx_t ctx) { return GRB(_is_d)(x, 1.0); }
static truth_t GRB(is_neg_one)(DB(srcptr) x, gr_ctx_t ctx) { return GRB(_is_d)(x, -1.0); }
static truth_t GRB(equal)(DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx)
{
    if (DB(_equal_inline)(x, y)) return T_TRUE;
    if (DB(_overlaps_inline)(x, y)) return T_UNKNOWN;
    return T_FALSE;
}

static int GRB(neg)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(_neg_inline)(res, x); return GR_SUCCESS; }
static int GRB(add)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx) { DB(_add_impl)(res, x, y); return GR_SUCCESS; }
static int GRB(sub)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx) { DB(_sub_impl)(res, x, y); return GR_SUCCESS; }
static int GRB(mul)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx) { DB(_mul_impl)(res, x, y); return GR_SUCCESS; }
static int GRB(sqr)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(_sqr_impl)(res, x); return GR_SUCCESS; }
static int GRB(mul_two)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(_mul_2exp_si_impl)(res, x, 1); return GR_SUCCESS; }
static int GRB(mul_2exp_si)(DB(ptr) res, DB(srcptr) x, slong e, gr_ctx_t ctx) { DB(_mul_2exp_si_impl)(res, x, e); return GR_SUCCESS; }
static int GRB(div)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx) { return DB(_div_ool)(res, x, y, DFLOAT_MODE_AUTO); }
static int GRB(inv)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(_inv_m)(res, x, DFLOAT_MODE_AUTO); }
static truth_t GRB(is_invertible)(DB(srcptr) x, gr_ctx_t ctx)
{
    if (DB(_is_zero_inline)(x)) return T_FALSE;
    if (DB(_contains_zero_inline)(x)) return T_UNKNOWN;
    return T_TRUE;
}
static int GRB(sqrt)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(_sqrt_m)(res, x, DFLOAT_MODE_AUTO); }
static int GRB(rsqrt)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(_rsqrt_ool)(res, x, DFLOAT_MODE_AUTO); }
static int GRB(abs)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(_abs_inline)(res, x); return GR_SUCCESS; }
static int GRB(sgn)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx)
{
    /* the sign function is defined everywhere: [0 +/- 1] encloses its
       values on a ball that contains zero (or is the whole line) */
    if (DB(_is_zero_inline)(x)) { DB(_zero_inline)(res); return GR_SUCCESS; }
    if (DB(_contains_zero_inline)(x)) { DB(_set_zero_pm)(res, 1.0); return GR_SUCCESS; }
    DB(_set_d_inline)(res, (x->d[0] > 0.0) ? 1.0 : -1.0);
    return GR_SUCCESS;
}
static int GRB(cmp)(int * res, DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx) { return DB(cmp)(res, x, y); }
static int GRB(cmpabs)(int * res, DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx)
{
    DB(t) a, b;
    DB(_abs_inline)(a, x);
    DB(_abs_inline)(b, y);
    return DB(cmp)(res, a, b);
}
static int GRB(exp)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(exp)(res, x); return GR_SUCCESS; }
static int GRB(log)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(log)(res, x); }
static int GRB(atan)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(atan)(res, x); return GR_SUCCESS; }
static int GRB(im)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(_zero_inline)(res); return GR_SUCCESS; }

static void GRB(vec_init)(DB(ptr) vec, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        DB(_zero_inline)(vec + i);
}
static void GRB(vec_clear)(DB(ptr) vec, slong len, gr_ctx_t ctx) { }
static int GRB(vec_set)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx)
{
    if (res != x)
        memcpy(res, x, len * sizeof(DB(struct)));
    return GR_SUCCESS;
}
static int GRB(vec_zero)(DB(ptr) vec, slong len, gr_ctx_t ctx) { GRB(vec_init)(vec, len, ctx); return GR_SUCCESS; }
static int GRB(vec_add)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx) { DBV(vec_add)(res, x, y, len); return GR_SUCCESS; }
static int GRB(vec_sub)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx) { DBV(vec_sub)(res, x, y, len); return GR_SUCCESS; }
static int GRB(vec_mul)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx) { DBV(vec_mul)(res, x, y, len); return GR_SUCCESS; }
static int GRB(vec_mul_scalar)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c, gr_ctx_t ctx) { DBV(vec_mul_scalar)(res, x, len, c); return GR_SUCCESS; }
static int GRB(vec_addmul_scalar)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c, gr_ctx_t ctx) { DBV(vec_addmul_scalar)(res, x, len, c); return GR_SUCCESS; }
static int GRB(vec_submul_scalar)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c, gr_ctx_t ctx) { DBV(vec_submul_scalar)(res, x, len, c); return GR_SUCCESS; }
static int GRB(_vec_div_status)(DB(srcptr) y, slong len)
{
    slong i;
    int status = GR_SUCCESS;
    for (i = 0; i < len; i++)
    {
        if (FLINT_UNLIKELY(DB(_contains_zero_inline)(y + i)))
        {
            if (DB(_is_zero_inline)(y + i))
                return GR_DOMAIN;
            status = GR_UNABLE;
        }
    }
    return status;
}
static int GRB(vec_div)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx)
{
    int status = GRB(_vec_div_status)(y, len);
    if (status != GR_SUCCESS)
        return status;
    DBV(vec_div)(res, x, y, len);
    return GR_SUCCESS;
}
static int GRB(vec_dot)(DB(ptr) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx) { return DBV_vec_dot_fwd(res, initial, subtract, x, y, len, DFLOAT_MODE_AUTO); }
static int GRB(vec_dot_rev)(DB(ptr) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx) { return DBV_vec_dot_bwd(res, initial, subtract, x, y, len, DFLOAT_MODE_AUTO); }

/* elementwise elementary functions, conversions and data movement */
static int GRB(vec_exp)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx) { DBV(vec_exp)(res, x, len); return GR_SUCCESS; }
static int GRB(vec_log)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx) { return DBV(vec_log)(res, x, len); }
static int GRB(vec_sin)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx) { DBV(vec_sin)(res, x, len); return GR_SUCCESS; }
static int GRB(vec_cos)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx) { DBV(vec_cos)(res, x, len); return GR_SUCCESS; }
static int GRB(vec_sin_cos)(DB(ptr) sn, DB(ptr) cs, DB(srcptr) x, slong len, gr_ctx_t ctx) { DBV(vec_sin_cos)(sn, cs, x, len); return GR_SUCCESS; }
static int GRB(vec_sqrt)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx) { return DBV(vec_sqrt)(res, x, len); }
static int GRB(vec_rsqrt)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx) { return DBV(vec_rsqrt)(res, x, len); }
static int GRB(get_interval_mid_rad)(DB(ptr) m, DB(ptr) r, DB(srcptr) x, gr_ctx_t ctx)
{
    int i;
    if (DB(_is_W)(x))
    {
        /* the radius is not a real number */
        DB(_zero_inline)(m);
        DB(_indeterminate_inline)(r);
        return GR_UNABLE;
    }
    for (i = 0; i < N; i++)
        m->d[i] = x->d[i];
    m->rad = 0.0;
    r->d[0] = x->rad;
    for (i = 1; i < N; i++)
        r->d[i] = 0.0;
    r->rad = 0.0;
    return GR_SUCCESS;
}
static int GRB(vec_get_interval_mid_rad)(DB(ptr) m, DB(ptr) r, DB(srcptr) x, slong len, gr_ctx_t ctx)
{
    slong i;
    int status = GR_SUCCESS;
    for (i = 0; i < len; i++)
        status |= GRB(get_interval_mid_rad)(m + i, r + i, x + i, ctx);
    return status;
}
/* from another dfloat ring: exact when it has at most N components;
   otherwise the leading N components with the tail in the radius when
   the source is canonical (a strong context), and a renormalisation
   with tracked error in general */
static int GRB(vec_set_other)(DB(ptr) res, gr_srcptr x, gr_ctx_t x_ctx, slong len, gr_ctx_t ctx)
{
    slong i;
    int k;
    if (x_ctx->which_ring == GR_CTX_DFLOAT || x_ctx->which_ring == GR_CTX_DFLOAT_BALL)
    {
        int nx = DFLOAT_CTX_N(x_ctx);
        int ball = (x_ctx->which_ring == GR_CTX_DFLOAT_BALL);
        slong sz = x_ctx->sizeof_elem;
        if (nx <= N)
        {
            for (i = 0; i < len; i++)
            {
                const double * xi = (const double *) GR_ENTRY(x, i, sz);
                for (k = 0; k < nx; k++)
                    res[i].d[k] = xi[k];
                for (; k < N; k++)
                    res[i].d[k] = 0.0;
                res[i].rad = ball ? xi[nx] : 0.0;
            }
        }
        else if (DFLOAT_CTX_STRONG(x_ctx))
        {
            for (i = 0; i < len; i++)
            {
                const double * xi = (const double *) GR_ENTRY(x, i, sz);
                double tail = ball ? xi[nx] : 0.0;
                for (k = 0; k < N; k++)
                    res[i].d[k] = xi[k];
                for (; k < nx; k++)
                    tail += fabs(xi[k]);
                res[i].rad = _dfloat_rad_finish(tail);
                DB(_fix_nonfinite)(res[i].d, &res[i].rad);
            }
        }
        else
        {
            double t[DFLOAT_MAX_N];
            for (i = 0; i < len; i++)
            {
                const double * xi = (const double *) GR_ENTRY(x, i, sz);
                double err = 0.0, rad = ball ? xi[nx] : 0.0;
                for (k = 0; k < nx; k++)
                    t[k] = xi[k];
                _dfloat_renorm(res[i].d, N, &err, t, nx);
                res[i].rad = _dfloat_rad_finish(rad + err);
                DB(_fix_nonfinite)(res[i].d, &res[i].rad);
            }
        }
        return GR_SUCCESS;
    }
    return gr_generic_vec_set_other(res, x, x_ctx, len, ctx);
}
static int GRB(vec_gather)(DB(ptr) res, DB(srcptr) vec, const slong * idx, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        if (i + 16 < len)
            DFLOAT_PREFETCH(vec + idx[i + 16]);
        res[i] = vec[idx[i]];
    }
    return GR_SUCCESS;
}
static int GRB(vec_scatter)(DB(ptr) vec, const slong * idx, DB(srcptr) src, slong len, gr_ctx_t ctx)
{
    slong i;
    for (i = 0; i < len; i++)
        vec[idx[i]] = src[i];
    return GR_SUCCESS;
}

/* ---- fast variants (DFLOAT_FAST): products bound their rounding
   errors instead of tracking them ---- */

static int GRB(mul_fast)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx) { DB(_mul_mode)(res, x, y, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRB(sqr_fast)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { DB(_sqr_mode)(res, x, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRB(div_fast)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx) { return DB(_div_ool)(res, x, y, DFLOAT_MODE_FAST); }
static int GRB(inv_fast)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(_inv_m)(res, x, DFLOAT_MODE_FAST); }
static int GRB(sqrt_fast)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(_sqrt_m)(res, x, DFLOAT_MODE_FAST); }
static int GRB(rsqrt_fast)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) { return DB(_rsqrt_ool)(res, x, DFLOAT_MODE_FAST); }
static int GRB(vec_sqrt_fast)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx) { return DBV_vec_sqrt_mode(res, x, len, DFLOAT_MODE_FAST); }
static int GRB(vec_rsqrt_fast)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx) { return DBV(_vec_rsqrt_m)(res, x, len, DFLOAT_MODE_FAST); }
static int GRB(vec_mul_fast)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx) { DBV(_vec_mul_m)(res, x, y, 1, len, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRB(vec_mul_scalar_fast)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c, gr_ctx_t ctx) { DBV(_vec_mul_scalar_m)(res, x, len, c, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRB(vec_addmul_scalar_fast)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c, gr_ctx_t ctx) { DBV(_vec_addmul_scalar_m)(res, x, len, c, 0, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRB(vec_submul_scalar_fast)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c, gr_ctx_t ctx) { DBV(_vec_addmul_scalar_m)(res, x, len, c, 1, DFLOAT_MODE_FAST); return GR_SUCCESS; }
static int GRB(vec_div_fast)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx)
{
    int status = GRB(_vec_div_status)(y, len);
    if (status != GR_SUCCESS)
        return status;
    DBV(_vec_div_m)(res, x, y, len, DFLOAT_MODE_FAST);
    return GR_SUCCESS;
}
static int GRB(vec_dot_fast)(DB(ptr) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx) { return DBV_vec_dot_fwd(res, initial, subtract, x, y, len, DFLOAT_MODE_FAST); }
static int GRB(vec_dot_rev_fast)(DB(ptr) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx) { return DBV_vec_dot_bwd(res, initial, subtract, x, y, len, DFLOAT_MODE_FAST); }

gr_method_tab_input DBV(gr_fast_methods_input)[] =
{
    {GR_METHOD_MUL,             (gr_funcptr) GRB(mul_fast)},
    {GR_METHOD_SQR,             (gr_funcptr) GRB(sqr_fast)},
    {GR_METHOD_DIV,             (gr_funcptr) GRB(div_fast)},
    {GR_METHOD_INV,             (gr_funcptr) GRB(inv_fast)},
    {GR_METHOD_SQRT,            (gr_funcptr) GRB(sqrt_fast)},
    {GR_METHOD_RSQRT,           (gr_funcptr) GRB(rsqrt_fast)},
    {GR_METHOD_VEC_SQRT,        (gr_funcptr) GRB(vec_sqrt_fast)},
    {GR_METHOD_VEC_RSQRT,       (gr_funcptr) GRB(vec_rsqrt_fast)},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) GRB(vec_mul_fast)},
    {GR_METHOD_VEC_DIV,         (gr_funcptr) GRB(vec_div_fast)},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) GRB(vec_mul_scalar_fast)},
    {GR_METHOD_VEC_ADDMUL_SCALAR, (gr_funcptr) GRB(vec_addmul_scalar_fast)},
    {GR_METHOD_VEC_SUBMUL_SCALAR, (gr_funcptr) GRB(vec_submul_scalar_fast)},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) GRB(vec_dot_fast)},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) GRB(vec_dot_rev_fast)},
    {0,                         (gr_funcptr) NULL},
};

/* ---- strong (canonical) variants: the same, followed by a
   canonicalisation of the result ---- */

#define GRB_FASTSEL(name) (DFLOAT_CTX_FAST(ctx) ? GRB(name##_fast) : GRB(name))
#define GRB_STRONG_BINARY(name, sel) \
static int GRB(name##_strong)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, gr_ctx_t ctx) \
{ int status = sel(res, x, y, ctx); DF(_canonicalise_inplace)(res->d); return status; }
#define GRB_STRONG_UNARY(name, sel) \
static int GRB(name##_strong)(DB(ptr) res, DB(srcptr) x, gr_ctx_t ctx) \
{ int status = sel(res, x, ctx); DF(_canonicalise_inplace)(res->d); return status; }
#define GRB_STRONG_VEC_BINARY(name, sel) \
static int GRB(name##_strong)(DB(ptr) res, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx) \
{ slong i; int status = sel(res, x, y, len, ctx); \
  for (i = 0; i < len; i++) { DF(_canonicalise_inplace)(res[i].d); } return status; }
#define GRB_STRONG_VEC_SCALAR(name, sel) \
static int GRB(name##_strong)(DB(ptr) res, DB(srcptr) x, slong len, DB(srcptr) c, gr_ctx_t ctx) \
{ slong i; int status = sel(res, x, len, c, ctx); \
  for (i = 0; i < len; i++) { DF(_canonicalise_inplace)(res[i].d); } return status; }
#define GRB_STRONG_DOT(name, sel) \
static int GRB(name##_strong)(DB(ptr) res, DB(srcptr) initial, int subtract, DB(srcptr) x, DB(srcptr) y, slong len, gr_ctx_t ctx) \
{ int status = sel(res, initial, subtract, x, y, len, ctx); DF(_canonicalise_inplace)(res->d); return status; }

GRB_STRONG_BINARY(add, GRB(add))
GRB_STRONG_BINARY(sub, GRB(sub))
GRB_STRONG_BINARY(mul, GRB_FASTSEL(mul))
GRB_STRONG_UNARY(sqr, GRB_FASTSEL(sqr))
GRB_STRONG_BINARY(div, GRB_FASTSEL(div))
GRB_STRONG_UNARY(inv, GRB_FASTSEL(inv))
GRB_STRONG_UNARY(sqrt, GRB_FASTSEL(sqrt))
GRB_STRONG_UNARY(rsqrt, GRB_FASTSEL(rsqrt))
GRB_STRONG_UNARY(exp, GRB(exp))
GRB_STRONG_UNARY(sin, GRB(sin))
GRB_STRONG_UNARY(cos, GRB(cos))
GRB_STRONG_UNARY(log, GRB(log))
GRB_STRONG_UNARY(atan, GRB(atan))
static int GRB(sin_cos_strong)(DB(ptr) sn, DB(ptr) cs, DB(srcptr) x, gr_ctx_t ctx)
{ int status = GRB(sin_cos)(sn, cs, x, ctx); DF(_canonicalise_inplace)(sn->d); DF(_canonicalise_inplace)(cs->d); return status; }
#define GRB_STRONG_VEC_UNARY(name) \
static int GRB(name##_strong)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx) \
{ slong i; int status = GRB(name)(res, x, len, ctx); \
  for (i = 0; i < len; i++) { DF(_canonicalise_inplace)(res[i].d); } return status; }
#define GRB_STRONG_VEC_UNARY_SEL(name, sel) \
static int GRB(name##_strong)(DB(ptr) res, DB(srcptr) x, slong len, gr_ctx_t ctx) \
{ slong i; int status = sel(res, x, len, ctx); \
  for (i = 0; i < len; i++) { DF(_canonicalise_inplace)(res[i].d); } return status; }
GRB_STRONG_VEC_UNARY_SEL(vec_sqrt, GRB_FASTSEL(vec_sqrt))
GRB_STRONG_VEC_UNARY_SEL(vec_rsqrt, GRB_FASTSEL(vec_rsqrt))
GRB_STRONG_VEC_UNARY(vec_exp)
GRB_STRONG_VEC_UNARY(vec_log)
GRB_STRONG_VEC_UNARY(vec_sin)
GRB_STRONG_VEC_UNARY(vec_cos)
static int GRB(vec_sin_cos_strong)(DB(ptr) sn, DB(ptr) cs, DB(srcptr) x, slong len, gr_ctx_t ctx)
{ slong i; int status = GRB(vec_sin_cos)(sn, cs, x, len, ctx);
  for (i = 0; i < len; i++) { DF(_canonicalise_inplace)(sn[i].d); DF(_canonicalise_inplace)(cs[i].d); } return status; }
GRB_STRONG_VEC_BINARY(vec_add, GRB(vec_add))
GRB_STRONG_VEC_BINARY(vec_sub, GRB(vec_sub))
GRB_STRONG_VEC_BINARY(vec_mul, GRB_FASTSEL(vec_mul))
GRB_STRONG_VEC_BINARY(vec_div, GRB_FASTSEL(vec_div))
GRB_STRONG_VEC_SCALAR(vec_mul_scalar, GRB_FASTSEL(vec_mul_scalar))
GRB_STRONG_VEC_SCALAR(vec_addmul_scalar, GRB_FASTSEL(vec_addmul_scalar))
GRB_STRONG_VEC_SCALAR(vec_submul_scalar, GRB_FASTSEL(vec_submul_scalar))
GRB_STRONG_DOT(vec_dot, GRB_FASTSEL(vec_dot))
GRB_STRONG_DOT(vec_dot_rev, GRB_FASTSEL(vec_dot_rev))

/* overrides applied on top of the weak table */
gr_method_tab_input DBV(gr_strong_methods_input)[] =
{
    {GR_METHOD_ADD,             (gr_funcptr) GRB(add_strong)},
    {GR_METHOD_SUB,             (gr_funcptr) GRB(sub_strong)},
    {GR_METHOD_MUL,             (gr_funcptr) GRB(mul_strong)},
    {GR_METHOD_SQR,             (gr_funcptr) GRB(sqr_strong)},
    {GR_METHOD_DIV,             (gr_funcptr) GRB(div_strong)},
    {GR_METHOD_INV,             (gr_funcptr) GRB(inv_strong)},
    {GR_METHOD_SQRT,            (gr_funcptr) GRB(sqrt_strong)},
    {GR_METHOD_RSQRT,           (gr_funcptr) GRB(rsqrt_strong)},
    {GR_METHOD_VEC_SQRT,        (gr_funcptr) GRB(vec_sqrt_strong)},
    {GR_METHOD_VEC_RSQRT,       (gr_funcptr) GRB(vec_rsqrt_strong)},
    {GR_METHOD_EXP,             (gr_funcptr) GRB(exp_strong)},
    {GR_METHOD_SIN,             (gr_funcptr) GRB(sin_strong)},
    {GR_METHOD_COS,             (gr_funcptr) GRB(cos_strong)},
    {GR_METHOD_SIN_COS,         (gr_funcptr) GRB(sin_cos_strong)},
    {GR_METHOD_LOG,             (gr_funcptr) GRB(log_strong)},
    {GR_METHOD_ATAN,            (gr_funcptr) GRB(atan_strong)},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) GRB(vec_add_strong)},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) GRB(vec_sub_strong)},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) GRB(vec_mul_strong)},
    {GR_METHOD_VEC_DIV,         (gr_funcptr) GRB(vec_div_strong)},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) GRB(vec_mul_scalar_strong)},
    {GR_METHOD_VEC_ADDMUL_SCALAR, (gr_funcptr) GRB(vec_addmul_scalar_strong)},
    {GR_METHOD_VEC_SUBMUL_SCALAR, (gr_funcptr) GRB(vec_submul_scalar_strong)},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) GRB(vec_dot_strong)},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) GRB(vec_dot_rev_strong)},
    {GR_METHOD_VEC_EXP,         (gr_funcptr) GRB(vec_exp_strong)},
    {GR_METHOD_VEC_LOG,         (gr_funcptr) GRB(vec_log_strong)},
    {GR_METHOD_VEC_SIN,         (gr_funcptr) GRB(vec_sin_strong)},
    {GR_METHOD_VEC_COS,         (gr_funcptr) GRB(vec_cos_strong)},
    {GR_METHOD_VEC_SIN_COS,     (gr_funcptr) GRB(vec_sin_cos_strong)},
    {0,                         (gr_funcptr) NULL},
};

gr_method_tab_input DBV(gr_methods_input)[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) GRB(ctx_write)},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,        (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_HAS_REAL_PREC,       (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_REAL_PREC,       (gr_funcptr) GRF(ctx_set_real_prec)},
    {GR_METHOD_CTX_GET_REAL_PREC,       (gr_funcptr) GRF(ctx_get_real_prec)},
    {GR_METHOD_INIT,            (gr_funcptr) GRB(init)},
    {GR_METHOD_CLEAR,           (gr_funcptr) GRB(clear)},
    {GR_METHOD_SWAP,            (gr_funcptr) GRB(swap)},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) GRB(set)},
    {GR_METHOD_RANDTEST,        (gr_funcptr) GRB(randtest)},
    {GR_METHOD_WRITE,           (gr_funcptr) GRB(write)},
    {GR_METHOD_WRITE_N,         (gr_funcptr) GRB(write_n)},
    {GR_METHOD_ZERO,            (gr_funcptr) GRB(zero)},
    {GR_METHOD_ONE,             (gr_funcptr) GRB(one)},
    {GR_METHOD_NEG_ONE,         (gr_funcptr) GRB(neg_one)},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) GRB(is_zero)},
    {GR_METHOD_IS_ONE,          (gr_funcptr) GRB(is_one)},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) GRB(is_neg_one)},
    {GR_METHOD_EQUAL,           (gr_funcptr) GRB(equal)},
    {GR_METHOD_SET,             (gr_funcptr) GRB(set)},
    {GR_METHOD_SET_SI,          (gr_funcptr) GRB(set_si)},
    {GR_METHOD_SET_UI,          (gr_funcptr) GRB(set_ui)},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) GRB(set_fmpz)},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) GRB(set_fmpq)},
    {GR_METHOD_SET_D,           (gr_funcptr) GRB(set_d)},
    {GR_METHOD_SET_STR,         (gr_funcptr) GRB(set_str)},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) GRB(set_other)},
    {GR_METHOD_SET_INTERVAL_MID_RAD, (gr_funcptr) GRB(set_interval_mid_rad)},
    {GR_METHOD_GET_SI,          (gr_funcptr) GRB(get_si)},
    {GR_METHOD_GET_UI,          (gr_funcptr) GRB(get_ui)},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) GRB(get_fmpz)},
    {GR_METHOD_GET_D,           (gr_funcptr) GRB(get_d)},
    {GR_METHOD_NEG,             (gr_funcptr) GRB(neg)},
    {GR_METHOD_ADD,             (gr_funcptr) GRB(add)},
    {GR_METHOD_SUB,             (gr_funcptr) GRB(sub)},
    {GR_METHOD_MUL,             (gr_funcptr) GRB(mul)},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) GRB(mul_two)},
    {GR_METHOD_SQR,             (gr_funcptr) GRB(sqr)},
    {GR_METHOD_MUL_2EXP_SI,     (gr_funcptr) GRB(mul_2exp_si)},
    {GR_METHOD_DIV,             (gr_funcptr) GRB(div)},
    {GR_METHOD_INV,             (gr_funcptr) GRB(inv)},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) GRB(is_invertible)},
    {GR_METHOD_SQRT,            (gr_funcptr) GRB(sqrt)},
    {GR_METHOD_RSQRT,           (gr_funcptr) GRB(rsqrt)},
    {GR_METHOD_ABS,             (gr_funcptr) GRB(abs)},
    {GR_METHOD_CONJ,            (gr_funcptr) GRB(set)},
    {GR_METHOD_RE,              (gr_funcptr) GRB(set)},
    {GR_METHOD_IM,              (gr_funcptr) GRB(im)},
    {GR_METHOD_SGN,             (gr_funcptr) GRB(sgn)},
    {GR_METHOD_CSGN,            (gr_funcptr) GRB(sgn)},
    {GR_METHOD_CMP,             (gr_funcptr) GRB(cmp)},
    {GR_METHOD_CMPABS,          (gr_funcptr) GRB(cmpabs)},
    {GR_METHOD_I,               (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_PI,              (gr_funcptr) GRB(pi)},
    {GR_METHOD_EULER,           (gr_funcptr) GRB(euler)},
    {GR_METHOD_CATALAN,         (gr_funcptr) GRB(catalan)},
    {GR_METHOD_EXP,             (gr_funcptr) GRB(exp)},
    {GR_METHOD_EXPM1,           (gr_funcptr) GRB(expm1)},
    {GR_METHOD_LOG,             (gr_funcptr) GRB(log)},
    {GR_METHOD_LOG1P,           (gr_funcptr) GRB(log1p)},
    {GR_METHOD_SIN,             (gr_funcptr) GRB(sin)},
    {GR_METHOD_COS,             (gr_funcptr) GRB(cos)},
    {GR_METHOD_SIN_COS,         (gr_funcptr) GRB(sin_cos)},
    {GR_METHOD_TAN,             (gr_funcptr) GRB(tan)},
    {GR_METHOD_SINH,            (gr_funcptr) GRB(sinh)},
    {GR_METHOD_COSH,            (gr_funcptr) GRB(cosh)},
    {GR_METHOD_TANH,            (gr_funcptr) GRB(tanh)},
    {GR_METHOD_ATAN,            (gr_funcptr) GRB(atan)},
    {GR_METHOD_ATAN2,           (gr_funcptr) GRB(atan2)},
    {GR_METHOD_POW,             (gr_funcptr) GRB(pow)},
    {GR_METHOD_ASIN,            (gr_funcptr) GRB(asin)},
    {GR_METHOD_ACOS,            (gr_funcptr) GRB(acos)},
    {GR_METHOD_VEC_INIT,        (gr_funcptr) GRB(vec_init)},
    {GR_METHOD_VEC_CLEAR,       (gr_funcptr) GRB(vec_clear)},
    {GR_METHOD_VEC_SET,         (gr_funcptr) GRB(vec_set)},
    {GR_METHOD_VEC_ZERO,        (gr_funcptr) GRB(vec_zero)},
    {GR_METHOD_VEC_ADD,         (gr_funcptr) GRB(vec_add)},
    {GR_METHOD_VEC_SUB,         (gr_funcptr) GRB(vec_sub)},
    {GR_METHOD_VEC_MUL,         (gr_funcptr) GRB(vec_mul)},
    {GR_METHOD_VEC_DIV,         (gr_funcptr) GRB(vec_div)},
    {GR_METHOD_VEC_MUL_SCALAR,  (gr_funcptr) GRB(vec_mul_scalar)},
    {GR_METHOD_VEC_ADDMUL_SCALAR, (gr_funcptr) GRB(vec_addmul_scalar)},
    {GR_METHOD_VEC_SUBMUL_SCALAR, (gr_funcptr) GRB(vec_submul_scalar)},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) GRB(vec_dot)},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) GRB(vec_dot_rev)},
    {GR_METHOD_VEC_SQRT,        (gr_funcptr) GRB(vec_sqrt)},
    {GR_METHOD_VEC_RSQRT,       (gr_funcptr) GRB(vec_rsqrt)},
    {GR_METHOD_VEC_EXP,         (gr_funcptr) GRB(vec_exp)},
    {GR_METHOD_VEC_LOG,         (gr_funcptr) GRB(vec_log)},
    {GR_METHOD_VEC_SIN,         (gr_funcptr) GRB(vec_sin)},
    {GR_METHOD_VEC_COS,         (gr_funcptr) GRB(vec_cos)},
    {GR_METHOD_VEC_SIN_COS,     (gr_funcptr) GRB(vec_sin_cos)},
    {GR_METHOD_VEC_SET_OTHER,   (gr_funcptr) GRB(vec_set_other)},
    {GR_METHOD_VEC_GATHER,      (gr_funcptr) GRB(vec_gather)},
    {GR_METHOD_VEC_SCATTER,     (gr_funcptr) GRB(vec_scatter)},
    {GR_METHOD_GET_INTERVAL_MID_RAD, (gr_funcptr) GRB(get_interval_mid_rad)},
    {GR_METHOD_VEC_GET_INTERVAL_MID_RAD, (gr_funcptr) GRB(vec_get_interval_mid_rad)},
    {0,                         (gr_funcptr) NULL},
};
