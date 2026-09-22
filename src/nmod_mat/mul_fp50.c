/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Matrix multiplication over Z/nZ for n < 2^50 in double precision, with
    the modular multiplication of fft_small.

    Every product is reduced as it is formed: the entries are lifted to
    doubles in the symmetric range |a| <= n/2, and one "k step" computes

        acc += mulmod(a, b) = a*b - n*rint(a*b/n)

    with the formula of machine_vectors.h (a multiply for the rounded
    product h, one for the quotient, its rounding, and two fused
    multiply-adds for h - q*n and h - a*b whose difference is the result,
    plus the accumulation, see mul_fp_vec.h); the mulmod result lies in
    (-9/8 n, 9/8 n), see fft_small/mulmod_satisfies_bounds.c. Up to 2^53/n
    such terms can then be summed exactly, so the accumulators are folded
    back to (-0.51 n, 0.51 n) every "cadence" steps with a reduce_pm1n
    (three operations), which for n near 2^50 is every 6 steps and below
    2^45 never inside a block. This is the arithmetic of fft_small applied to
    a matrix product: it covers moduli up to 50 bits in one pass, where
    nmod_mat_mul_blas needs several dgemm passes and a CRT, on any machine
    with a vector FMA. Its cost is 6-7 floating point operations per
    product against 1 for a dgemm, 2 for nmod_mat_mul_u32 and 3 for the
    integer two-limb kernel of mul_k52.c. The latter is the alternative on the
    same machines and how they compare depends on the architecture;
    use the profile p-mul_tune.c to compare them on complete multiplications.

    The packing, the register-tile microkernel, the blocked core and the
    thread split are those of mul_blocked_templ.h; the vector primitives
    are those of mul_fp_vec.h (AVX-512, AVX2, NEON, plain C).
*/

#include <string.h>
#include "nmod_mat.h"

#if FLINT_BITS == 64

#include "nmod.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "nmod_mat/mul_fp_vec.h"

/* tile geometry and blocking ************************************************/

/*
    MR*NACC accumulators, NACC B operands, the A operands of the rows,
    the constants n, 1/n (and the rounding constant on AVX2) and the
    temporaries of the mulmods in flight must fit the register file; the
    geometries below were checked on the generated code (no spill in the
    k loop on AVX-512 and AVX2, two on NEON), they are not tuned.
*/
#define FP50_VL FPV_VL
#if defined(FPV_AVX512)
# ifndef FP50_MR
#  define FP50_MR 8
# endif
# ifndef FP50_NACC
#  define FP50_NACC 2
# endif
#elif defined(FPV_AVX2)
# ifndef FP50_MR
#  define FP50_MR 3
# endif
# ifndef FP50_NACC
#  define FP50_NACC 2
# endif
#elif defined(FPV_NEON)
# ifndef FP50_MR
#  define FP50_MR 4
# endif
# ifndef FP50_NACC
#  define FP50_NACC 2
# endif
#else
# ifndef FP50_MR
#  define FP50_MR 4
# endif
# ifndef FP50_NACC
#  define FP50_NACC 2
# endif
#endif

#ifndef FP50_KC
# define FP50_KC 256
#endif
#ifndef FP50_MC
# define FP50_MC 96
#endif
#ifndef FP50_NC
# define FP50_NC 2048
#endif
/* products per thread, as U32_MT_MIN_WORK in mul_u32.c */
#ifndef FP50_MT_MIN_WORK
# define FP50_MT_MIN_WORK 500000.0
#endif

/* modulus-derived parameters ************************************************/

typedef struct
{
    ulong n;
    double nd;
    double ninv;
    slong cadence;
}
fp50_ctx_struct;

static void
fp50_ctx_init(fp50_ctx_struct * ctx, ulong n)
{
    double cad;

    FLINT_ASSERT(n >= 1 && n < (UWORD(1) << 50));

    ctx->n = n;
    ctx->nd = (double) n;
    ctx->ninv = 1.0 / (double) n;

    /*
        The accumulator starts from a canonical entry of C or from a fold,
        so |acc| <= n, and each step adds less than 9/8 n in absolute value:
        n + cadence * 9/8 n < 2^53 keeps every partial sum exact.
    */
    cad = (0x1.0p53 / (double) n - 1.0) * (8.0 / 9.0);
    if (cad >= (double) FP50_KC)
        ctx->cadence = FP50_KC;
    else
        ctx->cadence = FLINT_MAX((slong) cad, 1);
}

/* symmetric lift of a canonical residue: a - n if a > n/2, as a double */
FLINT_FORCE_INLINE double
fp50_lift(ulong a, ulong n)
{
    return (a > n / 2) ? (double) ((slong) a - (slong) n) : (double) a;
}

/* primitives ****************************************************************/

typedef struct
{
    fpv nv;
    fpv ninvv;
}
fp50_consts;

FLINT_FORCE_INLINE fp50_consts
fp50_consts_init(const fp50_ctx_struct * ctx)
{
    fp50_consts C;

    C.nv = fpv_set1(ctx->nd);
    C.ninvv = fpv_set1(ctx->ninv);

    return C;
}

FLINT_FORCE_INLINE void
fp50_load_bstep(fpv * bv, const double * p)
{
    slong v;

    for (v = 0; v < FP50_NACC; v++)
        bv[v] = fpv_load(p + v * FP50_VL);
}

FLINT_FORCE_INLINE fpv
fp50_mul_add(fpv acc, fpv a, fpv b, const fp50_consts * C)
{
    return fpv_add(acc, fpv_mulmod(a, b, C->nv, C->ninvv));
}

FLINT_FORCE_INLINE fpv
fp50_fold(fpv acc, const fp50_consts * C)
{
    return fpv_reduce_pm1n(acc, C->nv, C->ninvv);
}

FLINT_FORCE_INLINE fpv
fp50_finish(fpv acc, const fp50_consts * C)
{
    return fpv_reduce_0n(fpv_reduce_pm1n(acc, C->nv, C->ninvv), C->nv);
}

/* template instantiation ****************************************************/

#define BT_NAME(x) fp50_##x
#define BT_ENTRY ulong
#define BT_PACKED double
#define BT_CTX fp50_ctx_struct
#define BT_MR FP50_MR
#define BT_VL FP50_VL
#define BT_NACC FP50_NACC
#define BT_KC FP50_KC
#define BT_MC FP50_MC
#define BT_NC FP50_NC
#define BT_MT_MIN_WORK FP50_MT_MIN_WORK
#define BT_LIFT(x, ctx) fp50_lift((x), (ctx)->n)
#define BT_BSLOT(j) (j)
#define BT_ACC fpv
#define BT_BV fpv
#define BT_AV fpv
#define BT_CONSTS fp50_consts
#define BT_CONSTS_INIT(ctx) fp50_consts_init(ctx)
#define BT_ACC_ZERO() fpv_zero()
#define BT_LOAD_C(p) fpv_load_u64(p)
#define BT_STORE_C(p, acc) fpv_store_u64(p, acc)
#define BT_LOAD_BSTEP(bv, p) fp50_load_bstep(bv, p)
#define BT_LOAD_A(p) fpv_set1(*(p))
#define BT_MUL_ADD(acc, a, b, C) fp50_mul_add(acc, a, b, C)
#define BT_CADENCE(ctx) ((ctx)->cadence)
#define BT_FOLD(acc, C, ctx) fp50_fold(acc, C)
#define BT_FINISH(acc, C, ctx) fp50_finish(acc, C)
#include "mul_blocked_templ.h"

/* public entry *************************************************************/

int
nmod_mat_mul_fp50(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong m = A->r;
    slong k = A->c;
    slong n = B->c;
    ulong modn = C->mod.n;
    fp50_ctx_struct ctx;

    FLINT_ASSERT(C->r == A->r);
    FLINT_ASSERT(C->c == B->c);
    FLINT_ASSERT(A->c == B->r);

    if (modn >= (UWORD(1) << 50))
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
        nmod_mat_mul_fp50(T, A, B);
        nmod_mat_swap_entrywise(C, T);
        nmod_mat_clear(T);
        return 1;
    }

    fp50_ctx_init(&ctx, modn);

    fp50_core_mt(C->entries, C->stride, A->entries, A->stride,
                 B->entries, B->stride, m, k, n, &ctx,
                 flint_get_num_threads());

    return 1;
}

#else

/* the conversions here assume 64-bit words */
int
nmod_mat_mul_fp50(nmod_mat_t FLINT_UNUSED(C),
                  const nmod_mat_t FLINT_UNUSED(A),
                  const nmod_mat_t FLINT_UNUSED(B))
{
    return 0;
}

#endif
