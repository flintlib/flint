/*
    Copyright (C) 2026 Brian Heckel

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Speedup sweep: kernel-only 4x4 matmul (conversion excluded), batched
    Montgomery (lazy) vs gr_mat_mul (classical), across modulus bit sizes.

    The batched kernel is compile-time specialized on the 31-bit limb count
    s = ceil((bits+2)/31) (bm_matmul dispatches on ctx->s), so its cost steps up
    with s; classical steps up with the modulus's 64-bit limb count (MPN_MOD
    nlimbs).  The two step at different boundaries, so the speedup is largest
    just after an s-step and smallest just before the next one.  n < 2^246 (lazy
    single-subtract validity).
*/

#include <stdlib.h>
#include "fmpz.h"
#include "gr.h"
#include "gr_mat.h"
#include "mpn_mod.h"
#include "profiler.h"
#include "mpn_mod_batched_mont.h"

#define DIM  4
#define REPS 100000

static volatile ulong sink;

static void
randmat(gr_mat_t mat, const fmpz_t N, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_t t;
    slong i, j;
    fmpz_init(t);
    for (i = 0; i < DIM; i++)
        for (j = 0; j < DIM; j++)
        {
            fmpz_randm(t, state, N);
            GR_MUST_SUCCEED(gr_set_fmpz(gr_mat_entry_ptr(mat, i, j, ctx), t, ctx));
        }
    fmpz_clear(t);
}

FLINT_STATIC_NOINLINE ulong
run_classical(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, slong reps, gr_ctx_t ctx)
{
    slong k;
    for (k = 0; k < reps; k++)
        GR_MUST_SUCCEED(gr_mat_mul_classical(C, A, B, ctx));
    return ((nn_srcptr) gr_mat_entry_srcptr(C, 0, 0, ctx))[0];
}

#ifdef __AVX2__
FLINT_STATIC_NOINLINE ulong
run_kernel(bm_mat Cm, const bm_mat Am, const bm_mat Bm, slong reps, const bm_ctx_t * bc)
{
    slong k;
    for (k = 0; k < reps; k++)
        bm_matmul(Cm, Am, Bm, bc);
    return Cm[0][0][0];
}
#endif

static void
bench(slong bits, flint_rand_t state)
{
    gr_ctx_t ctx;
    fmpz_t N;
    gr_mat_t A, B, C;
    slong nlimbs;
    double tcpu, t_classical;
#ifdef __AVX2__
    double t_kernel;
    bm_ctx_t bc;
    bm_mat mA, mB, mC;
#endif

    fmpz_init(N);
    fmpz_randtest_unsigned(N, state, bits);
    fmpz_setbit(N, bits - 1);
    fmpz_setbit(N, 0);
    GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, N));
    nlimbs = MPN_MOD_CTX_NLIMBS(ctx);

    gr_mat_init(A, DIM, DIM, ctx);
    gr_mat_init(B, DIM, DIM, ctx);
    gr_mat_init(C, DIM, DIM, ctx);
    randmat(A, N, state, ctx);
    randmat(B, N, state, ctx);

    TIMEIT_START;
    sink ^= run_classical(C, A, B, REPS, ctx);
    TIMEIT_STOP_VALUES(tcpu, t_classical);

#ifdef __AVX2__
    bm_ctx_init(&bc, ctx);
    bm_load(mA, A, ctx, &bc);       /* converted once, outside the timed region */
    bm_load(mB, B, ctx, &bc);
    TIMEIT_START;
    sink ^= run_kernel(mC, (const ulong (*)[DIM][BM_MAXLIMB]) mA,
                           (const ulong (*)[DIM][BM_MAXLIMB]) mB, REPS, &bc);
    TIMEIT_STOP_VALUES(tcpu, t_kernel);
#endif

    flint_printf("%4wd bits  (nlimbs=%wd):  classical %7.1f ns", bits, nlimbs,
                 t_classical / REPS / 1e-9);
#ifdef __AVX2__
    flint_printf(",  batched %7.1f ns,  speedup %.2fx",
                 t_kernel / REPS / 1e-9, t_classical / t_kernel);
#else
    flint_printf(",  (no AVX2)");
#endif
    flint_printf("\n");

    (void) tcpu;
    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);
    gr_ctx_clear(ctx);
    fmpz_clear(N);
}

int main(void)
{
    flint_rand_t state;
    slong sizes[] = { 65, 80, 96, 112, 127, 129, 144, 160, 176, 191,
                      193, 208, 224, 240, 246 };
    int i, nsizes = sizeof(sizes) / sizeof(sizes[0]);

    flint_rand_init(state);
    flint_printf("kernel-only 4x4 matmul speedup vs modulus size "
                 "(per-s specialized lazy kernel, conversion excluded):\n\n");
    for (i = 0; i < nsizes; i++)
        bench(sizes[i], state);
    flint_rand_clear(state);
    return 0;
}
