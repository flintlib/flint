/*
    Copyright (C) 2026 Brian Heckel

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Chained multiplication: 4x4 matrix power A^exp over odd moduli < 2^246.
    Batched Montgomery (AVX2) converts to/from Montgomery form once and runs the
    whole square-and-multiply chain in that form, vs gr_mat_pow_ui (which for 4x4
    chains gr_mat_mul_classical on plain residues).
*/

#include <stdlib.h>
#include "fmpz.h"
#include "gr.h"
#include "gr_mat.h"
#include "mpn_mod.h"
#include "profiler.h"

#define DIM  4
#define REPS 500

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
run_gr(gr_mat_t C, const gr_mat_t A, ulong exp, slong reps, gr_ctx_t ctx)
{
    slong k;
    for (k = 0; k < reps; k++)
        GR_MUST_SUCCEED(gr_mat_pow_ui(C, A, exp, ctx));
    return ((nn_srcptr) gr_mat_entry_srcptr(C, 0, 0, ctx))[0];
}

#ifdef __AVX2__
FLINT_STATIC_NOINLINE ulong
run_batched(gr_mat_t C, const gr_mat_t A, ulong exp, slong reps, gr_ctx_t ctx)
{
    slong k;
    for (k = 0; k < reps; k++)
        GR_MUST_SUCCEED(mpn_mod_mat_pow_ui_batched_mont(C, A, exp, ctx));
    return ((nn_srcptr) gr_mat_entry_srcptr(C, 0, 0, ctx))[0];
}
#endif

static void
bench(slong bits, ulong exp, flint_rand_t state)
{
    gr_ctx_t ctx;
    fmpz_t N;
    gr_mat_t A, C;
    // TODO: Don't know how to get rid of the warning that this is unused
    double tcpu, t_gr;
#ifdef __AVX2__
    double t_batched;
#endif

    fmpz_init(N);
    fmpz_randtest_unsigned(N, state, bits);
    fmpz_setbit(N, bits - 1);
    fmpz_setbit(N, 0);
    GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, N));

    gr_mat_init(A, DIM, DIM, ctx);
    gr_mat_init(C, DIM, DIM, ctx);
    randmat(A, N, state, ctx);

    TIMEIT_START;
    sink ^= run_gr(C, A, exp, REPS, ctx);
    TIMEIT_STOP_VALUES(tcpu, t_gr);

#ifdef __AVX2__
    TIMEIT_START;
    sink ^= run_batched(C, A, exp, REPS, ctx);
    TIMEIT_STOP_VALUES(tcpu, t_batched);
#endif

    flint_printf("n = %wd bits (odd), A^exp, exp = %wu (%u-bit chain):\n",
                 bits, exp, (unsigned) FLINT_BIT_COUNT(exp));
    flint_printf("  gr_mat_pow_ui (classical): %8.2f us/pow\n", t_gr / REPS / 1e-6);
#ifdef __AVX2__
    flint_printf("  batched_mont             : %8.2f us/pow\n", t_batched / REPS / 1e-6);
    flint_printf("  speedup batched vs classical: %.2fx\n", t_gr / t_batched);
#else
    flint_printf("  (batched_mont skipped: built without AVX2)\n");
#endif
    flint_printf("\n");

    gr_mat_clear(A, ctx);
    gr_mat_clear(C, ctx);
    gr_ctx_clear(ctx);
    fmpz_clear(N);

    (void) tcpu;
}

int main(void)
{
    flint_rand_t state;
    slong sizes[] = { 96, 128, 192, 246 };
    // 39 squarings
    ulong exp = (UWORD(1) << 40) - 1;
    int i;

    flint_rand_init(state);
    for (i = 0; i < 4; i++)
        bench(sizes[i], exp, state);
    flint_rand_clear(state);
    return 0;
}
