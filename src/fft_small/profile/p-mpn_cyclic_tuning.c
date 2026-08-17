/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Times the powm fold and linear pipeline step for every admissible
   cyclic geometry at several modulus sizes; the selection rule in
   fft_small_plan_init_mpn_cyclic was set from this data. */

#include <stdio.h>
#include "flint.h"
#include "mpn_extras.h"
#include "fft_small.h"
#include "profiler.h"

static void run(mp_size_t mn, flint_rand_t state)
{
    mpn_ctx_struct * R = get_default_mpn_ctx();
    printf("== mn = %ld (%ld bits) ==\n", (slong) mn, (slong) mn * 64);

    for (ulong pi = 0; pi < R->profiles_size; pi++)
    {
        ulong bits = R->profiles[pi].bits, np = R->profiles[pi].np;
        ulong depth, N0, nn, zc, t, thi, nn2;
        fft_small_plan_t P, PL;
        fft_small_op_t Fm, FC, FA, FB, FZ;
        nn_ptr m, q, w, X;
        timeit_t tf, tl;
        slong reps, i;
        double fold_ms, lin_ms;

        if (np < 4) continue;
        depth = n_max(7, n_clog2(n_cdiv(64 * (ulong) mn, bits)));
        N0 = n_pow2(depth);
        nn = (bits << depth) / 64;
        umul_ppmm(thi, t, UWORD(1), N0);
        if (thi || flint_mpn_cmp_ui_2exp(
                crt_data_prod_primes(R->crts + np - 1),
                (R->crts + np - 1)->coeff_len, t, 2*bits) < 0)
            continue;

        P->R = R; P->sign = 0; P->depth = depth;
        P->zl = 0; P->zh = N0; P->zn = N0; P->ztrunc = N0;
        P->bits = bits; P->bound_c = 1; P->bound_e = 2*bits;
        P->np = np; P->offset = 0;
        P->ffts = R->ffts; P->crts = R->crts; P->use_direct_fft = 0;
        _fft_small_plan_set_normalizers(P);

        if (!fft_small_plan_init_mpn(PL, R, nn, nn, 1))
            continue;

        zc = nn + (2*bits + depth)/64 + 2;
        m = flint_malloc(mn * 8);
        q = flint_malloc(nn * 8);
        w = flint_malloc((nn + 8) * 8);
        X = flint_malloc(2 * nn * 8);
        flint_mpn_rrandom(m, state, mn);
        flint_mpn_rrandom(q, state, nn);

        fft_small_op_init(Fm, P);  fft_small_op_init(FC, P);
        fft_small_op_init(FA, PL); fft_small_op_init(FB, PL);
        fft_small_op_init(FZ, PL);
        fft_small_fft_mpn(Fm, m, mn, P);
        fft_small_fft_mpn(FB, q, nn, PL);

        reps = FLINT_MAX(10, 5000000 / (slong) (nn * np / 4));
        timeit_start(tf);
        for (i = 0; i < reps; i++)
        {
            fft_small_fft_mpn(FC, q, nn, P);
            fft_small_op_mul(FC, FC, Fm, P);
            fft_small_ifft(FC, P);
            fft_small_export_mpn(w, zc, FC, P);
        }
        timeit_stop(tf);
        fold_ms = (double) tf->cpu / reps;

        reps = FLINT_MAX(5, 1200000 / (slong) nn);
        timeit_start(tl);
        for (i = 0; i < reps; i++)
        {
            fft_small_fft_mpn(FA, q, nn, PL);
            fft_small_op_mul(FZ, FA, FB, PL);
            fft_small_ifft(FZ, PL);
            fft_small_export_mpn(X, 2*nn, FZ, PL);
        }
        timeit_stop(tl);
        lin_ms = (double) tl->cpu / reps;

        printf("  cyc(np %lu, bits %3lu, d %2lu) nn %6lu pad %5.1f%% | "
               "lin(np %lu, bits %3lu) | fold %7.3f  linstep %7.3f  "
               "total~ %7.3f ms\n",
               np, bits, depth, nn, 100.0*(nn-mn)/mn,
               PL->np, PL->bits, fold_ms, lin_ms, 2*lin_ms + fold_ms);
        fflush(stdout);

        fft_small_op_clear(Fm); fft_small_op_clear(FC);
        fft_small_op_clear(FA); fft_small_op_clear(FB);
        fft_small_op_clear(FZ);
        fft_small_plan_clear(P); fft_small_plan_clear(PL);
        flint_free(m); flint_free(q); flint_free(w); flint_free(X);
    }
}

int main(void)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_set_num_threads(1);
    run(1300, state);
    run(1600, state);
    run(5200, state);
    run(16000, state);
    flint_rand_clear(state);
    return 0;
}
