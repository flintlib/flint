/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "flint/profiler.h"

#define TIMEIT_END_REPEAT3(__timer, __reps, __min_time) \
            } \
            timeit_stop(__timer); \
            if (__timer->cpu >= __min_time) \
                break; \
            __reps *= 10; \
        } \
    } while (0);

#define TIMEIT_STOP_VALUES3(tcpu, twall, __min_time) \
        TIMEIT_END_REPEAT3(__timer, __reps, __min_time) \
        (tcpu) = __timer->cpu*0.001 / __reps; \
        (twall) = __timer->wall*0.001 / __reps; \
    } while (0);

#if 0
#define CASE_A GR_IGNORE(gr_poly_inv_series_basecase(B, A, len, ctx));
#define CASE_B GR_IGNORE(gr_poly_inv_series_newton(B, A, len, len, ctx));
#endif

#if 0
#define CASE_A GR_IGNORE(gr_poly_rsqrt_series_basecase(B, A, len, ctx));
#define CASE_B GR_IGNORE(gr_poly_rsqrt_series_newton(B, A, len, len, ctx));
#endif

#if 1
#define CASE_A GR_IGNORE(gr_poly_sqrt_series_basecase(B, A, len, ctx));
#define CASE_B GR_IGNORE(gr_poly_sqrt_series_newton(B, A, len, len, ctx));
#endif

double
get_profile(gr_ctx_t ctx, slong len)
{
    gr_poly_t A, B;
    slong i;
    double tcpu, twall, tbase, tnew;
    flint_rand_t state;

    gr_poly_init(A, ctx);
    gr_poly_init(B, ctx);

    flint_randinit(state);

    for (i = 0; i < len; i++)
        GR_IGNORE(gr_poly_set_coeff_si(A, i, n_randlimb(state), ctx));
    GR_IGNORE(gr_poly_set_coeff_si(A, 0, 1, ctx));

    TIMEIT_START
    CASE_A
    TIMEIT_STOP_VALUES3(tcpu, twall, 10.0)
    (void) tcpu;
    tbase = twall;

    TIMEIT_START
    CASE_B
    TIMEIT_STOP_VALUES3(tcpu, twall, 10.0)
    (void) tcpu;
    tnew = twall;

    printf("len %ld : %f\n", len, tbase / tnew);

    flint_randclear(state);

    gr_poly_clear(A, ctx);
    gr_poly_clear(B, ctx);

    return tbase / tnew;
}

int ok(double x)
{
    return x >= 0.9 && x <= 1.1;
}

slong
get_tuning(gr_ctx_t ctx, slong from)
{
    double speedup;
    slong cutoff = 0, len, consecutive = 0;

    do
    {
        for (len = from; len <= 32767; len = FLINT_MAX(len+1, len*1.01))
        {
            speedup = get_profile(ctx, len);

            if (speedup > 1.0)
            {
                consecutive++;

                if (consecutive == 1)
                    cutoff = len;

                if (consecutive == 3)
                    break;
            }
            else
            {
                consecutive = 0;
            }
        }
    }
    while (!ok(get_profile(ctx, len)));

    return cutoff;
}

int main()
{
    gr_ctx_t ctx;
    slong i, results[64];
    slong bits, cutoff, prev_cutoff = 0;

    for (bits = 2; bits <= 64; bits++)
    {
        gr_ctx_init_nmod(ctx, n_nextprime(UWORD(1) << (bits - 1), 0));
        cutoff = get_tuning(ctx, FLINT_MAX(prev_cutoff * 0.75, 2));
        results[bits - 1] = cutoff;
        prev_cutoff = cutoff;
        flint_printf("bits = %wd  cutoff = %wd  accuracy = %f\n", bits, cutoff, get_profile(ctx, cutoff));
        flint_printf("tab[] = {");
        for (i = 0; i < bits; i++)
            flint_printf("%wd, ", results[i]);
        flint_printf("};\n");
    }

    return EXIT_SUCCESS;
}
