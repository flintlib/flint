/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "nmod_poly.h"
#include "fq.h"
#include "fq_nmod.h"
#include "fq_zech.h"
#include "profiler.h"

slong next_powhalf2(slong n)
{
    if ((n & (n - 1)) == 0)
        return n * 1.414213562373095 + 0.5;
    else
        return WORD(1) << FLINT_BIT_COUNT(n);
}

void _nmod_poly_mul_mid_default_mpn_ctx(mp_ptr res, slong zl, slong zh, mp_srcptr a, slong an, mp_srcptr b, slong bn, nmod_t mod);

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
#define INIT_CTX gr_ctx_init_nmod(ctx, n_nextprime(UWORD(1) << (bits - 1), 0));
#define RANDCOEFF(t, ctx) GR_IGNORE(gr_set_ui(t, n_randlimb(state), ctx))
#define STEP_BITS for (bits = 1, j = 0; bits <= 64; bits++, j++)
#endif

#if 1
#define INIT_CTX fmpz_t t; fmpz_init(t); fmpz_ui_pow_ui(t, 2, bits - 1); fmpz_add_ui(t, t, 1); /* fmpz_nextprime(t, t, 0); */ gr_ctx_init_fmpz_mod(ctx, t); fmpz_clear(t);
#define RANDCOEFF(t, ctx) fmpz_mod_rand(t, state, gr_ctx_data_as_ptr(ctx));
#define STEP_BITS for (bits = 32, j = 0; bits <= 65536; bits = next_powhalf2(bits), j++)
#endif

#if 0
#define INIT_CTX fmpz_t t; fmpz_init(t); fmpz_ui_pow_ui(t, 2, bits - 1); fmpz_nextprime(t, t, 0); gr_ctx_init_fq(ctx, t, 2, "a"); fmpz_clear(t);
#define RANDCOEFF(t, ctx) fq_rand(t, state, gr_ctx_data_as_ptr(ctx));
#define STEP_BITS for (bits = 100, j = 1; bits <= 100; bits *= 2, j++)
#endif

#if 0
#define INIT_CTX fmpz_t t; fmpz_init(t); fmpz_set_ui(t, n_nextprime(1000000, 0)); gr_ctx_init_fq_nmod(ctx, t, bits, "a"); fmpz_clear(t);
#define RANDCOEFF(t, ctx) fq_nmod_rand(t, state, gr_ctx_data_as_ptr(ctx));
#define STEP_BITS for (bits = 4, j = 1; bits <= 64; bits *= 2, j++)
#endif

#if 0
#define INIT_CTX fmpz_t t; fmpz_init(t); fmpz_set_ui(t, 3); gr_ctx_init_fq_zech(ctx, t, bits, "a"); fmpz_clear(t);
#define RANDCOEFF(t, ctx) fq_zech_rand(t, state, gr_ctx_data_as_ptr(ctx));
#define STEP_BITS for (bits = 2, j = 1; bits <= 12; bits++, j++)
#endif


#if 0
#define INFO "inv_series"
#define SETUP random_input(A, state, len, ctx); \
              GR_IGNORE(gr_poly_set_coeff_si(A, 0, 1, ctx));
#define CASE_A GR_IGNORE(gr_poly_inv_series_basecase(B, A, len, ctx));
#define CASE_B GR_IGNORE(gr_poly_inv_series_newton(B, A, len, len, ctx));
#endif

#if 0
#define INFO "div_series"
#define SETUP random_input(A, state, len, ctx); \
              random_input(B, state, len, ctx); \
              GR_IGNORE(gr_poly_set_coeff_si(B, 0, 1, ctx));
#define CASE_A GR_IGNORE(gr_poly_div_series_basecase(C, A, B, len, ctx));
#define CASE_B GR_IGNORE(gr_poly_div_series_newton(C, A, B, len, len, ctx));
#endif

#if 0
#define INFO "rsqrt_series"
#define SETUP random_input(A, state, len, ctx); \
              GR_IGNORE(gr_poly_set_coeff_si(A, 0, 1, ctx));
#define CASE_A GR_IGNORE(gr_poly_rsqrt_series_basecase(B, A, len, ctx));
#define CASE_B GR_IGNORE(gr_poly_rsqrt_series_newton(B, A, len, len, ctx));
#endif

#if 0
#define INFO "sqrt_series"
#define SETUP random_input(A, state, len, ctx); \
              GR_IGNORE(gr_poly_set_coeff_si(A, 0, 1, ctx));
#define CASE_A GR_IGNORE(gr_poly_sqrt_series_basecase(B, A, len, ctx));
#define CASE_B GR_IGNORE(gr_poly_sqrt_series_newton(B, A, len, len, ctx));
#endif

#if 0
#define INFO "exp_series (basecase -> mul)"
#define SETUP random_input(A, state, len, ctx); \
              GR_IGNORE(gr_poly_set_coeff_si(A, 0, 0, ctx));
#define CASE_A GR_IGNORE(gr_poly_exp_series_basecase(B, A, len, ctx));
#define CASE_B GR_IGNORE(gr_poly_exp_series_basecase_mul(B, A, len, ctx));
#endif

#if 0
#define INFO "exp_series (mul-> newton)"
#define SETUP random_input(A, state, len, ctx); \
              GR_IGNORE(gr_poly_set_coeff_si(A, 0, 0, ctx));
#define CASE_A GR_IGNORE(gr_poly_exp_series_basecase_mul(B, A, len, ctx));
#define CASE_B GR_IGNORE(gr_poly_exp_series_newton(B, A, len, len, ctx));
#endif


#if 0
#define INFO "divrem"
#define SETUP random_input(A, state, 2 * len, ctx); \
              random_input(B, state, len, ctx); \
              GR_IGNORE(gr_poly_set_coeff_si(B, len - 1, 1, ctx));
#define CASE_A GR_IGNORE(gr_poly_divrem_basecase(C, D, A, B, ctx));
#define CASE_B GR_IGNORE(gr_poly_divrem_newton(C, D, A, B, ctx));
#endif

#if 0
#define INFO "divrem (nmod basecase)"
#define SETUP random_input(A, state, 2 * len, ctx); \
              random_input(B, state, len, ctx); \
              GR_IGNORE(gr_poly_set_coeff_si(B, len - 1, 1, ctx));
#define CASE_A gr_poly_fit_length(C, A->length - B->length + 1, ctx); gr_poly_fit_length(D, B->length - 1, ctx); \
               _nmod_poly_divrem(C->coeffs, D->coeffs, A->coeffs, A->length, B->coeffs, B->length, ((nmod_t *) gr_ctx_data_ptr(ctx))[0]);
#define CASE_B GR_IGNORE(gr_poly_divrem_newton(C, D, A, B, ctx));
#endif

#if 0
#define INFO "divrem (fmpz_mod basecase)"
#define SETUP random_input(A, state, 2 * len, ctx); \
              random_input(B, state, len, ctx); \
              GR_IGNORE(gr_poly_set_coeff_si(B, len - 1, 1, ctx));
#define CASE_A GR_IGNORE(gr_poly_set(D, A, ctx)); fmpz_mod_poly_divrem_basecase(C, D, D, B, gr_ctx_data_as_ptr(ctx));
#define CASE_B GR_IGNORE(gr_poly_set(D, A, ctx)); GR_IGNORE(gr_poly_divrem_newton(C, D, D, B, ctx));
#endif

#if 0
#define INFO "mul (nmod)"
#define SETUP random_input(A, state, len, ctx); \
              random_input(B, state, len, ctx);
#define CASE_A gr_poly_fit_length(C, A->length + B->length - 1, ctx); \
               _nmod_poly_mul(C->coeffs, A->coeffs, A->length, B->coeffs, B->length, ((nmod_t *) gr_ctx_data_ptr(ctx))[0]);
#define CASE_B gr_poly_fit_length(C, A->length + B->length - 1, ctx); \
               flint_nmod_poly_mul_mid_fft_small(C->coeffs, 0, A->length + B->length - 1, A->coeffs, A->length, B->coeffs, B->length, ((nmod_t *) gr_ctx_data_ptr(ctx))[0]);
#endif

#if 0
#define INFO "sqr (nmod)"
#define SETUP random_input(A, state, len, ctx); \
              random_input(B, state, len, ctx);
#define CASE_A gr_poly_fit_length(C, A->length + A->length - 1, ctx); \
               _nmod_poly_mul(C->coeffs, A->coeffs, A->length, A->coeffs, A->length, ((nmod_t *) gr_ctx_data_ptr(ctx))[0]);
#define CASE_B gr_poly_fit_length(C, A->length + A->length - 1, ctx); \
               flint_nmod_poly_mul_mid_fft_small(C->coeffs, 0, A->length + A->length - 1, A->coeffs, A->length, A->coeffs, A->length, ((nmod_t *) gr_ctx_data_ptr(ctx))[0]);
#endif

#if 0
#define INFO "mullow (nmod)"
#define SETUP random_input(A, state, len, ctx); \
              random_input(B, state, len, ctx);
#define CASE_A gr_poly_fit_length(C, A->length + B->length - 1, ctx); \
               _nmod_poly_mullow(C->coeffs, A->coeffs, A->length, B->coeffs, B->length, B->length, ((nmod_t *) gr_ctx_data_ptr(ctx))[0]);
#define CASE_B gr_poly_fit_length(C, A->length + B->length - 1, ctx); \
               _nmod_poly_mul_mid_default_mpn_ctx(C->coeffs, 0, B->length, A->coeffs, A->length, B->coeffs, B->length, ((nmod_t *) gr_ctx_data_ptr(ctx))[0]);
#endif


void random_input(gr_poly_t A, flint_rand_t state, slong len, gr_ctx_t ctx)
{
    gr_ptr t;
    slong i;
    GR_TMP_INIT(t, ctx);

    for (i = 0; i < len; i++)
    {
        RANDCOEFF(t, ctx);
        GR_IGNORE(gr_poly_set_coeff_scalar(A, i, t, ctx));
    }

    while (A->length != len)
    {
        RANDCOEFF(t, ctx);
        GR_IGNORE(gr_poly_set_coeff_scalar(A, len - 1, t, ctx));
    }

    GR_TMP_CLEAR(t, ctx);
}

double
get_profile(gr_ctx_t ctx, slong len)
{
    gr_poly_t A, B, C, D;
    double tcpu, twall, tbase, tnew;
    flint_rand_t state;

    gr_poly_init(A, ctx);
    gr_poly_init(B, ctx);
    gr_poly_init(C, ctx);
    gr_poly_init(D, ctx);

    flint_randinit(state);

    SETUP

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
    gr_poly_clear(C, ctx);
    gr_poly_clear(D, ctx);

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
        for (len = from; len <= 32767; len = FLINT_MAX(len+1, len*1.05))
        {
            speedup = get_profile(ctx, len);

            if (speedup > 1.0)
            {
                consecutive++;

                if (consecutive == 1)
                    cutoff = len;

                if (consecutive == 30)
                    break;
            }
            else
            {
                consecutive = 0;
            }
        }
    }
    while (0 /*!ok(get_profile(ctx, len))*/);

    return cutoff;
}

int main()
{
    gr_ctx_t ctx;
    slong i, results[64], j;
    slong bits, cutoff /*, prev_cutoff = 0 */;

    printf("%s\n", INFO);

    STEP_BITS
    {
        INIT_CTX
        /* cutoff = get_tuning(ctx, FLINT_MAX(prev_cutoff * 0.75 - 5, 2)); */
        cutoff = get_tuning(ctx, 2);
        results[j] = cutoff;
        /* prev_cutoff = cutoff; */
        flint_printf("bits = %wd  cutoff = %wd  accuracy = %f\n", bits, cutoff, get_profile(ctx, cutoff));
        flint_printf("tab[] = {");
        for (i = 0; i <= j; i++)
            flint_printf("%wd, ", results[i]);
        flint_printf("};\n");
    }

    return 0;
}
