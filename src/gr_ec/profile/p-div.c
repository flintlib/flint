/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    What knowing the order of the group buys, and what dividing a point
    costs.

    Columns:

      count      gr_ec_ctx_cardinality on a context that knows nothing:
                 the price of learning the order the hard way
      hint       gr_ec_ctx_set_order with the right answer, on a context
                 that knows nothing: the price of being told instead. Over
                 a field too large to walk this is twenty scalar
                 multiplications against random points, so it should track
                 the mul column and not the count one
      mul raw    n P for a scalar four times as long as the order, with
                 nothing cached
      mul cached the same multiplication once the order is known, so the
                 scalar is reduced first
      div cop    dividing by a small n coprime to the order: one modular
                 inverse and one scalar multiplication, always unique
      div tors   dividing by a small prime factor of the order, which is
                 the case the cheap path cannot do: a division polynomial
                 of degree n^2 and its roots

    The two div columns are the point of the file. They are the same
    operation and differ by orders of magnitude, which is why
    gr_ec_point_div_fmpz works as hard as it does to stay in the first one.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "gr.h"
#include "gr_ec.h"

#define TIME_US(dest, stmt) \
    do { \
        double _tc, _tw; \
        TIMEIT_START \
        stmt; \
        TIMEIT_STOP_VALUES(_tc, _tw); \
        (dest) = _tw * 1e6; \
        (void) _tc; \
    } while (0)

static int
cmp_double(const void * a, const void * b)
{
    double x = *(const double *) a, y = *(const double *) b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static double
median(double * v, slong n)
{
    if (n == 0)
        return 0.0;
    qsort(v, n, sizeof(double), cmp_double);
    return v[n / 2];
}

/* the smallest prime factor of N below the cap, or 0 */
static ulong
small_factor(const fmpz_t N, ulong cap)
{
    ulong p;

    for (p = 2; p <= cap; p = n_nextprime(p, 1))
        if (fmpz_divisible_si(N, (slong) p))
            return p;

    return 0;
}

#define MAX_CURVES 9

static void
run_size(flint_rand_t state, slong bits, slong ncurves)
{
    double tc[MAX_CURVES], th[MAX_CURVES], mr[MAX_CURVES], mc[MAX_CURVES];
    double dc[MAX_CURVES], dt[MAX_CURVES];
    slong nc = 0, nh = 0, nmr = 0, nmc = 0, ndc = 0, ndt = 0, i;
    fmpz_t p;

    fmpz_init(p);

    for (i = 0; i < ncurves && i < MAX_CURVES; i++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_point_t P, Q;
        gr_ptr a4, a6;
        fmpz_t N, big, n;
        slong tries;
        int built = 0;

        fmpz_randprime(p, state, bits, 0);

        if (bits < 62)
        {
            if (!fmpz_abs_fits_ui(p)
                    || gr_ctx_init_nmod(R, fmpz_get_ui(p)) != GR_SUCCESS)
                continue;
        }
        else
            gr_ctx_init_fmpz_mod(R, p);

        GR_IGNORE(gr_ctx_set_is_field(R, T_TRUE));

        GR_TMP_INIT2(a4, a6, R);

        for (tries = 0; tries < 20 && !built; tries++)
            if (gr_randtest(a4, state, R) == GR_SUCCESS
                    && gr_randtest(a6, state, R) == GR_SUCCESS
                    && gr_ec_ctx_init_short_weierstrass(E, R, a4, a6) == GR_SUCCESS)
                built = 1;

        GR_TMP_CLEAR2(a4, a6, R);

        if (!built)
        {
            gr_ctx_clear(R);
            continue;
        }

        fmpz_init(N);
        fmpz_init(big);
        fmpz_init(n);
        gr_ec_point_init(P, E);
        gr_ec_point_init(Q, E);

        if (gr_ec_ctx_cardinality(N, E) != GR_SUCCESS)
            goto next;

        /* the identity divides by anything at no cost, which would make
           the division columns meaningless */
        for (tries = 0; tries < 50; tries++)
            if (gr_ec_point_randtest(P, state, E) == GR_SUCCESS
                    && gr_ec_point_is_inf(P, E) == T_FALSE)
                break;

        if (tries == 50)
            goto next;

        /* counting, on a context that has been told nothing */
        gr_ec_ctx_clear_order(E);
        TIME_US(tc[nc], { fmpz_t t; fmpz_init(t);
                GR_IGNORE(gr_ec_ctx_cardinality(t, E)); fmpz_clear(t); });
        nc++;

        /* being told, on a context that has been told nothing */
        gr_ec_ctx_clear_order(E);
        TIME_US(th[nh], { gr_ec_ctx_clear_order(E);
                GR_IGNORE(gr_ec_ctx_set_order(E, N)); });
        nh++;

        /* a random scalar four times as long as the order, so that
           reducing it modulo the order genuinely shortens it (a multiple
           of N plus a constant would reduce to the constant and flatter
           the cached column out of all recognition) */
        fmpz_randbits(big, state, 4 * bits);
        fmpz_abs(big, big);

        gr_ec_ctx_clear_order(E);
        TIME_US(mr[nmr], { gr_ec_ctx_clear_order(E);
                GR_IGNORE(gr_ec_point_mul_fmpz(Q, P, big, E)); });
        nmr++;

        GR_IGNORE(gr_ec_ctx_set_order(E, N));
        TIME_US(mc[nmc], { GR_IGNORE(gr_ec_point_mul_fmpz(Q, P, big, E)); });
        nmc++;

        /* dividing by something coprime to the order */
        {
            ulong d;
            fmpz_t g;
            fmpz_init(g);

            for (d = 3; d < 200; d = n_nextprime(d, 1))
            {
                fmpz_set_ui(n, d);
                fmpz_gcd(g, n, N);

                if (fmpz_is_one(g))
                    break;
            }

            fmpz_clear(g);

            if (d < 200 && gr_ec_point_div_fmpz(Q, P, n, E) == GR_SUCCESS)
            {
                TIME_US(dc[ndc], { GR_IGNORE(gr_ec_point_div_fmpz(Q, P, n, E)); });
                ndc++;
            }
        }

        /* and by a small prime that does divide it */
        {
            ulong d = small_factor(N, 40);

            if (d != 0)
            {
                fmpz_set_ui(n, d);

                if (gr_ec_point_div_fmpz_nonunique(Q, P, n, E) == GR_SUCCESS)
                {
                    TIME_US(dt[ndt], {
                        GR_IGNORE(gr_ec_point_div_fmpz_nonunique(Q, P, n, E)); });
                    ndt++;
                }
            }
        }

next:
        gr_ec_point_clear(Q, E);
        gr_ec_point_clear(P, E);
        fmpz_clear(N);
        fmpz_clear(big);
        fmpz_clear(n);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);
    }

    flint_printf("%5wd ", bits);

    if (nc)  flint_printf("%12.1f ", median(tc, nc));  else flint_printf("%12s ", "-");
    if (nh)  flint_printf("%10.1f ", median(th, nh));  else flint_printf("%10s ", "-");
    if (nmr) flint_printf("%10.1f ", median(mr, nmr)); else flint_printf("%10s ", "-");
    if (nmc) flint_printf("%11.1f ", median(mc, nmc)); else flint_printf("%11s ", "-");
    if (ndc) flint_printf("%10.1f ", median(dc, ndc)); else flint_printf("%10s ", "-");
    if (ndt) flint_printf("%10.1f ", median(dt, ndt)); else flint_printf("%10s ", "-");

    flint_printf("\n");

    fmpz_clear(p);
}

int main(void)
{
    flint_rand_t state;

    flint_rand_init(state);

    flint_printf("gr_ec: the order cache, and division by a scalar\n");
    flint_printf("microseconds, median over the curves of each size\n\n");

    flint_printf("%5s %12s %10s %10s %11s %10s %10s\n",
            "bits", "count", "hint", "mul raw", "mul cached", "div cop", "div tors");

    run_size(state, 16, 9);
    run_size(state, 24, 9);
    run_size(state, 32, 7);
    run_size(state, 40, 5);
    run_size(state, 48, 5);
    run_size(state, 56, 5);
    run_size(state, 64, 3);

    flint_printf("\n");
    flint_printf("hint against count is what gr_ec_ctx_set_order saves a caller who\n");
    flint_printf("already knows the answer; mul cached against mul raw is what the\n");
    flint_printf("context saves itself once it does.\n");

    flint_rand_clear(state);

    return 0;
}
