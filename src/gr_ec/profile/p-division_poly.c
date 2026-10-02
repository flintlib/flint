/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    The n-th division polynomial Psi_n, which has degree about n^2 / 2,
    over a word sized prime field, a 256 bit prime field and Q.

    The last two columns are the point of the exercise. "table" builds
    Psi_0, ..., Psi_n with gr_ec_ctx_division_poly_vec, which is the
    obvious way to run the recursion and costs O(n M(n^2)); the ladder in
    gr_ec_ctx_division_poly evaluates only the transitive closure of the
    recursion, about ten indices per halving of n, and costs O(M(n^2)).
    The ratio is what that buys.

    For reference, PARI 2.15 elldivpol, measured with

        default(parisizemax, 8*10^9);
        NS = [10,20,50,100,200,400];
        pp = randomprime([2^61, 2^62]); aa = random(pp); bb = random(pp);
        EE = ellinit([aa,bb], pp);
        for(i=1, length(NS), nn = NS[i]; rr = 1;
            while(1, gettime(); for(j=1,rr, elldivpol(EE,nn)); tt = gettime();
                  if(tt >= 300, break); rr = rr*4);
            printf("%5d %12.1f\n", nn, tt*1000.0/rr));

    gave, on the machine this was written on (microseconds, next to the
    numbers this profile reports there):

        n                      10      20       50       100       200       400
        PARI  F_p 2^62       23.5   133.3   1519.5    8578.1   44437.5  209000.0
        gr_ec F_p 2^62        4.6    36.8    325.0    1360.0    5580.0   24000.0

        PARI  F_p 2^256     101.6   721.7   8562.5   43875.0  189750.0  915000.0
        gr_ec F_p 2^256      43.7   369.0   3230.0   14000.0   66300.0  356000.0

        PARI  over Q         15.5   195.8  14500.0  365000.0        --        --
        gr_ec over Q          9.0   144.0   5720.0        --        --        --

    The 256 bit column runs on mpn_mod, which is what ecpp_gr_ctx_init
    selects at four limbs.

    Note that PARI's elldivpol uses the opposite normalisation for even n:
    it returns psi_n * psi_2 = Psi_n * psi_2^2, of degree (n^2 + 2) / 2,
    where we return Psi_n = psi_n / psi_2, of degree (n^2 - 4) / 2. The two
    differ by the cubic gr_ec_ctx_psi2_sqr, so for even n the timings above
    are for a slightly larger object than ours.
*/

#include <stdio.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_ec.h"
#include "ecpp.h"

#define TIME_US(dest, stmt) \
    do { \
        double _tc, _tw; \
        TIMEIT_START \
        stmt; \
        TIMEIT_STOP_VALUES(_tc, _tw); \
        (dest) = _tw * 1e6; \
        (void) _tc; \
    } while (0)

/* a nonsingular short Weierstrass curve over R */
static int
random_curve(gr_ec_ctx_t E, gr_ctx_t R, flint_rand_t state)
{
    gr_ptr a4, a6;
    int status;
    slong iter;

    GR_TMP_INIT2(a4, a6, R);

    status = GR_UNABLE;

    for (iter = 0; iter < 20; iter++)
    {
        if (gr_randtest(a4, state, R) != GR_SUCCESS
                || gr_randtest(a6, state, R) != GR_SUCCESS)
            continue;

        status = gr_ec_ctx_init_short_weierstrass(E, R, a4, a6);

        if (status == GR_SUCCESS)
            break;
    }

    GR_TMP_CLEAR2(a4, a6, R);

    return status;
}

static double
time_ladder(gr_ec_ctx_t E, ulong n)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(E);
    gr_poly_t p;
    double t;

    gr_poly_init(p, R);
    TIME_US(t, { GR_IGNORE(gr_ec_ctx_division_poly(p, n, E)); });
    gr_poly_clear(p, R);

    return t;
}

static double
time_table(gr_ec_ctx_t E, ulong n)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(E);
    gr_poly_struct * tab;
    slong len = (slong) n + 1, k;
    double t;

    tab = flint_malloc(len * sizeof(gr_poly_struct));
    for (k = 0; k < len; k++)
        gr_poly_init(tab + k, R);

    TIME_US(t, { GR_IGNORE(gr_ec_ctx_division_poly_vec(tab, len, E)); });

    for (k = 0; k < len; k++)
        gr_poly_clear(tab + k, R);
    flint_free(tab);

    return t;
}

int main(void)
{
    const ulong ns[] = { 10, 20, 50, 100, 200, 400 };
    const slong num_ns = sizeof(ns) / sizeof(ulong);
    flint_rand_t state;
    gr_ctx_t Rword, Rbig, Rq;
    gr_ec_ctx_t Eword, Ebig, Eq;
    fmpz_t p;
    fmpz_mod_ctx_t mod;
    slong i;
    int have_word, have_big, have_q;

    flint_rand_init(state);

    /* word sized prime field */
    GR_MUST_SUCCEED(gr_ctx_init_nmod(Rword, n_randprime(state, 62, 1)));
    have_word = (random_curve(Eword, Rword, state) == GR_SUCCESS);

    /* 256 bit prime field, on whichever ring FLINT prefers at that size */
    fmpz_init(p);
    fmpz_randprime(p, state, 256, 0);
    fmpz_mod_ctx_init(mod, p);
    ecpp_gr_ctx_init(Rbig, mod);
    have_big = (random_curve(Ebig, Rbig, state) == GR_SUCCESS);

    /* over Q */
    gr_ctx_init_fmpq(Rq);
    have_q = (gr_ec_ctx_init_si(Eq, Rq, 1, 2, 3, 4, 5) == GR_SUCCESS);

    flint_printf("gr_ec: n-th division polynomial Psi_n (degree about n^2/2)\n");
    flint_printf("microseconds (lower is better)\n\n");
    flint_printf("%6s %12s %12s %12s %12s %8s\n",
            "n", "F_p 2^62", "F_p 2^256", "over Q", "table/F_p", "ratio");

    for (i = 0; i < num_ns; i++)
    {
        ulong n = ns[i];
        double tw = 0, tb = 0, tq = 0, tt = 0;

        if (have_word)
        {
            tw = time_ladder(Eword, n);
            tt = time_table(Eword, n);
        }

        if (have_big)
            tb = time_ladder(Ebig, n);

        /* the coefficients over Q grow like n^2 digits, so stop early */
        if (have_q && n <= 50)
            tq = time_ladder(Eq, n);

        flint_printf("%6wu %12.1f %12.1f ", n, tw, tb);

        if (tq != 0)
            flint_printf("%12.1f ", tq);
        else
            flint_printf("%12s ", "--");

        flint_printf("%12.1f %8.1f\n", tt, (tw != 0) ? tt / tw : 0.0);
    }

    if (have_q)
        gr_ec_ctx_clear(Eq);
    gr_ctx_clear(Rq);

    if (have_big)
        gr_ec_ctx_clear(Ebig);
    gr_ctx_clear(Rbig);
    fmpz_mod_ctx_clear(mod);
    fmpz_clear(p);

    if (have_word)
        gr_ec_ctx_clear(Eword);
    gr_ctx_clear(Rword);

    flint_rand_clear(state);

    return 0;
}
