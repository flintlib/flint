/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Scalar multiplication k P on y^2 = x^3 + a x + b over F_p, for a scalar
    of the same size as p. Compares the three point representations of
    gr_ec against each other and against ecpp_point_mul_gr, which uses the
    same modular arithmetic and the same Jacobian formulas but a width-4
    NAF window instead of a plain binary ladder.

    The unit is microseconds per k P, so the numbers can be compared
    directly with what other libraries report at the same modulus size.
    For reference, PARI 2.15 ellmul over the same kind of curve, measured
    with

        S = [32,64,128,192,256,384,512,1024,2048];
        for(idx=1, length(S), bb = S[idx]; pp = randomprime([2^(bb-1), 2^bb]);
            aa = random(pp); xx = random(pp); yy = random(pp);
            cc = lift(Mod(yy^2 - xx^3 - aa*xx, pp)); EE = ellinit([aa, cc], pp);
            PP = [Mod(xx,pp), Mod(yy,pp)]; kk = random(2^bb); nn = 1;
            while(1, gettime(); for(j = 1, nn, ellmul(EE, PP, kk)); tt = gettime();
                  if(tt >= 200, break); nn = nn * 4);
            printf("%5d %10.2f\n", bb, tt*1000.0/nn));

    gave, on the machine this was written on (in microseconds, median of
    three runs, next to the numbers this profile reports there):

        bits      32     64    128    192    256    384    512   1024    2048
        PARI   21.18  40.41 133.54 252.93 372.07 718.75 1296.9 6125.0 36562.5
        ecpp    7.35  25.40  24.70  91.20 137.00 332.00  591.0 3700.0 26800.0
        gr_ec   6.89  26.10  22.10  84.30 132.00 299.00  540.0 3280.0 24000.0

    The jump between 64 and 128 bits is the modulus crossing into mpn_mod,
    which handles 2 to 16 limbs; outside that range the ring is fmpz_mod.
*/

#include <stdio.h>
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "gr.h"
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

int main(void)
{
    const slong bitsizes[] = { 32, 64, 128, 192, 256, 384, 512, 1024, 2048 };
    const slong num_sizes = sizeof(bitsizes) / sizeof(slong);
    flint_rand_t state;
    slong bi;

    flint_rand_init(state);

    flint_printf("gr_ec: scalar multiplication k P over F_p, k of the size of p\n");
    flint_printf("microseconds per k P (lower is better)\n\n");
    flint_printf("%9s %10s %10s %10s %10s %12s\n",
            "bits(p)", "proj", "affine", "jacobian", "jac-binary", "ecpp-naf4");

    for (bi = 0; bi < num_sizes; bi++)
    {
        slong bits = bitsizes[bi];
        fmpz_t p, k;
        fmpz_mod_ctx_t mod;
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_point_t P, A;
        gr_ec_aff_point_t Pa, Aa;
        gr_ec_jac_point_t Pj, Aj;
        gr_ptr a4, a6, x, y, t;
        gr_ptr gP, gR, ga, gacc;
        double t_proj, t_aff, t_jac, t_jacbin, t_ecpp;
        slong sz;
        int status;

        fmpz_init(p);
        fmpz_init(k);
        fmpz_randprime(p, state, bits, 0);
        fmpz_mod_ctx_init(mod, p);

        /* the same ring ecpp would pick: mpn_mod when it applies, else fmpz_mod */
        ecpp_gr_ctx_init(R, mod);
        sz = R->sizeof_elem;

        GR_TMP_INIT5(a4, a6, x, y, t, R);

        /* a random point (x, y) on a random curve: pick a4, x, y and solve for a6 */
        status = gr_randtest(a4, state, R);
        status |= gr_randtest(x, state, R);
        status |= gr_randtest(y, state, R);

        /* a6 = y^2 - x^3 - a4 x */
        status |= gr_sqr(a6, y, R);
        status |= gr_sqr(t, x, R);
        status |= gr_mul(t, t, x, R);
        status |= gr_sub(a6, a6, t, R);
        status |= gr_mul(t, a4, x, R);
        status |= gr_sub(a6, a6, t, R);

        if (status != GR_SUCCESS
                || gr_ec_ctx_init_short_weierstrass(E, R, a4, a6) != GR_SUCCESS)
        {
            flint_printf("%9wd  (singular or unusable curve, skipped)\n", bits);
            GR_TMP_CLEAR5(a4, a6, x, y, t, R);
            gr_ctx_clear(R);
            fmpz_mod_ctx_clear(mod);
            fmpz_clear(p);
            fmpz_clear(k);
            continue;
        }

        gr_ec_point_init(P, E);
        gr_ec_point_init(A, E);
        gr_ec_aff_point_init(Pa, E);
        gr_ec_aff_point_init(Aa, E);
        gr_ec_jac_point_init(Pj, E);
        gr_ec_jac_point_init(Aj, E);

        GR_MUST_SUCCEED(gr_ec_point_set_affine(P, x, y, E));
        GR_MUST_SUCCEED(gr_ec_aff_point_set_affine(Pa, x, y, E));
        GR_MUST_SUCCEED(gr_ec_jac_point_set_affine(Pj, x, y, E));

        fmpz_randbits(k, state, bits);
        fmpz_abs(k, k);

        /* the ecpp reference works on a raw (X, Y, Z) triple with Z = 1 */
        gP = gr_heap_init_vec(3, R);
        gR = gr_heap_init_vec(3, R);
        ga = gr_heap_init(R);
        gacc = gr_heap_init(R);
        GR_MUST_SUCCEED(gr_set(gP, x, R));
        GR_MUST_SUCCEED(gr_set(GR_ENTRY(gP, 1, sz), y, R));
        GR_MUST_SUCCEED(gr_one(GR_ENTRY(gP, 2, sz), R));
        GR_MUST_SUCCEED(gr_set(ga, a4, R));
        GR_MUST_SUCCEED(gr_one(gacc, R));

        /* check that all five compute the same point before timing them */
        GR_MUST_SUCCEED(gr_ec_point_mul_fmpz(A, P, k, E));
        GR_MUST_SUCCEED(gr_ec_aff_point_mul_fmpz(Aa, Pa, k, E));
        GR_MUST_SUCCEED(gr_ec_jac_point_mul_fmpz(Aj, Pj, k, E));

        {
            gr_ec_point_t B;
            gr_ec_jac_point_t Bj;

            gr_ec_point_init(B, E);
            gr_ec_jac_point_init(Bj, E);

            GR_MUST_SUCCEED(gr_ec_point_set_aff_point(B, Aa, E));

            if (gr_ec_point_equal(A, B, E) != T_TRUE)
                flint_throw(FLINT_ERROR, "affine result differs at %wd bits\n", bits);

            GR_MUST_SUCCEED(gr_ec_point_set_jac_point(B, Aj, E));

            if (gr_ec_point_equal(A, B, E) != T_TRUE)
                flint_throw(FLINT_ERROR, "jacobian result differs at %wd bits\n", bits);

            GR_MUST_SUCCEED(_gr_ec_jac_point_mul_fmpz_binary(Bj, Pj, k, E));
            GR_MUST_SUCCEED(gr_ec_point_set_jac_point(B, Bj, E));

            if (gr_ec_point_equal(A, B, E) != T_TRUE)
                flint_throw(FLINT_ERROR, "jacobian ladder differs at %wd bits\n", bits);

            if (!ecpp_point_mul_gr(gR, gP, k, ga, gacc, R))
                flint_throw(FLINT_ERROR, "ecpp reference failed at %wd bits\n", bits);

            GR_MUST_SUCCEED(_gr_ec_jac_point_set_jacobian(Bj, gR,
                        GR_ENTRY(gR, 1, sz), GR_ENTRY(gR, 2, sz), E));
            GR_MUST_SUCCEED(gr_ec_point_set_jac_point(B, Bj, E));

            if (gr_ec_point_equal(A, B, E) != T_TRUE)
                flint_throw(FLINT_ERROR, "ecpp result differs at %wd bits\n", bits);

            gr_ec_jac_point_clear(Bj, E);
            gr_ec_point_clear(B, E);
        }

        TIME_US(t_proj, GR_MUST_SUCCEED(gr_ec_point_mul_fmpz(A, P, k, E)));
        TIME_US(t_aff, GR_MUST_SUCCEED(gr_ec_aff_point_mul_fmpz(Aa, Pa, k, E)));
        TIME_US(t_jac, GR_MUST_SUCCEED(gr_ec_jac_point_mul_fmpz(Aj, Pj, k, E)));
        TIME_US(t_jacbin, GR_MUST_SUCCEED(_gr_ec_jac_point_mul_fmpz_binary(Aj, Pj, k, E)));
        TIME_US(t_ecpp, (void) ecpp_point_mul_gr(gR, gP, k, ga, gacc, R));

        flint_printf("%9wd %10.2f %10.2f %10.2f %10.2f %12.2f\n",
                bits, t_proj, t_aff, t_jac, t_jacbin, t_ecpp);
        fflush(stdout);

        gr_heap_clear_vec(gP, 3, R);
        gr_heap_clear_vec(gR, 3, R);
        gr_heap_clear(ga, R);
        gr_heap_clear(gacc, R);
        gr_ec_jac_point_clear(Pj, E);
        gr_ec_jac_point_clear(Aj, E);
        gr_ec_aff_point_clear(Pa, E);
        gr_ec_aff_point_clear(Aa, E);
        gr_ec_point_clear(P, E);
        gr_ec_point_clear(A, E);
        gr_ec_ctx_clear(E);
        GR_TMP_CLEAR5(a4, a6, x, y, t, R);
        gr_ctx_clear(R);
        fmpz_mod_ctx_clear(mod);
        fmpz_clear(p);
        fmpz_clear(k);
    }

    flint_rand_clear(state);
    flint_cleanup_master();

    return 0;
}
