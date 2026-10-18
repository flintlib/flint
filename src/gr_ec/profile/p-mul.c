/*
    Copyright (C) 2026 Maël Hostettler
    Developed using Claude Opus 5 and Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Scalar multiplication k P on y^2 = x^3 + a x + b over F_p, with k of
    the size of p: affine, Jacobian, the plain binary Jacobian ladder and
    the witness ladder used by ECPP. Prints microseconds per k P.
*/

#include <stdio.h>
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "gr.h"
#include "gr_ec.h"
#include "gr_ec/impl.h"
#include "mpn_mod.h"

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
    flint_printf("%9s %10s %10s %10s %12s\n",
            "bits(p)", "affine", "jacobian", "jac-binary", "jac-witness");

    for (bi = 0; bi < num_sizes; bi++)
    {
        slong bits = bitsizes[bi];
        fmpz_t p, k;
        fmpz_mod_ctx_t mod;
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_jac_point_t P, A;
        gr_ec_aff_point_t Pa, Aa;
        gr_ec_jac_point_t Pj, Aj;
        gr_ptr a4, a6, x, y, t;
        gr_ptr gacc;
        double t_aff, t_jac, t_jacbin, t_wit;
        int status;

        fmpz_init(p);
        fmpz_init(k);
        fmpz_randprime(p, state, bits, 0);
        fmpz_mod_ctx_init(mod, p);

        /* mpn_mod when it applies, else fmpz_mod */
        if (fmpz_size(p) < MPN_MOD_MIN_LIMBS || fmpz_size(p) > MPN_MOD_MAX_LIMBS
                || gr_ctx_init_mpn_mod(R, p) != GR_SUCCESS)
            gr_ctx_init_fmpz_mod(R, p);
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

        gr_ec_jac_point_init(P, E);
        gr_ec_jac_point_init(A, E);
        gr_ec_aff_point_init(Pa, E);
        gr_ec_aff_point_init(Aa, E);
        gr_ec_jac_point_init(Pj, E);
        gr_ec_jac_point_init(Aj, E);

        GR_MUST_SUCCEED(gr_ec_jac_point_set_affine(P, x, y, E));
        GR_MUST_SUCCEED(gr_ec_aff_point_set_affine(Pa, x, y, E));
        GR_MUST_SUCCEED(gr_ec_jac_point_set_affine(Pj, x, y, E));

        fmpz_randbits(k, state, bits);
        fmpz_abs(k, k);

        /* the ecpp reference works on a raw (X, Y, Z) triple with Z = 1 */
        gacc = gr_heap_init(R);

        /* check that all four compute the same point before timing them */
        GR_MUST_SUCCEED(gr_ec_jac_point_mul_fmpz(A, P, k, E));
        GR_MUST_SUCCEED(gr_ec_aff_point_mul_fmpz(Aa, Pa, k, E));
        GR_MUST_SUCCEED(gr_ec_jac_point_mul_fmpz(Aj, Pj, k, E));

        {
            gr_ec_jac_point_t B;
            gr_ec_jac_point_t Bj;

            gr_ec_jac_point_init(B, E);
            gr_ec_jac_point_init(Bj, E);

            GR_MUST_SUCCEED(gr_ec_jac_point_set_aff_point(B, Aa, E));

            if (gr_ec_jac_point_equal(A, B, E) != T_TRUE)
                flint_throw(FLINT_ERROR, "affine result differs at %wd bits\n", bits);

            GR_MUST_SUCCEED(gr_ec_jac_point_set(B, Aj, E));

            if (gr_ec_jac_point_equal(A, B, E) != T_TRUE)
                flint_throw(FLINT_ERROR, "jacobian result differs at %wd bits\n", bits);

            GR_MUST_SUCCEED(_gr_ec_jac_point_mul_fmpz_binary(Bj, Pj, k, E));
            GR_MUST_SUCCEED(gr_ec_jac_point_set(B, Bj, E));

            if (gr_ec_jac_point_equal(A, B, E) != T_TRUE)
                flint_throw(FLINT_ERROR, "jacobian ladder differs at %wd bits\n", bits);

            /* the witness ladder has its own fallback, so it is worth
               checking that it lands on the same point */
            GR_MUST_SUCCEED(gr_one(gacc, R));
            GR_MUST_SUCCEED(gr_ec_jac_point_mul_fmpz_witness(Bj, Pj, k, gacc, E));
            GR_MUST_SUCCEED(gr_ec_jac_point_set(B, Bj, E));

            if (gr_ec_jac_point_equal(A, B, E) != T_TRUE)
                flint_throw(FLINT_ERROR, "witness result differs at %wd bits\n", bits);

            /* over a prime modulus no branch can have been ambiguous */
            if (gr_is_zero(gacc, R) == T_TRUE)
                flint_throw(FLINT_ERROR, "witness vanished at %wd bits\n", bits);

            gr_ec_jac_point_clear(Bj, E);
            gr_ec_jac_point_clear(B, E);
        }

        TIME_US(t_aff, GR_MUST_SUCCEED(gr_ec_aff_point_mul_fmpz(Aa, Pa, k, E)));
        TIME_US(t_jac, GR_MUST_SUCCEED(gr_ec_jac_point_mul_fmpz(Aj, Pj, k, E)));
        TIME_US(t_jacbin, GR_MUST_SUCCEED(_gr_ec_jac_point_mul_fmpz_binary(Aj, Pj, k, E)));
        TIME_US(t_wit, { GR_MUST_SUCCEED(gr_one(gacc, R));
                GR_MUST_SUCCEED(gr_ec_jac_point_mul_fmpz_witness(Aj, Pj, k, gacc, E)); });

        flint_printf("%9wd %10.2f %10.2f %10.2f %12.2f\n",
                bits, t_aff, t_jac, t_jacbin, t_wit);
        fflush(stdout);

        gr_heap_clear(gacc, R);
        gr_ec_jac_point_clear(Pj, E);
        gr_ec_jac_point_clear(Aj, E);
        gr_ec_aff_point_clear(Pa, E);
        gr_ec_aff_point_clear(Aa, E);
        gr_ec_jac_point_clear(P, E);
        gr_ec_jac_point_clear(A, E);
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
