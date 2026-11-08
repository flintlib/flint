/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Code (claude-opus-5-5) */

#include "test_helpers.h"
#include "acb_types.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"
#include "gr_vec.h"

#define LOOP_STEPS 16

/* Compare the local monodromy matrix at a regular singular point at the origin
   with the transition matrix along the path 0 → a → (loop around 0 in
   LOOP_STEPS steps) → a → 0. */
TEST_GR_FUNCTION_START(acb_ode_local_monodromy, state, count_success, count_domain, count_unable)
{
    gr_ctx_t ZZ;
    gr_ctx_init_fmpz(ZZ);

    for (slong iter = 0; iter < 25 * flint_test_multiplier(); iter++)
    {
        ulong seed1, seed2;
        flint_rand_get_seed(&seed1, &seed2, state);

        int status = GR_SUCCESS;
        slong prec = 30 + n_randint(state, 100);

        gr_ctx_t Pol, Euler, CC, CCx, Dop;
        gr_ctx_init_gr_poly(Pol, ZZ);
        gr_ore_poly_ctx_init(Euler, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);
        gr_ctx_init_complex_acb(CC, prec);

        /* Euler operator dop = Σ_k x^k P_k(θ) where the indicial polynomial
           P_0(θ) = ∏ (d_i θ - n_i) has roots with repetitions and integer
           differences; the leading coefficient is constant or has a single
           root at ±1 */

        slong order = 1 + n_randint(state, 4);
        slong xdeg = 1 + n_randint(state, 3);

        gr_ore_poly_t dop, dop_mono;
        gr_poly_t P0, fac;
        gr_ore_poly_init(dop, Euler);
        gr_poly_init(P0, ZZ);
        gr_poly_init(fac, ZZ);

        status |= gr_poly_one(P0, ZZ);
        slong n = 0, d = 1, lc0 = 1;
        for (slong i = 0; i < order; i++)
        {
            ulong r = n_randint(state, 4);
            if (i == 0 || r == 0)
            {
                d = 1 + n_randint(state, 3);
                n = (slong) n_randint(state, 9) - 4;
            }
            else if (r == 1)
                n += d * (1 + n_randint(state, 3));
            /* otherwise repeat the previous root */
            status |= gr_poly_set_coeff_si(fac, 1, d, ZZ);
            status |= gr_poly_set_coeff_si(fac, 0, -n, ZZ);
            status |= gr_poly_mul(P0, P0, fac, ZZ);
            lc0 *= d;
        }

        gr_ore_poly_fit_length(dop, order + 1, Euler);
        _gr_ore_poly_set_length(dop, order + 1, Euler);
        for (slong i = 0; i <= order; i++)
        {
            gr_poly_struct * c = gr_ore_poly_coeff_ptr(dop, i, Euler);
            status |= gr_poly_set_coeff_scalar(c, 0,
                    gr_poly_coeff_srcptr(P0, i, ZZ), ZZ);
            if (i < order)
                for (slong k = 1; k <= xdeg; k++)
                    status |= gr_poly_set_coeff_si(c, k,
                            (slong) n_randint(state, 7) - 3, ZZ);
            else if (n_randint(state, 2))
                status |= gr_poly_set_coeff_si(c, 1,
                        n_randint(state, 2) ? lc0 : -lc0, ZZ);
        }

        /* exponents, shared by both computations to fix the basis */

        acb_ode_exponents_t expos;
        acb_ode_exponents_init(expos);
        status |= acb_ode_exponents(expos, dop, Euler, CC);

        /* path 0 → a → a ω → ... → a ω^(LOOP_STEPS-1) → a → 0 where |a| = 1/4
           and ω = exp(2πi/LOOP_STEPS) (with the intermediate points rounded to
           dyadic numbers) */

        gr_vec_t path;
        gr_vec_init(path, LOOP_STEPS + 3, CC);
        acb_ptr pts = path->entries;
        acb_t a, w;
        acb_init(a);
        acb_init(w);

        switch (n_randint(state, 4))
        {
            case 0: acb_set_d_d(a,  0.25,  0.  ); break;
            case 1: acb_set_d_d(a,  0.  ,  0.25); break;
            case 2: acb_set_d_d(a, -0.25,  0.  ); break;  /* on the cut */
            default: acb_set_d_d(a, 0.  , -0.25); break;
        }

        acb_zero(pts + 0);
        acb_set(pts + 1, a);
        for (slong k = 1; k < LOOP_STEPS; k++)
        {
            acb_set_si(w, 2 * k);
            acb_div_si(w, w, LOOP_STEPS, prec);
            acb_exp_pi_i(w, w, prec);
            acb_mul(pts + 1 + k, a, w, prec);
            acb_get_mid(pts + 1 + k, pts + 1 + k);
        }
        acb_set(pts + LOOP_STEPS + 1, a);
        acb_zero(pts + LOOP_STEPS + 2);

        /* operator passed to local_monodromy: same as dop, or converted to
           ball coefficients and/or to the standard derivation */

        int mono_cc = n_randint(state, 2);
        ore_algebra_t mono_alg = n_randint(state, 2)
                                 ? ORE_ALGEBRA_DERIVATIVE
                                 : ORE_ALGEBRA_EULER_DERIVATIVE;
        gr_ctx_init_gr_poly(CCx, CC);
        gr_ore_poly_ctx_init(Dop, mono_cc ? CCx : Pol, 0, mono_alg);
        gr_ore_poly_init(dop_mono, Dop);

        if (mono_cc)
        {
            gr_ore_poly_t tmp;
            gr_ctx_t CCxT;
            gr_ore_poly_ctx_init(CCxT, CCx, 0, ORE_ALGEBRA_EULER_DERIVATIVE);
            gr_ore_poly_init(tmp, CCxT);
            status |= gr_poly_set_gr_poly_other((gr_ptr) tmp, (gr_ptr) dop,
                                                Pol, CCx);
            slong sh;
            status |= gr_ore_poly_convert(dop_mono, &sh, tmp, Dop, CCxT);
            gr_ore_poly_clear(tmp, CCxT);
            gr_ctx_clear(CCxT);
        }
        else
        {
            slong sh;
            status |= gr_ore_poly_convert(dop_mono, &sh, dop, Dop, Euler);
        }

        /* test */

        acb_mat_t mono, ref;
        acb_mat_init(mono, order, order);
        acb_mat_init(ref, order, order);

        if (status != GR_SUCCESS)
            goto end;

        status |= acb_ode_local_monodromy(mono, dop_mono, Dop, expos, 0,
                                          prec);
        status |= acb_ode_transition_matrix_fixed(ref, dop, Euler, path, CC,
                                                  expos, expos, 0, prec);

        if (status != GR_SUCCESS)
            goto end;

        if (!acb_mat_overlaps(mono, ref) || !acb_mat_is_finite(mono))
        {
            flint_printf("FAIL\n\n");
            flint_printf("flint_rand_set_seed(state, 0x%wx, 0x%wx);\n",
                         seed1, seed2);
            flint_printf("dop = %{gr}\n", dop, Euler);
            flint_printf("Dop = %{gr_ctx}\n", Dop);
            flint_printf("a = %{acb}, prec = %wd\n", a, prec);
            acb_ode_exponents_println(expos);
            flint_printf("mono = %{acb_mat}\n", mono);
            flint_printf("ref = %{acb_mat}\n", ref);
            flint_abort();
        }

end:
        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        acb_mat_clear(ref);
        acb_mat_clear(mono);
        gr_ore_poly_clear(dop_mono, Dop);
        acb_clear(w);
        acb_clear(a);
        gr_vec_clear(path, CC);
        acb_ode_exponents_clear(expos);
        gr_poly_clear(fac, ZZ);
        gr_poly_clear(P0, ZZ);
        gr_ore_poly_clear(dop, Euler);
        gr_ctx_clear(Dop);
        gr_ctx_clear(CCx);
        gr_ctx_clear(CC);
        gr_ctx_clear(Euler);
        gr_ctx_clear(Pol);
    }

    gr_ctx_clear(ZZ);

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
