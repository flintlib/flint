/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Code (Claude Fable 5) */

#include "test_helpers.h"
#include "acb_types.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_ode.h"
#include "acb_poly.h"
#include "arb.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_vec.h"
#include "mag.h"

/* Evaluation of pFq(a; b; z) for z in the disk |z - 1| < 1 based on the
 * hypergeometric differential equation, compared with acb_hypgeom_pfq. */

TEST_FUNCTION_START(acb_ode_hypgeom_pfq, state)
{
    gr_ctx_t QQ, ZZi, CaCC;

    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_fmpzi(ZZi);
    gr_ctx_init_complex_ca(CaCC);

    for (slong iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        ulong seed1, seed2;
        flint_rand_get_seed(&seed1, &seed2, state);

        int status = GR_SUCCESS;

        slong prec = 64 + n_randint(state, 192);

        /* parameter ring */

        gr_ctx_t CC;
        gr_ctx_init_complex_acb(CC, prec + 64);

        gr_ctx_struct * R;
        int r_acb = 0, r_ca = 0;
        switch (n_randint(state, 8))
        {
            case 0: case 1: case 2:
                R = QQ;
                break;
            case 3: case 4:
                R = ZZi;
                break;
            case 5:
                R = CaCC;
                r_ca = 1;
                break;
            default:
                R = CC;
                r_acb = 1;
                break;
        }

        /* parameters, with degenerate configurations */

        slong q = n_randint(state, r_ca ? 3 : 6);
        slong p = n_randint(state, 3) ? q + 1 : n_randint(state, q + 2);
        slong order = q + 1;

        gr_vec_t a, b;
        gr_vec_init(a, p, R);
        gr_vec_init(b, q, R);

        for (slong i = 0; i < p + q; i++)
        {
            gr_ptr x = i < p ? gr_vec_entry_ptr(a, i, R)
                             : gr_vec_entry_ptr(b, i - p, R);
            status |= gr_randtest_small(x, state, R);
            /* exact parameters make the structure of the equation at 1
               decidable; keep genuine balls sometimes */
            if (r_acb && n_randint(state, 4))
                acb_get_mid((acb_ptr) x, (acb_ptr) x);
        }

        if (p > 0 && !n_randint(state, 4))  /* terminating series */
            status |= gr_set_si(gr_vec_entry_ptr(a, n_randint(state, p), R),
                                -(slong) n_randint(state, 4), R);

        if (p > 0 && q > 0 && !n_randint(state, 4))
        {
            /* parameters differing by integers */
            slong i = n_randint(state, p), j = n_randint(state, q);
            status |= gr_add_si(gr_vec_entry_ptr(a, i, R),
                                gr_vec_entry_srcptr(b, j, R),
                                (slong) n_randint(state, 5) - 2, R);
        }

        if (q > 0 && !n_randint(state, 4))
        {
            /* make s = Σb - Σa a (small) integer; 0 <= s <= q-1 gives
               logarithmic solutions at 1 */
            gr_ptr bq = gr_vec_entry_ptr(b, q - 1, R);
            status |= gr_set_si(bq, (slong) n_randint(state, q + 3) - 2, R);
            for (slong i = 0; i < p; i++)
                status |= gr_add(bq, bq, gr_vec_entry_srcptr(a, i, R), R);
            for (slong j = 0; j < q - 1; j++)
                status |= gr_sub(bq, bq, gr_vec_entry_srcptr(b, j, R), R);
        }

        /* pFq is undefined for nonpositive integer lower parameters */
        for (slong j = 0; j < q; j++)
        {
            gr_ptr x = gr_vec_entry_ptr(b, j, R);
            fmpz_t m;
            fmpz_init(m);
            if (gr_get_fmpz(m, x, R) == GR_SUCCESS && fmpz_sgn(m) <= 0
                    && fmpz_fits_si(m))
                status |= gr_add_si(x, x, 1 - fmpz_get_si(m), R);
            fmpz_clear(m);
        }

        /* evaluation point z = 1 + u, 0 < |u| <= 1/2 */

        acb_t u, z;
        acb_init(u);
        acb_init(z);

        switch (n_randint(state, 8))
        {
            case 0:  /* on the branch cut */
            case 1:  /* real, left of 1 */
                arb_urandom(acb_realref(u), state, prec);
                acb_mul_2exp_si(u, u, -1);
                acb_get_mid(u, u);
                if (n_randint(state, 2))
                    acb_neg(u, u);
                break;
            case 2:  /* z = pivot used internally */
                acb_set_d(u, -0.25);
                break;
            case 3:  /* ball straddling the branch cut */
                acb_urandom(u, state, prec);
                acb_mul_2exp_si(u, u, -1);
                acb_get_mid(u, u);
                arb_abs(acb_realref(u), acb_realref(u));
                arb_zero(acb_imagref(u));
                mag_one(arb_radref(acb_imagref(u)));
                mag_mul_2exp_si(arb_radref(acb_imagref(u)),
                                arb_radref(acb_imagref(u)), -20);
                break;
            case 4:  /* lower half-disk with |z| < 1: the reference series
                        converges and the value is branch-sensitive */
            {
                mag_t zabs;
                mag_init(zabs);
                do
                {
                    acb_urandom(u, state, prec);
                    acb_mul_2exp_si(u, u, -1);
                    acb_get_mid(u, u);
                    arb_abs(acb_imagref(u), acb_imagref(u));
                    arb_neg(acb_imagref(u), acb_imagref(u));
                    acb_add_ui(z, u, 1, prec);
                    acb_get_mag(zabs, z);
                }
                while (arb_is_zero(acb_imagref(u))
                       || mag_cmp_2exp_si(zabs, 0) >= 0);
                mag_clear(zabs);
                break;
            }
            default:
                acb_urandom(u, state, prec);
                acb_mul_2exp_si(u, u, -1);
                acb_get_mid(u, u);
                break;
        }
        if (acb_contains_zero(u))
            acb_set_d(u, 0.25);
        acb_add_ui(z, u, 1, prec + 32);

        /* evaluate using the differential equation, then compare */

        slong len = (order > 1 && n_randint(state, 2)) ? 2 : 1;
        acb_ptr ours = _acb_vec_init(len);
        acb_ptr ref = _acb_vec_init(len);
        acb_ptr ab_acb = _acb_vec_init(p + q);

        status |= acb_ode_hypgeom_pfq(ours, a, b, R, z, len, prec);
        status |= _gr_vec_set_other(ab_acb, a->entries, R, p, CC);
        status |= _gr_vec_set_other(ab_acb + p, b->entries, R, q, CC);

        if (status == GR_SUCCESS)
        {
            acb_hypgeom_pfq(ref, ab_acb, p, ab_acb + p, q, z, 0, prec);

            if (len > 1)
            {
                /* reference derivative (the series converges only for |z| < 1
                   or in terminating cases; otherwise, the output is not finite
                   but may have an exact imaginary part, so we replace it by an
                   indeterminate value); the series-level pfq functions need an
                   extra unit lower parameter for the k! factor */
                acb_poly_struct * pq = _acb_poly_vec_init(p + q + 2);
                acb_poly_t refjet;
                acb_poly_init(refjet);

                for (slong i = 0; i < p + q; i++)
                    acb_poly_set_acb(pq + i, ab_acb + i);
                acb_poly_one(pq + p + q);
                acb_poly_set_coeff_acb(pq + p + q + 1, 0, z);
                acb_poly_set_coeff_si(pq + p + q + 1, 1, 1);
                acb_hypgeom_pfq_series_direct(refjet, pq, p, pq + p, q + 1,
                                              pq + p + q + 1, 0, -1, len,
                                              prec);
                for (slong i = 1; i < len; i++)
                {
                    acb_poly_get_coeff_acb(ref + i, refjet, i);
                    if (!acb_is_finite(ref + i))
                        acb_indeterminate(ref + i);
                }

                acb_poly_clear(refjet);
                _acb_poly_vec_clear(pq, p + q + 2);
            }

            for (slong i = 0; i < len; i++)
            {
                if (acb_overlaps(ours + i, ref + i))
                    continue;

                flint_printf("FAIL: coefficient %wd\n\n", i);
                flint_printf("seed = 0x%wx, 0x%wx\n", seed1, seed2);
                flint_printf("R = %{gr_ctx}, p = %wd, q = %wd, prec = %wd\n",
                             R, p, q, prec);
                for (slong i2 = 0; i2 < p; i2++)
                    flint_printf("a[%wd] = %{gr}\n", i2,
                                 gr_vec_entry_srcptr(a, i2, R), R);
                for (slong j = 0; j < q; j++)
                    flint_printf("b[%wd] = %{gr}\n", j,
                                 gr_vec_entry_srcptr(b, j, R), R);
                flint_printf("z = %{acb}\n", z);
                flint_printf("ours = %{acb}\n", ours + i);
                flint_printf("ref  = %{acb}\n", ref + i);
                flint_abort();
            }
        }

        _acb_vec_clear(ab_acb, p + q);
        _acb_vec_clear(ref, len);
        _acb_vec_clear(ours, len);
        acb_clear(z);
        acb_clear(u);
        gr_vec_clear(b, R);
        gr_vec_clear(a, R);
        gr_ctx_clear(CC);
    }

    gr_ctx_clear(CaCC);
    gr_ctx_clear(ZZi);
    gr_ctx_clear(QQ);

    TEST_FUNCTION_END(state);
}
