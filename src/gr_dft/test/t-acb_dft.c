/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "acb.h"
#include "gr.h"
#include "gr_dft.h"

TEST_FUNCTION_START(gr_dft_acb, state)
{
    slong iter;

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, 60);
        slong prec = 53 + n_randint(state, 300);
        slong refprec = prec + 128;
        acb_ptr v, w1, w2, ref;
        slong j;
        int inverse = n_randint(state, 2);
        int status;

        if (n_randint(state, 4) == 0)
            n = UWORD(1) << (1 + n_randint(state, 8));

        v = _acb_vec_init(n);
        w1 = _acb_vec_init(n);
        w2 = _acb_vec_init(n);
        ref = _acb_vec_init(n);

        for (j = 0; j < n; j++)
        {
            acb_randtest(v + j, state, prec, 4);
            if (n_randint(state, 2))
                acb_get_mid(v + j, v + j);
        }

        /* reference: this module's ball path (plain gr_dft over an
           acb context, itself tested against naive references in
           t-dft) at elevated precision on the midpoints */
        for (j = 0; j < n; j++)
            acb_get_mid(ref + j, v + j);
        if (_gr_dft_acb(ref, ref, n, inverse, 1, refprec) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("reference failed, n = %wd\n", n);

        /* both internal paths must contain the reference and agree */
        status = _gr_dft_acb(w1, v, n, inverse, 2, prec);
        if (status != GR_SUCCESS)
            TEST_FUNCTION_FAIL("nfixed path failed, n = %wd, prec = %wd\n",
                    n, prec);
        status = _gr_dft_acb(w2, v, n, inverse, 1, prec);
        if (status != GR_SUCCESS)
            TEST_FUNCTION_FAIL("ball path failed, n = %wd, prec = %wd\n",
                    n, prec);

        for (j = 0; j < n; j++)
        {
            if (!acb_overlaps(w1 + j, ref + j))
                TEST_FUNCTION_FAIL("nfixed path containment, n = %wd, "
                        "prec = %wd, j = %wd\n", n, prec, j);
            if (!acb_overlaps(w2 + j, ref + j))
                TEST_FUNCTION_FAIL("ball path containment, n = %wd, "
                        "prec = %wd, j = %wd\n", n, prec, j);
            if (!acb_overlaps(w1 + j, w2 + j))
                TEST_FUNCTION_FAIL("paths disagree, n = %wd, j = %wd\n",
                        n, j);
        }

        /* public interface roundtrip */
        gr_dft_acb(w1, v, n, prec);
        gr_dft_acb_inverse(w2, w1, n, prec);
        for (j = 0; j < n; j++)
            if (!acb_overlaps(w2 + j, v + j))
                TEST_FUNCTION_FAIL("roundtrip, n = %wd, j = %wd\n", n, j);

        /* precomputed interface: agrees with the one-shot interface */
        {
            gr_dft_acb_pre_t Q;

            if (gr_dft_acb_precomp_init(Q, n, prec) != GR_SUCCESS)
                TEST_FUNCTION_FAIL("precomp init, n = %wd, prec = %wd\n",
                        n, prec);

            if (inverse)
                gr_dft_acb_inverse_precomp(w2, v, Q, prec);
            else
                gr_dft_acb_precomp(w2, v, Q, prec);

            for (j = 0; j < n; j++)
                if (!acb_overlaps(w2 + j, ref + j))
                    TEST_FUNCTION_FAIL("precomp containment, n = %wd, "
                            "j = %wd\n", n, j);

            gr_dft_acb_precomp_clear(Q);
        }

        _acb_vec_clear(v, n);
        _acb_vec_clear(w1, n);
        _acb_vec_clear(w2, n);
        _acb_vec_clear(ref, n);
    }


    /* Edge cases of the drop-in interface: non-finite input, all-zero
       midpoints with radii, exponents beyond the fixed-point range
       (forcing the ball fallback), and a length large enough to
       engage the threaded conversion loops. */
    {
        slong n = 64, j;
        slong prec = 64;
        acb_ptr v = _acb_vec_init(n), w = _acb_vec_init(n);
        int status;

        /* non-finite input: the fixed-point path must refuse */
        for (j = 0; j < n; j++)
            acb_one(v + j);
        acb_indeterminate(v + 17);
        status = _gr_dft_acb(w, v, n, 0, 2, prec);
        if (status == GR_SUCCESS)
            TEST_FUNCTION_FAIL("non-finite input accepted\n");
        gr_dft_acb(w, v, n, prec);   /* public interface falls back */

        /* all-zero midpoints: output midpoints are exactly zero and
           the radii come from the perturbation bound */
        for (j = 0; j < n; j++)
        {
            acb_zero(v + j);
            mag_set_ui_2exp_si(arb_radref(acb_realref(v + j)), 1, -30);
        }
        status = _gr_dft_acb(w, v, n, 1, 2, prec);
        if (status != GR_SUCCESS)
            TEST_FUNCTION_FAIL("all-zero midpoints failed\n");
        for (j = 0; j < n; j++)
            if (!arf_is_zero(arb_midref(acb_realref(w + j))) ||
                !arf_is_zero(arb_midref(acb_imagref(w + j))) ||
                mag_is_zero(arb_radref(acb_realref(w + j))))
                TEST_FUNCTION_FAIL("all-zero midpoints output, j = %wd\n", j);

        /* bignum exponent out of the fixed-point range: fixed-point
           path refuses, ball fallback still works */
        {
            fmpz_t big;
            fmpz_init(big);
            fmpz_one(big);
            fmpz_mul_2exp(big, big, 100);

            for (j = 0; j < n; j++)
                acb_set_si(v + j, 1 + (j % 5));
            acb_mul_2exp_fmpz(v + 0, v + 0, big);
            status = _gr_dft_acb(w, v, n, 0, 2, prec);
            if (status == GR_SUCCESS)
                TEST_FUNCTION_FAIL("extreme exponent accepted\n");
            status = _gr_dft_acb(w, v, n, 0, 1, prec);
            if (status != GR_SUCCESS)
                TEST_FUNCTION_FAIL("ball path on extreme exponent\n");

            /* every midpoint hugely tiny: the magnitude bound
               saturates negatively and the fixed-point path refuses */
            fmpz_neg(big, big);
            for (j = 0; j < n; j++)
            {
                acb_set_si(v + j, 1 + (j % 5));
                acb_mul_2exp_fmpz(v + j, v + j, big);
            }
            status = _gr_dft_acb(w, v, n, 0, 2, prec);
            if (status == GR_SUCCESS)
                TEST_FUNCTION_FAIL("all-tiny bignum exponents accepted\n");
            status = _gr_dft_acb(w, v, n, 0, 1, prec);
            if (status != GR_SUCCESS)
                TEST_FUNCTION_FAIL("ball path on tiny exponents\n");

            fmpz_clear(big);
        }

        /* mixed magnitudes: a coefficient with a bignum-negative
           exponent underflows to an exact zero inside the fixed-point
           conversion (with the loss covered by the error bound), and
           a large in-range slong exponent is handled through the
           global scaling; both enclosures must overlap the ball
           result */
        {
            acb_ptr w2 = _acb_vec_init(n);
            fmpz_t big;
            fmpz_init(big);
            fmpz_one(big);
            fmpz_mul_2exp(big, big, 100);
            fmpz_neg(big, big);

            for (j = 0; j < n; j++)
                acb_set_si_si(v + j, (slong) j - 7, 3 - (slong) (j % 6));
            acb_mul_2exp_fmpz(v + 1, v + 1, big);
            acb_mul_2exp_si(v + 2, v + 2, WORD_MAX / 8);
            acb_mul_2exp_si(v + 3, v + 3, -(WORD_MAX / 16));

            status = _gr_dft_acb(w, v, n, 0, 2, prec);
            if (status != GR_SUCCESS)
                TEST_FUNCTION_FAIL("mixed magnitudes refused\n");
            status = _gr_dft_acb(w2, v, n, 0, 1, prec);
            if (status != GR_SUCCESS)
                TEST_FUNCTION_FAIL("ball path on mixed magnitudes\n");
            for (j = 0; j < n; j++)
                if (!acb_overlaps(w + j, w2 + j))
                    TEST_FUNCTION_FAIL("mixed magnitudes containment, "
                            "j = %wd\n", j);

            fmpz_clear(big);
            _acb_vec_clear(w2, n);
        }

        _acb_vec_clear(v, n);
        _acb_vec_clear(w, n);
    }

    /* very high precision exceeds the fixed-point limb cap, so the
       precomputed interface falls back to ball arithmetic */
    {
        slong n = 8, j;
        slong prec = 140000;
        acb_ptr v = _acb_vec_init(n), w = _acb_vec_init(n);
        gr_dft_acb_pre_t Q;

        for (j = 0; j < n; j++)
            acb_set_si(v + j, j - 3);

        if (gr_dft_acb_precomp_init(Q, n, prec) != GR_SUCCESS)
            TEST_FUNCTION_FAIL("ball-mode precomp init\n");
        gr_dft_acb_precomp(w, v, Q, prec);
        gr_dft_acb_inverse_precomp(w, w, Q, prec);
        for (j = 0; j < n; j++)
            if (!acb_contains(w + j, v + j))
                TEST_FUNCTION_FAIL("ball-mode roundtrip, j = %wd\n", j);
        gr_dft_acb_precomp_clear(Q);

        _acb_vec_clear(v, n);
        _acb_vec_clear(w, n);
    }

    {
        slong n = 20000, j;
        slong prec = 53;
        acb_ptr v = _acb_vec_init(n), w1 = _acb_vec_init(n),
                w2 = _acb_vec_init(n);

        for (j = 0; j < n; j++)
            acb_set_si_si(v + j, (slong) (j % 97) - 48,
                    (slong) (j % 53) - 26);

        flint_set_num_threads(3);
        gr_dft_acb(w1, v, n, prec);
        gr_dft_acb_inverse(w2, w1, n, prec);
        flint_set_num_threads(1);

        for (j = 0; j < n; j++)
            if (!acb_contains(w2 + j, v + j))
                TEST_FUNCTION_FAIL("threaded roundtrip, j = %wd\n", j);

        _acb_vec_clear(v, n);
        _acb_vec_clear(w1, n);
        _acb_vec_clear(w2, n);
    }

    TEST_FUNCTION_END(state);
}
