/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Accuracy of the complex DFT: for various lengths, precisions and
   input distributions, run both internal paths of gr_dft_acb (acb
   ball arithmetic and fixed-point arithmetic) on the same input and
   report the largest output radius together with the precision loss.

   The loss is measured in bits at the scale of the output: an ideal
   enclosure of an output vector of magnitude M at precision p would
   have radii around M 2^-p, so

       loss = prec - (log2 max_k |w_k| - log2 max_k rad(w_k)).

   Components that are tiny due to cancellation carry fewer relative
   bits than this scale-level measure indicates, in both modes alike.

   Input distributions:

       uniform    acb_urandom: unit-scale midpoints (a benign, flat
                  magnitude profile)
       randtest   acb_randtest with mag_bits = 12: non-uniform
                  magnitudes with exponents spanning roughly
                  +-2^12, including zero and pure-real/imaginary
                  special values (the fixed-point path handles the
                  wild exponents through its global scaling)

   each in an exact variant (zero radii) and an inexact variant where
   every component is perturbed by a relative error of 2^(4 - prec),
   so that the propagated input radii and the arithmetic error are of
   comparable size. */

#include <stdio.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "arf.h"
#include "acb.h"
#include "gr.h"
#include "gr_dft.h"

static void
max_rad_and_mid(mag_t maxrad, mag_t maxmid, acb_srcptr w, slong n)
{
    mag_t t;
    slong j;

    mag_init(t);
    mag_zero(maxrad);
    mag_zero(maxmid);

    for (j = 0; j < n; j++)
    {
        mag_max(maxrad, maxrad, arb_radref(acb_realref(w + j)));
        mag_max(maxrad, maxrad, arb_radref(acb_imagref(w + j)));

        arf_get_mag(t, arb_midref(acb_realref(w + j)));
        mag_max(maxmid, maxmid, t);
        arf_get_mag(t, arb_midref(acb_imagref(w + j)));
        mag_max(maxmid, maxmid, t);
    }

    mag_clear(t);
}

static void
report(const char * mode, acb_srcptr w, slong n, slong prec, int status)
{
    mag_t maxrad, maxmid;

    if (status != GR_SUCCESS)
    {
        flint_printf("  %-7s unable\n", mode);
        return;
    }

    mag_init(maxrad);
    mag_init(maxmid);
    max_rad_and_mid(maxrad, maxmid, w, n);

    if (mag_is_zero(maxrad))
    {
        flint_printf("  %-7s max rad = 0 (exact)\n", mode);
    }
    else
    {
        double lrad = mag_get_d_log2_approx(maxrad);
        double lmid = mag_get_d_log2_approx(maxmid);
        double loss = (double) prec - (lmid - lrad);

        flint_printf("  %-7s max rad = 2^%.1f   loss = %.1f bits\n",
                mode, lrad, loss);
    }

    mag_clear(maxrad);
    mag_clear(maxmid);
}

int main(void)
{
    flint_rand_t state;
    ulong n_tab[] = { 64, 256, 1024, 4096, 1000, 1009 };
    slong prec_tab[] = { 53, 128, 512, 2048 };
    slong ni, pi, style;

    flint_rand_init(state);

    for (ni = 0; ni < 6; ni++)
    for (pi = 0; pi < 4; pi++)
    {
        ulong n = n_tab[ni];
        slong prec = prec_tab[pi];

        flint_printf("n = %wu, prec = %wd\n", n, prec);

        for (style = 0; style < 4; style++)
        {
            int uniform = (style < 2);
            int inexact = (style & 1);
            acb_ptr v = _acb_vec_init(n);
            acb_ptr w1 = _acb_vec_init(n);
            acb_ptr w2 = _acb_vec_init(n);
            slong j;
            int status1, status2;

            for (j = 0; j < (slong) n; j++)
            {
                if (uniform)
                {
                    acb_urandom(v + j, state, prec);
                }
                else
                {
                    acb_randtest(v + j, state, prec, 12);
                    acb_get_mid(v + j, v + j);
                }

                if (inexact)
                {
                    /* relative perturbation 2^(4 - prec) per part */
                    mag_t e;
                    mag_init(e);
                    arf_get_mag(e, arb_midref(acb_realref(v + j)));
                    mag_mul_2exp_si(e, e, 4 - prec);
                    mag_add(arb_radref(acb_realref(v + j)),
                            arb_radref(acb_realref(v + j)), e);
                    arf_get_mag(e, arb_midref(acb_imagref(v + j)));
                    mag_mul_2exp_si(e, e, 4 - prec);
                    mag_add(arb_radref(acb_imagref(v + j)),
                            arb_radref(acb_imagref(v + j)), e);
                    mag_clear(e);
                }
            }

            flint_printf(" %s, %s:\n", uniform ? "uniform" : "randtest",
                    inexact ? "inexact" : "exact");

            status1 = _gr_dft_acb(w1, v, n, 0, 1, prec);
            report("acb", w1, n, prec, status1);

            status2 = _gr_dft_acb(w2, v, n, 0, 2, prec);
            report("nfixed", w2, n, prec, status2);

            _acb_vec_clear(v, n);
            _acb_vec_clear(w1, n);
            _acb_vec_clear(w2, n);
        }

        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
