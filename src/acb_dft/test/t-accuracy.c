/*
    Copyright (C) 2026 Edgar Costa
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_dft.h"

/* Accuracy guard for the acb_dft plans (historically issue #2709,
   where loosely built root tables made the output radius grow about
   linearly in n). The module is now a wrapper around gr_dft, whose
   fixed-point path guarantees an absolute output error below
   2^-(prec+2) relative to the input scale; with exact integer input
   (zero input radius) the output radius is purely the transform
   error, so we assert directly that the radius stays below
   2^-(prec-MARGIN) times the largest output at every length. This
   guards the entire limb-selection and error-bound machinery: any
   regression that lets the error grow with n shows up as a shrinking
   margin. A plan built at prec+64 must also stay consistent (overlap)
   with the plan built at prec. */

#define MARGIN 16

static void
fill_pattern(acb_ptr v, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
    {
        /* exact integers: zero input radius */
        arb_set_si(acb_realref(v + i), (i % 19) - 9);
        arb_set_si(acb_imagref(v + i), (i % 13) - 6);
    }
}

static void
check_overlap(acb_srcptr a, acb_srcptr b, slong n, const char * what)
{
    slong i;

    for (i = 0; i < n; i++)
    {
        if (!acb_overlaps(a + i, b + i))
        {
            flint_printf("FAIL (overlap): %s, index %wd\n", what, i);
            flint_abort();
        }
    }
}

/* max radius <= 2^-(prec-MARGIN) max |midpoint|, in mag arithmetic
   (robust at high precision, where the radii underflow doubles) */
static void
check_radius(acb_srcptr w, slong n, slong prec, const char * what)
{
    mag_t rad, mid, t;
    slong i;

    mag_init(rad);
    mag_init(mid);
    mag_init(t);

    for (i = 0; i < n; i++)
    {
        mag_max(rad, rad, arb_radref(acb_realref(w + i)));
        mag_max(rad, rad, arb_radref(acb_imagref(w + i)));
        arf_get_mag(t, arb_midref(acb_realref(w + i)));
        mag_max(mid, mid, t);
        arf_get_mag(t, arb_midref(acb_imagref(w + i)));
        mag_max(mid, mid, t);
    }

    mag_mul_2exp_si(mid, mid, -(prec - MARGIN));
    if (mag_cmp(rad, mid) > 0)
    {
        flint_printf("FAIL (radius): %s, n = %wd, prec = %wd\n",
                what, n, prec);
        flint_printf("max rad = "); mag_printd(rad, 5);
        flint_printf(", bound = "); mag_printd(mid, 5);
        flint_printf("\n");
        flint_abort();
    }

    mag_clear(rad);
    mag_clear(mid);
    mag_clear(t);
}

/* one length: the one-shot and the plans at prec and prec + 64 all
   agree, and every output radius meets the precision bound */
static void
check_length(slong n, slong prec)
{
    acb_ptr in = _acb_vec_init(n);
    acb_ptr w = _acb_vec_init(n);
    acb_ptr wp = _acb_vec_init(n);
    acb_dft_pre_t plo, phi;

    fill_pattern(in, n);

    acb_dft(w, in, n, prec);
    check_radius(w, n, prec, "acb_dft");

    acb_dft_precomp_init(plo, n, prec);
    acb_dft_precomp_init(phi, n, prec + 64);
    acb_dft_precomp(wp, in, plo, prec);
    check_overlap(w, wp, n, "one-shot vs plan");
    check_radius(wp, n, prec, "plan at prec");
    acb_dft_precomp(wp, in, phi, prec);
    check_overlap(w, wp, n, "plan prec vs prec+64");
    acb_dft_precomp_clear(plo);
    acb_dft_precomp_clear(phi);

    _acb_vec_clear(in, n);
    _acb_vec_clear(w, n);
    _acb_vec_clear(wp, n);
}

TEST_FUNCTION_START(acb_dft_accuracy, state)
{
    slong prec = 128;
    int es[5] = { 8, 11, 12, 14, 16 };
    int ne = 4;                 /* default tops out at 2^14              */
    int i;
    slong nb;

    if (flint_test_multiplier() > 1)
        ne = 5;                 /* exercise n = 2^16 too                 */

    /* powers of two across the size range: the radius bound holding
       uniformly in n is the growth guard */
    for (i = 0; i < ne; i++)
        check_length(WORD(1) << es[i], prec);

    /* tiny n */
    for (i = 1; i <= 3; i++)
        check_length(WORD(1) << i, prec);

    /* small and large working precision */
    check_length(1024, 8);
    check_length(1009, 256);

    /* prime lengths */
    nb = (flint_test_multiplier() > 1) ? 16411 : 1009;
    check_length(nb, prec);
    check_length(97, prec);
    check_length(29, prec);

    TEST_FUNCTION_END(state);
}
