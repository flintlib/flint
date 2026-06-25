/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_dft.h"

/* Regression guard for the acb_dft plan root/factor tables (issue #2709).

   Every acb_dft plan builds its roots of unity by raising a base root to
   powers. When the base root was built at the call precision, the powering
   amplified its error roughly in proportion to the root index, so the stored
   roots near the top of the table carried about (index) * 2^-prec of error and
   the transform inherited it: the output radius grew about linearly in n
   (about 2490x looser than achievable at n = 2^16, prec = 128 for rad2, and
   more for the bluestein factor table used at prime lengths).

   With exact integer input (zero input radius) the output radius is purely the
   transform error, so we isolate the plan table by holding the transform
   arithmetic at prec while building the plan once at prec and once at prec+64,
   and comparing the output radii. A loose table makes rad(plan=prec) grow with
   n relative to rad(plan=prec+64); a tight table keeps the ratio a small
   constant. Comparing acb_dft at prec vs prec+64 would instead vary the
   arithmetic and is dominated by the working-precision gap, so we always vary
   only the plan precision. acb_dft_naive is intentionally not checked: it has
   its own O(n) summation widening, independent of the root table. */

#define RATIO_MAX  8.0   /* fixed: rad2 ~1.5x, bluestein ~2-4x; broken: >= 24x */
#define GROWTH_MAX 4.0   /* fixed: ~1x; broken: ~54x over 2^8..2^14           */

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

static double
max_rad(acb_srcptr v, slong n)
{
    arb_t m, r;
    double d;
    slong i;

    arb_init(m);
    arb_init(r);
    for (i = 0; i < n; i++)
    {
        arb_get_rad_arb(r, acb_realref(v + i));
        if (arb_gt(r, m))
            arb_set(m, r);
        arb_get_rad_arb(r, acb_imagref(v + i));
        if (arb_gt(r, m))
            arb_set(m, r);
    }
    d = arf_get_d(arb_midref(m), ARF_RND_NEAR);
    arb_clear(m);
    arb_clear(r);
    return d;
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

static void
check_ratio(double ratio, double bound, const char * what)
{
    if (ratio > bound)
    {
        flint_printf("FAIL (looseness): %s : ratio %g > %g (see issue #2709)\n",
                     what, ratio, bound);
        flint_abort();
    }
}

/* rad2 plan built at prec vs prec+64, both transforms run at prec. */
static double
rad2_plan_ratio(int e, slong prec)
{
    slong n = WORD(1) << e;
    acb_ptr in = _acb_vec_init(n);
    acb_ptr lo = _acb_vec_init(n);
    acb_ptr hi = _acb_vec_init(n);
    acb_dft_rad2_t plo, phi;
    double rlo, rhi;

    fill_pattern(in, n);
    acb_dft_rad2_init(plo, e, prec);
    acb_dft_rad2_init(phi, e, prec + 64);
    acb_dft_rad2_precomp(lo, in, plo, prec);
    acb_dft_rad2_precomp(hi, in, phi, prec);
    acb_dft_rad2_clear(plo);
    acb_dft_rad2_clear(phi);

    check_overlap(lo, hi, n, "rad2 plan prec vs prec+64");
    rlo = max_rad(lo, n);
    rhi = max_rad(hi, n);

    _acb_vec_clear(in, n);
    _acb_vec_clear(lo, n);
    _acb_vec_clear(hi, n);
    return (rhi > 0.0) ? rlo / rhi : 1.0;
}

/* bluestein plan (any length >= 30) built at prec vs prec+64, run at prec. */
static double
bluestein_plan_ratio(slong n, slong prec)
{
    acb_ptr in = _acb_vec_init(n);
    acb_ptr lo = _acb_vec_init(n);
    acb_ptr hi = _acb_vec_init(n);
    acb_dft_bluestein_t plo, phi;
    double rlo, rhi;

    fill_pattern(in, n);
    acb_dft_bluestein_init(plo, n, prec);
    acb_dft_bluestein_init(phi, n, prec + 64);
    acb_dft_bluestein_precomp(lo, in, plo, prec);
    acb_dft_bluestein_precomp(hi, in, phi, prec);
    acb_dft_bluestein_clear(plo);
    acb_dft_bluestein_clear(phi);

    check_overlap(lo, hi, n, "bluestein plan prec vs prec+64");
    rlo = max_rad(lo, n);
    rhi = max_rad(hi, n);

    _acb_vec_clear(in, n);
    _acb_vec_clear(lo, n);
    _acb_vec_clear(hi, n);
    return (rhi > 0.0) ? rlo / rhi : 1.0;
}

/* the public acb_dft dispatcher at prec vs a tight reference plan (same kind,
   built at prec+64, run at prec). kind 0 = rad2 (power of two), 1 = bluestein. */
static double
dispatcher_ratio(slong n, int e, slong prec, int kind)
{
    acb_ptr in = _acb_vec_init(n);
    acb_ptr wd = _acb_vec_init(n);
    acb_ptr ref = _acb_vec_init(n);
    double rd, rr;

    fill_pattern(in, n);
    acb_dft(wd, in, n, prec);

    if (kind == 0)
    {
        acb_dft_rad2_t p;
        acb_dft_rad2_init(p, e, prec + 64);
        acb_dft_rad2_precomp(ref, in, p, prec);
        acb_dft_rad2_clear(p);
    }
    else
    {
        acb_dft_bluestein_t p;
        acb_dft_bluestein_init(p, n, prec + 64);
        acb_dft_bluestein_precomp(ref, in, p, prec);
        acb_dft_bluestein_clear(p);
    }

    check_overlap(wd, ref, n, "acb_dft dispatcher vs tight reference");
    rd = max_rad(wd, n);
    rr = max_rad(ref, n);

    _acb_vec_clear(in, n);
    _acb_vec_clear(wd, n);
    _acb_vec_clear(ref, n);
    return (rr > 0.0) ? rd / rr : 1.0;
}

TEST_FUNCTION_START(acb_dft_accuracy, state)
{
    slong prec = 128;
    int es[5] = { 8, 11, 12, 14, 16 };
    int ne = 4;                 /* default tops out at 2^14 (~0.1s)           */
    int ehi, e_disp, i;
    slong nb;
    double r_lo, r_hi = 0.0;

    if (flint_test_multiplier() > 1)
        ne = 5;                 /* exercise the headline n = 2^16 too         */
    ehi = es[ne - 1];

    /* (A) rad2 plan tightness, isolated, plus (B) the ratio must not grow
       with n (a ratio of ratios, so per-platform rounding cancels). */
    r_lo = rad2_plan_ratio(es[0], prec);
    for (i = 0; i < ne; i++)
    {
        double ratio = rad2_plan_ratio(es[i], prec);

        check_ratio(ratio, RATIO_MAX, "rad2 plan");
        if (es[i] == ehi)
            r_hi = ratio;
    }
    if (r_hi > GROWTH_MAX * r_lo)
    {
        flint_printf("FAIL (growth): rad2 ratio(2^%d)=%g vs ratio(2^%d)=%g "
                     "(factor %g > %g) (see issue #2709)\n",
                     ehi, r_hi, es[0], r_lo,
                     r_lo > 0 ? r_hi / r_lo : 0.0, GROWTH_MAX);
        flint_abort();
    }

    /* tiny n and small precision exercise the worst guard margin */
    for (i = 1; i <= 3; i++)
        check_ratio(rad2_plan_ratio(i, prec), RATIO_MAX, "rad2 plan tiny n");
    check_ratio(rad2_plan_ratio(10, 8), RATIO_MAX, "rad2 plan small prec");

    /* (.5) bluestein factor table: a prime >= 100 and a prime in [30, 100) */
    nb = (flint_test_multiplier() > 1) ? 16411 : 1009;
    check_ratio(bluestein_plan_ratio(nb, prec), RATIO_MAX, "bluestein large prime");
    check_ratio(bluestein_plan_ratio(97, prec), RATIO_MAX, "bluestein n=97");
    check_ratio(bluestein_plan_ratio(29, prec), RATIO_MAX, "bluestein n=29");

    /* prec > 128 forces heap-allocated mantissas: this case additionally lets
       `make valgrind` catch a regression of the np-length scratch clear in
       acb_dft_bluestein_precomp (which only leaks above the inline limb cap). */
    check_ratio(bluestein_plan_ratio(1009, 256), RATIO_MAX, "bluestein prec 256");

    /* the public acb_dft dispatcher is tight on both routed paths */
    e_disp = (flint_test_multiplier() > 1) ? 16 : 14;
    check_ratio(dispatcher_ratio(WORD(1) << e_disp, e_disp, prec, 0), RATIO_MAX,
                "acb_dft power of two");
    check_ratio(dispatcher_ratio(nb, 0, prec, 1), RATIO_MAX, "acb_dft prime");

    TEST_FUNCTION_END(state);
}
