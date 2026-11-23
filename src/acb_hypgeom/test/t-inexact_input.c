/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

/* the Bessel functions (0 to 5) and the incomplete gamma functions (6 to 11) */
static void
_inexact_func(acb_t r, int kind, const acb_t a, const acb_t z, slong prec)
{
    switch (kind)
    {
        case 0: acb_hypgeom_bessel_j(r, a, z, prec); break;
        case 1: acb_hypgeom_bessel_y(r, a, z, prec); break;
        case 2: acb_hypgeom_bessel_i(r, a, z, prec); break;
        case 3: acb_hypgeom_bessel_k(r, a, z, prec); break;
        case 4: acb_hypgeom_bessel_i_scaled(r, a, z, prec); break;
        case 5: acb_hypgeom_bessel_k_scaled(r, a, z, prec); break;
        case 6: case 7: case 8: acb_hypgeom_gamma_upper(r, a, z, kind - 6, prec); break;
        default: acb_hypgeom_gamma_lower(r, a, z, kind - 9, prec); break;
    }
}

/* a point of the ball: an end point (corner) or a random point */
static void
_inexact_pick(arb_t x, const arb_t b, flint_rand_t state, int corner)
{
    arf_t e, r;

    arf_init(e);
    arf_init(r);
    arf_set_mag(r, arb_radref(b));

    if (corner)
    {
        arf_set(e, r);
        if (n_randint(state, 2))
            arf_neg(e, e);
    }
    else
    {
        arf_mul_ui(e, r, n_randint(state, 2001), 80, ARF_RND_DOWN);
        arf_div_ui(e, e, 1000, 80, ARF_RND_DOWN);
        arf_sub(e, e, r, 80, ARF_RND_DOWN);
    }

    arb_set_arf(x, arb_midref(b));
    arb_add_arf(x, x, e, ARF_PREC_EXACT);
    if (!arb_contains(b, x))
        arb_get_mid_arb(x, b);

    arf_clear(e);
    arf_clear(r);
}

TEST_FUNCTION_START(acb_hypgeom_inexact_input, state)
{
    slong iter;

    /* The values at points of the input balls are contained in the
       enclosures for the balls. */
    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, z, r, p, q, v;
        mag_t m;
        slong prec;
        int kind, cplx, j;
        double sa, sz;

        acb_init(a); acb_init(z); acb_init(r); acb_init(p); acb_init(q); acb_init(v);
        mag_init(m);

        kind = n_randint(state, 12);
        cplx = (n_randint(state, 3) == 0);
        prec = 30 + n_randint(state, n_randint(state, 4) ? 60 : 300);
        sa = 0.3 * exp(log(300.0) * n_randint(state, 1000) / 1000.0);
        sz = 0.3 * exp(log(300.0) * n_randint(state, 1000) / 1000.0);

        arb_set_d(acb_realref(a), sa * (n_randint(state, 5) ? 1 : -1) * (0.05 + n_randint(state, 1000) / 1000.0));
        arb_set_d(acb_realref(z), sz * (0.05 + n_randint(state, 1000) / 1000.0));
        if (cplx)
        {
            arb_set_d(acb_imagref(a), sa * (n_randint(state, 2001) / 1000.0 - 1.0) * (n_randint(state, 2) ? 1 : 0.01));
            arb_set_d(acb_imagref(z), sz * (n_randint(state, 2001) / 1000.0 - 1.0));
        }

        if (n_randint(state, 4))
        {
            mag_set_d(m, sa * ldexp(1.0, -(int) (2 + n_randint(state, 90))));
            arb_add_error_mag(acb_realref(a), m);
            if (cplx && n_randint(state, 2))
                arb_add_error_mag(acb_imagref(a), m);
        }
        if (n_randint(state, 3))
        {
            mag_set_d(m, sz * ldexp(1.0, -(int) (2 + n_randint(state, 90))));
            arb_add_error_mag(acb_realref(z), m);
            if (cplx && n_randint(state, 2))
                arb_add_error_mag(acb_imagref(z), m);
        }

        if (!acb_contains_zero(z))
        {
            _inexact_func(r, kind, a, z, prec);

            for (j = 0; j < 3 && acb_is_finite(r); j++)
            {
                _inexact_pick(acb_realref(p), acb_realref(a), state, j < 2);
                _inexact_pick(acb_imagref(p), acb_imagref(a), state, j < 2);
                _inexact_pick(acb_realref(q), acb_realref(z), state, j < 2);
                _inexact_pick(acb_imagref(q), acb_imagref(z), state, j < 2);

                _inexact_func(v, kind, p, q, prec + 60);

                if (acb_is_finite(v) && acb_rel_accuracy_bits(v) > prec && !acb_overlaps(v, r))
                {
                    flint_printf("FAIL: containment\n\n");
                    flint_printf("kind = %d, prec = %wd\n\n", kind, prec);
                    flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                    flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
                    flint_printf("p = "); acb_printd(p, 30); flint_printf("\n\n");
                    flint_printf("q = "); acb_printd(q, 30); flint_printf("\n\n");
                    flint_printf("v = "); acb_printd(v, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        acb_clear(a); acb_clear(z); acb_clear(r); acb_clear(p); acb_clear(q); acb_clear(v);
        mag_clear(m);
    }

    /* Monotone cases: the enclosure contains the values at the two
       corners where the function is largest and smallest, and is not
       much wider than the range between them. */
    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, z, r, lo, hi, p, q;
        mag_t m, w;
        int kind;
        static const int kinds[6] = { 2, 3, 5, 7, 8, 10 };
        /* the corner of the largest value: the sign of the direction in a and z */
        static const int da[6] = { -1, 1, 1, 1, 1, -1 }, dz[6] = { 1, -1, -1, -1, -1, 1 };

        acb_init(a); acb_init(z); acb_init(r); acb_init(lo); acb_init(hi); acb_init(p); acb_init(q);
        mag_init(m); mag_init(w);

        j_again:
        kind = n_randint(state, 6);
        arb_set_d(acb_realref(a), 0.5 + 0.01 * n_randint(state, 3000));
        arb_set_d(acb_realref(z), 0.5 + 0.01 * n_randint(state, 3000));
        mag_set_d(m, ldexp(1.0, -(int) (3 + n_randint(state, 30))));
        arb_add_error_mag(acb_realref(a), m);
        mag_set_d(m, ldexp(1.0, -(int) (3 + n_randint(state, 30))));
        arb_add_error_mag(acb_realref(z), m);
        if (!arb_is_positive(acb_realref(a)) || !arb_is_positive(acb_realref(z)))
            goto j_again;

        _inexact_func(r, kinds[kind], a, z, 64);

        arb_set(acb_realref(p), acb_realref(a));
        arb_set(acb_realref(q), acb_realref(z));
        arf_set_mag(arb_midref(acb_realref(lo)), arb_radref(acb_realref(a)));
        arf_set_mag(arb_midref(acb_realref(hi)), arb_radref(acb_realref(z)));
        arb_get_mid_arb(acb_realref(p), acb_realref(p));
        arb_get_mid_arb(acb_realref(q), acb_realref(q));
        if (da[kind] > 0) arb_add(acb_realref(p), acb_realref(p), acb_realref(lo), ARF_PREC_EXACT);
        else arb_sub(acb_realref(p), acb_realref(p), acb_realref(lo), ARF_PREC_EXACT);
        if (dz[kind] > 0) arb_add(acb_realref(q), acb_realref(q), acb_realref(hi), ARF_PREC_EXACT);
        else arb_sub(acb_realref(q), acb_realref(q), acb_realref(hi), ARF_PREC_EXACT);
        /* the largest value, and the smallest one at the opposite corner */
        {
            acb_t p2, q2;
            acb_init(p2); acb_init(q2);
            arb_get_mid_arb(acb_realref(p2), acb_realref(a));
            arb_get_mid_arb(acb_realref(q2), acb_realref(z));
            if (da[kind] > 0) arb_sub(acb_realref(p2), acb_realref(p2), acb_realref(lo), ARF_PREC_EXACT);
            else arb_add(acb_realref(p2), acb_realref(p2), acb_realref(lo), ARF_PREC_EXACT);
            if (dz[kind] > 0) arb_sub(acb_realref(q2), acb_realref(q2), acb_realref(hi), ARF_PREC_EXACT);
            else arb_add(acb_realref(q2), acb_realref(q2), acb_realref(hi), ARF_PREC_EXACT);
            _inexact_func(hi, kinds[kind], p, q, 200);
            _inexact_func(lo, kinds[kind], p2, q2, 200);
            acb_clear(p2); acb_clear(q2);
        }

        /* r contains both, and its radius is at most the half-width of
           the range (and rounding) */
        acb_sub(p, hi, lo, 200);
        acb_get_mag(w, p);
        mag_mul_2exp_si(w, w, -1);
        mag_mul_ui(w, w, 3);
        mag_div_ui(w, w, 2);
        /* (an evaluation on the balls is accepted if it is within
           some bits of the attainable accuracy) */
        mag_mul_2exp_si(w, w, 10);
        /* (and the rounding at 64 bits) */
        acb_get_mag(m, hi);
        mag_mul_2exp_si(m, m, -58);
        mag_add(w, w, m);

        if (!acb_is_finite(r) || !acb_overlaps(r, lo) || !acb_overlaps(r, hi)
            || !arb_is_zero(acb_imagref(r))
            || mag_cmp(arb_radref(acb_realref(r)), w) > 0
            || arb_gt(acb_realref(lo), acb_realref(hi)))
        {
            flint_printf("FAIL: monotone\n\n");
            flint_printf("kind = %d\n\n", kinds[kind]);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
            flint_printf("lo = "); acb_printd(lo, 30); flint_printf("\n\n");
            flint_printf("hi = "); acb_printd(hi, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a); acb_clear(z); acb_clear(r); acb_clear(lo); acb_clear(hi); acb_clear(p); acb_clear(q);
        mag_clear(m); mag_clear(w);
    }

    /* K_nu(z) with an exact order on a ball z: the value at the lower
       end point (a regression test for the bound through the
       differential equation) */
    {
        acb_t nu, z, r, q, v;

        acb_init(nu); acb_init(z); acb_init(r); acb_init(q); acb_init(v);

        arb_load_str(acb_realref(nu), "104eebeff8b5df -35 0 0");
        arb_load_str(acb_realref(z), "18a9f1b7f9611b -32 2a291129 -24");
        acb_hypgeom_bessel_k(r, nu, z, 65);
        arb_get_lbound_arf(arb_midref(acb_realref(q)), acb_realref(z), ARF_PREC_EXACT);
        acb_hypgeom_bessel_k(v, nu, q, 125);

        if (!acb_overlaps(r, v))
        {
            flint_printf("FAIL: K on a ball\n\n");
            flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
            flint_printf("v = "); acb_printd(v, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(nu); acb_clear(z); acb_clear(r); acb_clear(q); acb_clear(v);
    }

    TEST_FUNCTION_END(state);
}
