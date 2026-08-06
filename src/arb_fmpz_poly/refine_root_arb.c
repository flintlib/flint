/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "fmpz_poly.h"
#include "arb_fmpz_poly.h"

#define VERBOSE 0

static void
_arb_fmpz_poly_evaluate_accurately(arb_t res, const fmpz_poly_t poly, const arb_t x, slong prec, slong initial_wp)
{
    for (slong wp = initial_wp; ; wp *= 2)
    {
        arb_fmpz_poly_evaluate_arb(res, poly, x, wp);

        if (arb_rel_accuracy_bits(res) >= prec)
            break;
    }
}

static int
_arb_fmpz_poly_evaluate_sign(arb_t res, const fmpz_poly_t poly, const arb_t x, slong initial_wp)
{
    for (slong wp = initial_wp; ; wp *= 2)
    {
        arb_fmpz_poly_evaluate_arb(res, poly, x, wp);

        if (arb_is_zero(res) || !arb_contains_zero(res))
            break;
    }

    return arf_sgn(arb_midref(res));
}

void
arb_fmpz_poly_refine_root_arb(arb_t res, const fmpz_poly_t poly, const arb_t initial, slong prec)
{
    slong d, wp, step, wp_new;
    fmpz_poly_t deriv;
    arb_t z, m, a, b, fdz, fm, fa, fb, t;
    mag_t err;
    int sign_a, sign_b, sign_m;
    slong guard = 10;

    slong num_interval, num_newton, num_bisect;

    num_interval = 0;
    num_newton = 0;
    num_bisect = 0;

    if (arb_rel_accuracy_bits(initial) >= prec)
    {
        arb_set(res, initial);
        return;
    }

    d = fmpz_poly_degree(poly);

    if (d == 1)
    {
        arb_set_fmpz(res, poly->coeffs);
        arb_div_fmpz(res, res, poly->coeffs + 1, prec + 2);
        arb_neg(res, res);
        return;
    }

    fmpz_poly_init(deriv);
    fmpz_poly_derivative(deriv, poly);
    arb_init(z);
    arb_init(m);
    arb_init(a);
    arb_init(b);
    arb_init(fdz);
    arb_init(fm);
    arb_init(fa);
    arb_init(fb);
    arb_init(t);
    mag_init(err);

    arb_set(z, initial);

    for (step = 0; ; step++)
    {
        guard += 1;

        wp_new = 2 * arb_rel_accuracy_bits(z);
        wp_new = FLINT_MIN(wp_new, 1.5 * prec);

        wp = wp_new + guard;
        wp = FLINT_MAX(wp, 64);

#if VERBOSE
        flint_printf("Step %10wd, wp = %10wd      ", step, wp);
        arb_printn(z, wp, ARB_STR_CONDENSE * 10);
        flint_printf("\n");
#endif

        arb_set_arf(m, arb_midref(z));

        num_interval++;

        /* Try an interval Newton step. This should succeed if the
           root is sufficiently well isolated. */

        arb_fmpz_poly_evaluate_arb(fdz, deriv, z, wp);
        if (!arb_contains_zero(fdz))
        {
            arb_fmpz_poly_evaluate_arb(fm, poly, m, wp);

            /* Lucky exact zero */
            if (arb_is_zero(fm))
            {
                arb_set(res, m);
                break;
            }

            arb_div(t, fm, fdz, wp);
            arb_sub(t, m, t, wp);

            if (arb_contains_interior(z, t) && arb_rel_accuracy_bits(t) >= 1.1 * prec)
            {
                arb_set(res, t);
                break;
            }

            /* Accept the refined value for the next iteration. */
            if (arb_contains_interior(z, t) && arb_rel_accuracy_bits(t) >= 1.1 * arb_rel_accuracy_bits(z))
            {
                arb_set(z, t);
                continue;
            }
        }

        /* Try standard Newton iteration, estimate error as size of correction. */
        num_newton++;
        _arb_fmpz_poly_evaluate_accurately(fm, poly, m, wp, 2 * wp);
        /* The sign of m is determined; it may be used by the bisection fallback below */
        sign_m = arf_sgn(arb_midref(fm));
        _arb_fmpz_poly_evaluate_accurately(fdz, deriv, m, wp, 2 * wp);
        arb_div(t, fm, fdz, wp);
        arb_get_mag(err, t);
        arb_sub(t, m, t, wp);
        arb_add_error_mag(t, err);

        /* The candidate enclosure maintains the isolation property. */
        if (arb_contains_interior(z, t))
        {
            arb_get_lbound_arf(arb_midref(a), t, ARF_PREC_EXACT);
            arb_get_ubound_arf(arb_midref(b), t, ARF_PREC_EXACT);
            mag_zero(arb_radref(a));
            mag_zero(arb_radref(b));

            sign_a = _arb_fmpz_poly_evaluate_sign(fa, poly, a, wp);
            sign_b = _arb_fmpz_poly_evaluate_sign(fb, poly, b, wp);

            /* Lucky exact zero */
            if (sign_a == 0)
            {
                arb_set(res, a);
                break;
            }

            /* Lucky exact zero */
            if (sign_b == 0)
            {
                arb_set(res, b);
                break;
            }

            /* The candidate enclosure brackets a root. */
            if (sign_a != sign_b)
            {
                arb_set_interval_arf(t, arb_midref(a), arb_midref(b), wp);

                /* Continue if we've made progress. */
                if (arb_contains_interior(z, t))
                {
                    arb_set(z, t);

                    if (arb_rel_accuracy_bits(t) >= 1.1 * prec)
                    {
                        arb_set(res, z);
                        break;
                    }

                    continue;
                }
            }
        }

        /* Fallback bisections if both Newton attempts failed. */
        num_bisect++;
        arb_get_lbound_arf(arb_midref(a), z, ARF_PREC_EXACT);
        arb_get_ubound_arf(arb_midref(b), z, ARF_PREC_EXACT);
        mag_zero(arb_radref(a));
        mag_zero(arb_radref(b));

        sign_a = _arb_fmpz_poly_evaluate_sign(fa, poly, a, wp);

        if (sign_a == 0)
            flint_throw(FLINT_ERROR, "arb_fmpz_poly_refine_root_arb: root encountered at endpoint\n");

        if (sign_a == sign_m)
            arb_set(a, m);   /* Sign change is on (m, b) */
        else
            arb_set(b, m);   /* Sign change is on (a, m) */

        arb_set_interval_arf(t, arb_midref(a), arb_midref(b), wp);

        /* Do further bisections until we've certainly made progress. */
        while (!(arb_contains_interior(z, t)))
        {
            num_bisect++;

            wp = arb_rel_accuracy_bits(t) + 10;
            wp = FLINT_MAX(wp, 64);

            arb_add(m, a, b, ARF_PREC_EXACT);
            arb_mul_2exp_si(m, m, -1);

            sign_m = _arb_fmpz_poly_evaluate_sign(fm, poly, m, wp);

            if (sign_a == sign_m)
                arb_set(a, m);   /* Sign change is on (m, b) */
            else
                arb_set(b, m);   /* Sign change is on (a, m) */

            arb_set_interval_arf(t, arb_midref(a), arb_midref(b), wp);
        }

        arb_set(z, t);
    }

#if VERBOSE
    flint_printf("final: "); arb_printd(res, 10); flint_printf("\n");
    flint_printf("interval %4wd  newton %4wd  bisect %4wd\n", num_interval, num_newton, num_bisect);
#endif
    (void) num_interval;
    (void) num_newton;
    (void) num_bisect;

    fmpz_poly_clear(deriv);
    arb_clear(z);
    arb_clear(m);
    arb_clear(a);
    arb_clear(b);
    arb_clear(fdz);
    arb_clear(fm);
    arb_clear(fa);
    arb_clear(fb);
    arb_clear(t);
    mag_clear(err);
}
