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

/* Given lo < hi where hi - lo = N * 2^e with 0 < N < 2^30 (so that
   [lo, hi] is exactly representable as a ball), and a, b with
   lo <= a < b <= hi, enlarge [a, b] to an interval with endpoints
   on the grid lo + 2^f Z for some f <= e, so that b - a = N' * 2^f with
   0 < N' < 2^30. This maintains the invariant that the current
   interval is exactly representable as a ball, which ensures that the
   output is contained in the initial interval. */
static void
_round_to_grid(arf_t a, arf_t b, const arf_t lo, const arf_t hi)
{
    arf_t w;
    fmpz_t man, exp, q;
    slong e, f;

    arf_init(w);
    fmpz_init(man);
    fmpz_init(exp);
    fmpz_init(q);

    arf_sub(w, hi, lo, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_get_fmpz_2exp(man, exp, w);
    e = fmpz_get_si(exp);

    arf_sub(w, b, a, ARF_PREC_EXACT, ARF_RND_DOWN);
    f = arf_abs_bound_lt_2exp_si(w) - 29;
    f = FLINT_MIN(f, e);

    /* a = lo + floor((a - lo) / 2^f) 2^f */
    arf_sub(w, a, lo, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul_2exp_si(w, w, -f);
    arf_get_fmpz(q, w, ARF_RND_FLOOR);
    fmpz_set_si(exp, f);
    arf_set_fmpz_2exp(w, q, exp);
    arf_add(a, lo, w, ARF_PREC_EXACT, ARF_RND_DOWN);

    /* b = lo + ceil((b - lo) / 2^f) 2^f */
    arf_sub(w, b, lo, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul_2exp_si(w, w, -f);
    arf_get_fmpz(q, w, ARF_RND_CEIL);
    arf_set_fmpz_2exp(w, q, exp);
    arf_add(b, lo, w, ARF_PREC_EXACT, ARF_RND_DOWN);

    arf_clear(w);
    fmpz_clear(man);
    fmpz_clear(exp);
    fmpz_clear(q);
}

/* Set z to the ball [lo, hi], exactly if hi - lo has at most MAG_BITS
   significant bits (which _round_to_grid ensures). */
static void
_set_ball(arb_t z, const arf_t lo, const arf_t hi)
{
    arf_t w;
    fmpz_t man, exp;

    arf_init(w);
    fmpz_init(man);
    fmpz_init(exp);

    arf_sub(w, hi, lo, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_get_fmpz_2exp(man, exp, w);

    if (fmpz_sgn(man) > 0 && fmpz_bits(man) <= MAG_BITS && COEFF_IS_MPZ(*exp) == 0)
    {
        arf_add(arb_midref(z), lo, hi, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul_2exp_si(arb_midref(z), arb_midref(z), -1);
        mag_set_ui_2exp_si(arb_radref(z), fmpz_get_ui(man), fmpz_get_si(exp) - 1);
    }
    else
    {
        arb_set_interval_arf(z, lo, hi, ARF_PREC_EXACT);
    }

    arf_clear(w);
    fmpz_clear(man);
    fmpz_clear(exp);
}

static int
_refine_done(const arb_t x, slong prec)
{
    return arb_rel_accuracy_bits(x) >= 1.1 * prec;
}

void
arb_fmpz_poly_refine_root_arb(arb_t res, const fmpz_poly_t poly, const arb_t initial, slong prec)
{
    slong d, wp, step, wp_new;
    fmpz_poly_t deriv;
    arb_t z, m, a, b, fdz, fm, fa, t, u;
    arf_t lo, hi;
    mag_t err, err2, err3, width0;
    int attempt, progress, done;
    int sign_a, sign_b, sign_m, sign_lo;
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
    arb_init(t);
    arb_init(u);
    arf_init(lo);
    arf_init(hi);
    mag_init(err);
    mag_init(err2);
    mag_init(err3);
    mag_init(width0);

    /* We maintain an exact interval [lo, hi] containing the root in its
       interior (and no other root), together with the sign of poly at lo
       when known (sign_lo = 2 means unknown). The sign at hi is then
       -sign_lo. The ball z is an enclosure of [lo, hi]. */
    arb_get_lbound_arf(lo, initial, ARF_PREC_EXACT);
    arb_get_ubound_arf(hi, initial, ARF_PREC_EXACT);
    sign_lo = 2;
    _set_ball(z, lo, hi);
    mag_mul_2exp_si(width0, arb_radref(z), 1);

    for (step = 0; ; step++)
    {
        /* The enclosure may already be accurate enough, e.g. after
           bisections (which are not followed by an accuracy check). */
        if (step != 0 && _refine_done(z, prec))
        {
            arb_set(res, z);
            break;
        }

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

        /* exact midpoint of [lo, hi] */
        arf_add(arb_midref(m), lo, hi, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul_2exp_si(arb_midref(m), arb_midref(m), -1);
        mag_zero(arb_radref(m));

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

            /* If N(z) is contained in z, then z contains a unique root,
               and it lies in N(z). */
            if (arb_contains_interior(z, t))
            {
                if (_refine_done(t, prec))
                {
                    arb_set(res, t);
                    break;
                }

                /* Accept the refined value for the next iteration. */
                if (arb_rel_accuracy_bits(t) >= 1.1 * arb_rel_accuracy_bits(z))
                {
                    arb_get_lbound_arf(arb_midref(a), t, ARF_PREC_EXACT);
                    arb_get_ubound_arf(arb_midref(b), t, ARF_PREC_EXACT);

                    if (arf_cmp(arb_midref(a), lo) < 0)
                        arf_set(arb_midref(a), lo);
                    if (arf_cmp(arb_midref(b), hi) > 0)
                        arf_set(arb_midref(b), hi);

                    _round_to_grid(arb_midref(a), arb_midref(b), lo, hi);

                    if (!arf_equal(arb_midref(a), lo))
                        sign_lo = 2;

                    arf_set(lo, arb_midref(a));
                    arf_set(hi, arb_midref(b));

                    _set_ball(z, lo, hi);
                    continue;
                }
            }
        }

        /* Try standard Newton iteration. */
        num_newton++;
        _arb_fmpz_poly_evaluate_accurately(fm, poly, m, wp, 2 * wp);
        /* The sign of m is determined; it may be used by the bisection fallback below */
        sign_m = arf_sgn(arb_midref(fm));

        /* Lucky exact zero */
        if (sign_m == 0)
        {
            arb_set(res, m);
            break;
        }

        _arb_fmpz_poly_evaluate_accurately(fdz, deriv, m, wp, 2 * wp);
        arb_div(t, fm, fdz, wp);
        arb_get_mag(err, t);
        arb_sub(t, m, t, wp);

        /* Candidate enclosures around the Newton iterate: first a tight one
           assuming quadratic convergence, then one with the size of
           the correction as radius. The candidate is intersected with
           [lo, hi] and accepted if the signs at its endpoints differ. */
        progress = 0;
        done = 0;
        for (attempt = 0; attempt < 2 && !progress && !done; attempt++)
        {
            arb_set(u, t);

            if (attempt == 0)
            {
                /* 2^16 err^2 / max(|m|, w) + 2^-wp |m| where w is the
                   width of the initial interval */
                arb_get_mag_lower(err2, m);
                mag_max(err2, err2, width0);
                if (mag_is_zero(err2))
                    continue;
                mag_mul(err3, err, err);
                mag_div(err3, err3, err2);
                mag_mul_2exp_si(err3, err3, 16);
                arb_get_mag(err2, m);
                mag_mul_2exp_si(err2, err2, -wp);
                mag_add(err3, err3, err2);
                if (mag_cmp(err3, err) >= 0)
                    continue;
                arb_add_error_mag(u, err3);
            }
            else
            {
                arb_add_error_mag(u, err);
            }

            arb_get_lbound_arf(arb_midref(a), u, ARF_PREC_EXACT);
            arb_get_ubound_arf(arb_midref(b), u, ARF_PREC_EXACT);
            if (arf_cmp(arb_midref(a), lo) < 0)
                arf_set(arb_midref(a), lo);
            if (arf_cmp(arb_midref(b), hi) > 0)
                arf_set(arb_midref(b), hi);
            mag_zero(arb_radref(a));
            mag_zero(arb_radref(b));

            if (arf_cmp(arb_midref(a), arb_midref(b)) >= 0)
                continue;

            _round_to_grid(arb_midref(a), arb_midref(b), lo, hi);

            /* no progress */
            if (arf_equal(arb_midref(a), lo) && arf_equal(arb_midref(b), hi))
                continue;

            if (arf_equal(arb_midref(a), lo) && sign_lo != 2)
                sign_a = sign_lo;
            else
                sign_a = _arb_fmpz_poly_evaluate_sign(fa, poly, a, wp);

            /* Lucky exact zero */
            if (sign_a == 0)
            {
                arb_set(res, a);
                done = 1;
                break;
            }

            if (arf_equal(arb_midref(b), hi) && sign_lo != 2)
                sign_b = -sign_lo;
            else
                sign_b = _arb_fmpz_poly_evaluate_sign(fa, poly, b, wp);

            /* Lucky exact zero */
            if (sign_b == 0)
            {
                arb_set(res, b);
                done = 1;
                break;
            }

            /* The candidate interval brackets the root. */
            if (sign_a != sign_b)
            {
                arf_set(lo, arb_midref(a));
                arf_set(hi, arb_midref(b));
                sign_lo = sign_a;
                progress = 1;
            }
        }

        if (done)
            break;

        if (progress)
        {
            _set_ball(z, lo, hi);

            if (_refine_done(z, prec))
            {
                arb_set(res, z);
                break;
            }

            continue;
        }

        /* Fallback bisection if both Newton attempts failed. */
        num_bisect++;

        if (sign_lo == 2)
        {
            arb_set_arf(a, lo);
            sign_lo = _arb_fmpz_poly_evaluate_sign(fa, poly, a, wp);

            /* The unique root in [lo, hi] is lo. */
            if (sign_lo == 0)
            {
                arb_set(res, a);
                break;
            }
        }

        if (sign_lo == sign_m)
            arf_set(lo, arb_midref(m));   /* Sign change is on (m, hi) */
        else
            arf_set(hi, arb_midref(m));   /* Sign change is on (lo, m) */

        _set_ball(z, lo, hi);
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
    arb_clear(t);
    arb_clear(u);
    arf_clear(lo);
    arf_clear(hi);
    mag_clear(err);
    mag_clear(err2);
    mag_clear(err3);
    mag_clear(width0);
}
