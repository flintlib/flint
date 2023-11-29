/*
    Copyright (C) 2019 D.H.J. Polymath
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "arb_calc.h"

static void
_acb_set_arf(acb_t res, const arf_t t)
{
    acb_zero(res);
    arb_set_arf(acb_realref(res), t);
}

int
_acb_dirichlet_definite_hardy_z(arb_t res, const arf_t t, slong *pprec)
{
    int msign;
    acb_t z;
    acb_init(z);
    while (1)
    {
        _acb_set_arf(z, t);
        acb_dirichlet_hardy_z(z, z, NULL, NULL, 1, *pprec);
        msign = arb_sgn_nonzero(acb_realref(z));
        if (msign)
        {
            break;
        }
        *pprec *= 2;
    }
    acb_get_real(res, z);
    acb_clear(z);
    return msign;
}

void
_refine_hardy_z_zero_illinois(arb_t res, const arf_t ra, const arf_t rb, slong prec)
{
    arf_t a, b, fa, fb, c, fc, t;
    arb_t z;
    slong k, nmag, abs_tol, wp;
    int asign, bsign, csign;

    arf_init(a);
    arf_init(b);
    arf_init(c);
    arf_init(fa);
    arf_init(fb);
    arf_init(fc);
    arf_init(t);
    arb_init(z);

    arf_set(a, ra);
    arf_set(b, rb);

    nmag = arf_abs_bound_lt_2exp_si(b);
    abs_tol = nmag - prec - 4;

    wp = prec + nmag + 8;
    asign = _acb_dirichlet_definite_hardy_z(z, a, &wp);
    arf_set(fa, arb_midref(z));
    bsign = _acb_dirichlet_definite_hardy_z(z, b, &wp);
    arf_set(fb, arb_midref(z));

    if (asign == bsign)
    {
        flint_throw(FLINT_ERROR, "isolate a zero before bisecting the interval\n");
    }

    for (k = 0; k < 40; k++)
    {
        /* c = a - fa * (b - a) / (fb - fa) */
        arf_sub(c, b, a, wp, ARF_RND_NEAR);
        arf_sub(t, fb, fa, wp, ARF_RND_NEAR);
        arf_div(c, c, t, wp, ARF_RND_NEAR);
        arf_mul(c, c, fa, wp, ARF_RND_NEAR);
        arf_sub(c, a, c, wp, ARF_RND_NEAR);

        /* if c is not sandwiched between a and b, improve precision
           and fall back to one bisection step */
        if (!arf_is_finite(c) ||
            !((arf_cmp(a, c) < 0 && arf_cmp(c, b) < 0) ||
              (arf_cmp(b, c) < 0 && arf_cmp(c, a) < 0)))
        {
            /* flint_printf("no sandwich (k = %wd)\n", k); */
            wp += 32;
            arf_add(c, a, b, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(c, c, -1);
        }

        csign = _acb_dirichlet_definite_hardy_z(z, c, &wp);
        arf_set(fc, arb_midref(z));

        if (csign != bsign)
        {
            arf_set(a, b);
            arf_set(fa, fb);
            asign = bsign;

            arf_set(b, c);
            arf_set(fb, fc);
            bsign = csign;
        }
        else
        {
            arf_set(b, c);
            arf_set(fb, fc);
            bsign = csign;

            arf_mul_2exp_si(fa, fa, -1);
        }

        arf_sub(t, a, b, wp, ARF_RND_DOWN);
        arf_abs(t, t);

        if (arf_cmpabs_2exp_si(t, abs_tol) < 0)
            break;
    }

    /* a and b may have changed places */
    if (arf_cmp(a, b) > 0)
        arf_swap(a, b);

    arb_set_interval_arf(res, a, b, prec);

    arf_clear(a);
    arf_clear(b);
    arf_clear(c);
    arf_clear(fa);
    arf_clear(fb);
    arf_clear(fc);
    arf_clear(t);
    arb_clear(z);
}

void
_refine_hardy_z_zero_newton(arb_t res, const arf_t ra, const arf_t rb, slong prec)
{
    acb_t z, zstart;
    acb_ptr v;
    mag_t der1, der2, err;
    slong nbits, initial_prec, extraprec, wp, step;
    slong * steps;

    acb_init(z);
    acb_init(zstart);
    v = _acb_vec_init(2);
    mag_init(der1);
    mag_init(der2);
    mag_init(err);

    nbits = arf_abs_bound_lt_2exp_si(rb);
    extraprec = nbits + 10;
    initial_prec = 3 * nbits + 30;

    _refine_hardy_z_zero_illinois(acb_imagref(zstart), ra, rb, initial_prec);
    arb_set_d(acb_realref(zstart), 0.5);
    /* Real part is exactly 1/2, but need an epsilon-enclosure (for bounds)
       since we work with the complex function. */
    mag_set_ui_2exp_si(arb_radref(acb_realref(zstart)), 1, nbits - initial_prec - 4);

    /* Bound |zeta''(zstart)| for Newton error bound. */
    acb_dirichlet_zeta_deriv_bound(der1, der2, zstart);

    steps = flint_malloc(sizeof(slong) * FLINT_BITS);

    step = 0;
    steps[step] = prec;

    while (steps[step] / 2 + extraprec > initial_prec)
    {
        steps[step + 1] = steps[step] / 2 + extraprec;
        step++;
    }

    acb_set(z, zstart);

    for ( ; step >= 0; step--)
    {
        wp = steps[step] + extraprec;

        mag_set(err, arb_radref(acb_imagref(z)));
        acb_get_mid(z, z);
        acb_dirichlet_zeta_jet(v, z, 0, 2, wp);
        mag_mul(err, err, der2);
        acb_add_error_mag(v + 1, err);
        acb_div(v, v, v + 1, wp);
        acb_sub(v, z, v, wp);

        if (acb_contains(zstart, v))
        {
            acb_set(z, v);
            arb_set_d(acb_realref(z), 0.5);
        }
        else
        {
            /* can this happen? should we fallback to illinois? */
            flint_throw(FLINT_ERROR, "no inclusion for interval newton!\n");
        }
    }

    arb_set(res, acb_imagref(z));

    flint_free(steps);
    acb_clear(z);
    acb_clear(zstart);
    _acb_vec_clear(v, 2);
    mag_clear(der1);
    mag_clear(der2);
    mag_clear(err);
}

void
_acb_dirichlet_refine_hardy_z_zero(arb_t res,
        const arf_t a, const arf_t b, slong prec)
{
    slong bits;

    arb_set_interval_arf(res, a, b, prec + 8);
    bits = arb_rel_accuracy_bits(res);

    if (bits < prec)
    {
        if (prec < 4 * arf_abs_bound_lt_2exp_si(b) + 40)
            _refine_hardy_z_zero_illinois(res, a, b, prec);
        else
            _refine_hardy_z_zero_newton(res, a, b, prec);
    }

    arb_set_round(res, res, prec);
}

