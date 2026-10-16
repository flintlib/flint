/*
    Copyright (C) 2019 D.H.J. Polymath
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_dirichlet.h"
#include "acb_dirichlet/impl.h"
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

/* Z(t) at prec bits; the sign, or 0 if undetermined */
static int
_hardy_z_sign(arb_t res, const arf_t t, slong prec)
{
    acb_t z;
    int msign;
    acb_init(z);
    _acb_set_arf(z, t);
    acb_dirichlet_hardy_z(z, z, NULL, NULL, 1, prec);
    msign = arb_sgn_nonzero(acb_realref(z));
    acb_get_real(res, z);
    acb_clear(z);
    return msign;
}

/*
    Refines the zero in the interval with endpoints ra, rb (Z having
    opposite signs there, and a unique zero in between, necessarily of
    odd multiplicity) to
    about prec bits by the Illinois (modified regula falsi) method.

    The iterates c are computed at wp = prec + nmag + 8 bits (cheap
    arf arithmetic, enough to represent points spaced far below the
    tolerance 2^abs_tol), but Z is evaluated at only ep = prec + 12
    bits: at the height t ~ 2^nmag, the phase of Z is determined to
    about ep - nmag bits after the binary point, and the values of Z
    near the zero (about |Z'| 2^abs_tol at the end) need only a few
    bits beyond abs_tol for the signs and the secant steps. This keeps
    the Riemann-Siegel main sum at about prec - nmag + 18 bits after
    the binary point (e.g. 84 bits at the 10^15-th zero to 114 bits,
    instead of prec + 14 = 128).

    When a secant iterate c lands so close to the zero that the sign
    of Z(c) is undetermined at ep, raising the precision (doubling it)
    would leave the fast evaluation range; instead, the points
    c -+ h with h = 2^(abs_tol - 2) are evaluated (where |Z| is about
    |Z'| h, far above the error): if the signs differ, [c - h, c + h]
    brackets the zero and is within the tolerance; otherwise (the zero
    between a and b being unique, and a sign change) Z has that sign on
    the whole of [c - h, c + h], and in particular at c. The bracket is
    clipped to [a, b], whose signs are known. Only if a sign stays
    undetermined is the precision raised (by 32 bits, then doubled).
*/
static void
_refine_hardy_z_zero_illinois_direct(arb_t res, const arf_t ra, const arf_t rb, slong prec)
{
    arf_t a, b, fa, fb, c, fc, t, h, lo, hi;
    arb_t z, z2, z3;
    slong k, nmag, abs_tol, wp, ep;
    int asign, bsign, csign, done = 0;

    arf_init(a);
    arf_init(b);
    arf_init(c);
    arf_init(fa);
    arf_init(fb);
    arf_init(fc);
    arf_init(t);
    arf_init(h);
    arf_init(lo);
    arf_init(hi);
    arb_init(z);
    arb_init(z2);
    arb_init(z3);

    arf_set(a, ra);
    arf_set(b, rb);

    nmag = arf_abs_bound_lt_2exp_si(b);
    abs_tol = nmag - prec - 4;

    wp = prec + nmag + 8;
    ep = prec + 12;
    asign = _acb_dirichlet_definite_hardy_z(z, a, &ep);
    arf_set(fa, arb_midref(z));
    bsign = _acb_dirichlet_definite_hardy_z(z, b, &ep);
    arf_set(fb, arb_midref(z));

    if (asign == bsign)
    {
        flint_throw(FLINT_ERROR, "isolate a zero before bisecting the interval\n");
    }

    arf_one(h);
    arf_mul_2exp_si(h, h, abs_tol - 2);

    for (k = 0; k < 40 && !done; k++)
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
            ep += 32;
            arf_add(c, a, b, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(c, c, -1);
        }

        csign = _hardy_z_sign(z, c, ep);

        if (csign == 0)
        {
            /* bracket the zero near c by [c - h, c + h], clipped to
               [a, b] (whose signs are known) */
            int aleft = (arf_cmp(a, b) < 0);
            arf_srcptr left = aleft ? a : b, right = aleft ? b : a;
            arf_srcptr fleft = aleft ? fa : fb, fright = aleft ? fb : fa;
            int sleft = aleft ? asign : bsign, sright = aleft ? bsign : asign;
            int lsign, hsign;

            arf_sub(lo, c, h, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_add(hi, c, h, ARF_PREC_EXACT, ARF_RND_DOWN);

            if (arf_cmp(lo, left) <= 0)
            {
                arf_set(lo, left);
                arb_set_arf(z2, fleft);
                lsign = sleft;
            }
            else
                lsign = _hardy_z_sign(z2, lo, ep);

            if (arf_cmp(hi, right) >= 0)
            {
                arf_set(hi, right);
                arb_set_arf(z3, fright);
                hsign = sright;
            }
            else
                hsign = _hardy_z_sign(z3, hi, ep);

            if (lsign != 0 && hsign != 0 && lsign != hsign)
            {
                arf_set(a, lo);
                arf_set(fa, arb_midref(z2));
                asign = lsign;
                arf_set(b, hi);
                arf_set(fb, arb_midref(z3));
                bsign = hsign;
                done = 1;
                break;
            }
            else if (lsign != 0 && lsign == hsign)
            {
                /* no zero in [lo, hi]: Z(c) has this sign; its value
                   is below the error, so take for the secant steps a
                   small value of that sign */
                csign = lsign;
                arf_set_mag(fc, arb_radref(z));
                arf_mul_2exp_si(fc, fc, -1);
                if (csign < 0)
                    arf_neg(fc, fc);
            }

            if (csign == 0)
            {
                ep += 32;
                csign = _acb_dirichlet_definite_hardy_z(z, c, &ep);
                arf_set(fc, arb_midref(z));
            }
        }
        else
        {
            arf_set(fc, arb_midref(z));
        }

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
    arf_clear(h);
    arf_clear(lo);
    arf_clear(hi);
    arb_clear(z);
    arb_clear(z2);
    arb_clear(z3);
}

/*
    The last stage of the refinement, from an interval [a, b] containing
    a unique zero (of odd multiplicity) that is already accurate to about
    half the target precision: a single secant step through a and b (Z at
    ep = prec + 12 bits) gives c within about |Z''/Z'| (b - a)^2 of the
    zero (a Newton step would need Z', whose evaluation is much more
    expensive), and the signs of Z at c -+ h, h = 2^(abs_tol - 2) (where
    |Z| is about |Z'| h, far above the error at ep) certify the bracket
    [c - h, c + h], which is within the tolerance.  This takes four
    evaluations at the target precision instead of the dozen or so of
    the Illinois iteration.  If the step or the certification fails, the
    Illinois iteration takes over on [a, b].
*/
static void
_refine_hardy_z_zero_final(arb_t res, const arf_t a, const arf_t b, slong prec)
{
    arf_t m, c, h, lo, hi;
    arb_t za, zb, z2;
    slong nmag, abs_tol, ep, wp;
    int sa, sb, sl, sh, ok = 0;

    arf_init(m);
    arf_init(c);
    arf_init(h);
    arf_init(lo);
    arf_init(hi);
    arb_init(za);
    arb_init(zb);
    arb_init(z2);

    nmag = arf_abs_bound_lt_2exp_si(b);
    abs_tol = nmag - prec - 4;
    ep = prec + 12;
    wp = prec + nmag + 8;

    sa = _hardy_z_sign(za, a, ep);
    sb = (sa != 0) ? _hardy_z_sign(zb, b, ep) : 0;
    if (sa != 0 && sb != 0 && sa != sb)
    {
        /* c = a - fa (b - a) / (fb - fa) */
        arf_sub(c, b, a, wp, ARF_RND_NEAR);
        arf_sub(m, arb_midref(zb), arb_midref(za), wp, ARF_RND_NEAR);
        arf_div(c, c, m, wp, ARF_RND_NEAR);
        arf_mul(c, c, arb_midref(za), wp, ARF_RND_NEAR);
        arf_sub(c, a, c, wp, ARF_RND_NEAR);

        arf_one(h);
        arf_mul_2exp_si(h, h, abs_tol - 2);
        arf_sub(lo, c, h, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_add(hi, c, h, ARF_PREC_EXACT, ARF_RND_DOWN);

        /* the bracket must lie in [a, b], where the zero is unique */
        if (arf_cmp(a, lo) < 0 && arf_cmp(hi, b) < 0)
        {
            sl = _hardy_z_sign(z2, lo, ep);
            sh = (sl != 0) ? _hardy_z_sign(z2, hi, ep) : 0;
            if (sl != 0 && sh != 0 && sl != sh)
            {
                arb_set_interval_arf(res, lo, hi, prec);
                ok = 1;
            }
        }
    }

    if (!ok)
        _refine_hardy_z_zero_illinois_direct(res, a, b, prec);

    arf_clear(m);
    arf_clear(c);
    arf_clear(h);
    arf_clear(lo);
    arf_clear(hi);
    arb_clear(za);
    arb_clear(zb);
    arb_clear(z2);
}

/* the target precision p1 of the first stage for the target precision
   prec at height 2^nmag: then (b - a)^2 is far below the tolerance */
#define REFINE_P1(prec, nmag) (((prec) + (nmag)) / 2 + 16)

static void _refine_hardy_z_zero_illinois(arb_t res, const arf_t ra, const arf_t rb, slong prec);

/*
    Two-stage refinement for a high target precision: the zero is first
    refined (recursively) to the precision p1 = (prec + nmag) / 2 + 16,
    where the evaluations of Z are much cheaper (in particular, they
    stay in the dfloat range of the Riemann-Siegel main sum up to higher
    target precisions), giving an interval [a, b] of width below
    2^(nmag - p1 - 4) that contains the zero; then the last stage.
*/
static void
_refine_hardy_z_zero_two_stage(arb_t res, const arf_t ra, const arf_t rb,
        slong prec, slong p1)
{
    arf_t a, b;
    arb_t z;
    slong wp = prec + arf_abs_bound_lt_2exp_si(rb) + 8;

    arf_init(a);
    arf_init(b);
    arb_init(z);

    _refine_hardy_z_zero_illinois(z, ra, rb, p1);
    arb_get_lbound_arf(a, z, wp);
    arb_get_ubound_arf(b, z, wp);
    /* (within the original bracket, whose endpoints may come in
       either order) */
    {
        arf_srcptr lo0 = (arf_cmp(ra, rb) < 0) ? ra : rb;
        arf_srcptr hi0 = (arf_cmp(ra, rb) < 0) ? rb : ra;
        if (arf_cmp(a, lo0) < 0)
            arf_set(a, lo0);
        if (arf_cmp(b, hi0) > 0)
            arf_set(b, hi0);
    }

    _refine_hardy_z_zero_final(res, a, b, prec);

    arf_clear(a);
    arf_clear(b);
    arb_clear(z);
}

/* the Illinois iteration, in two stages for a high target precision */
static void
_refine_hardy_z_zero_illinois(arb_t res, const arf_t ra, const arf_t rb, slong prec)
{
    slong nmag = arf_abs_bound_lt_2exp_si(rb);
    slong p1 = REFINE_P1(prec, nmag);

    if (prec >= p1 + 32)
        _refine_hardy_z_zero_two_stage(res, ra, rb, prec, p1);
    else
        _refine_hardy_z_zero_illinois_direct(res, ra, rb, prec);
}

static void
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

/* the zero in the ball z (rigorous, containing a unique zero, e.g. from
   the large height method) to prec bits; z may alias res */
void
_acb_dirichlet_refine_hardy_z_zero_ball(arb_t res, const arb_t z, slong prec)
{
    arf_t a, b;
    slong nmag, wp;

    if (arb_rel_accuracy_bits(z) >= prec - 2)
    {
        arb_set_round(res, z, prec);
        return;
    }

    arf_init(a);
    arf_init(b);
    nmag = arf_abs_bound_lt_2exp_si(arb_midref(z)) + 1;
    wp = prec + nmag + 8;
    arb_get_lbound_arf(a, z, wp);
    arb_get_ubound_arf(b, z, wp);

    /* accurate to about half the precision: the last stage directly
       (signs at the endpoints are determined at prec + 12 bits); else
       the full refinement (in two stages if worthwhile) */
    if (arb_rel_accuracy_bits(z) >= REFINE_P1(prec, nmag) - 8 &&
        prec >= REFINE_P1(prec, nmag) + 32)
        _refine_hardy_z_zero_final(res, a, b, prec);
    else
        _acb_dirichlet_refine_hardy_z_zero(res, a, b, prec);

    arb_set_round(res, res, prec);
    arf_clear(a);
    arf_clear(b);
}
