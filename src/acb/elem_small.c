/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb/impl.h"
#include "arb/elem.h"
#include "acf/elem.h"

/*
    The complex elementary functions at one and two limbs by the register
    paths shared with acf and nfloat (acf/elem.h): the midpoint through
    the frames, each part rounded with its error as the radius
    (_arb_small_set_1 / _2), plus the propagation of the input radii
    through a bound on the derivative over the ball, for a radius below
    2^-8 (wider balls take the general compositions, whose own handling
    of wide intervals is better):

        exp: |exp(z + d) - exp(z)| <= |exp z| (e^|d| - 1) <= 2 |d| |exp z|
        exp(pi i z): pi |d| |exp(pi i z)| e^(pi |d|) <= 2^(2 + 1) |d| |E|
        sin, cos: |d| cosh(|b| + |d|) <= 2 |d| cosh b
        sin(pi z), cos(pi z): 2 pi |d| cosh(pi b) <= 2^3 |d| cosh(pi b)
        log: |d| / (|z| - |d|) with |z| at least the larger part, for
             |d| < 2^(E - 3), E the larger exponent

    with |d| <= the sum of the radii, and the magnitudes from the
    exponents of the frames (|exp z| < 2^(k + 1) for the frame of E with
    k = re->k; cosh b < 2^(k + 2) for the frame of the real part of the
    cosine, k = cre->k).
*/

#if FLINT_BITS == 64

static void
_frame_set(arb_t res, const _acf_frame_struct * f, int two, slong prec)
{
    if (two)
        _arb_small_set_2(res, f->y[3], f->y[2], f->y[1], f->y[0], f->k, f->err, f->sgnbit, prec);
    else
        _arb_small_set_1(res, f->y[2], f->y[1], f->y[0], f->k, f->err, f->sgnbit, prec);
}

static int
_frame_nonzero(const _acf_frame_struct * f, int two)
{
    return two ? ((f->y[3] | f->y[2]) != 0) : ((f->y[2] | f->y[1]) != 0);
}

int
_acb_elem_small(acb_t r1, acb_t r2, const acb_t z, int f, slong prec)
{
    arf_srcptr a = arb_midref(acb_realref(z)), b = arb_midref(acb_imagref(z));
    _acf_frame_struct F1, F2, F3, F4;
    nn_srcptr ad, bd;
    slong an, bn, ea, eb;
    ulong da1, da0, db1, db0;
    int two, asgn, bsgn, st, pi;
    mag_t r, m;
    arb_t t1, t2, t3, t4;

    if (prec > 2 * FLINT_BITS || arf_is_special(a) || arf_is_special(b)
        || ARF_SIZE(a) > 2 || ARF_SIZE(b) > 2
        || COEFF_IS_MPZ(ARF_EXP(a)) || COEFF_IS_MPZ(ARF_EXP(b)))
        return 0;

    ea = ARF_EXP(a);
    eb = ARF_EXP(b);
    if (FLINT_ABS(ea) >= (WORD(1) << 40) || FLINT_ABS(eb) >= (WORD(1) << 40))
        return 0;

    /* the radius: zero, or small; the logarithm's ball must stay off
       the branch cut (an imaginary part containing zero with a negative
       real part) */
    mag_init(r);
    mag_init(m);
    mag_add(r, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
    if (!mag_is_finite(r) || (!mag_is_zero(r) && mag_cmp_2exp_si(r, -8) >= 0)
        || (f == ACB_ELEM_LOG && ARF_SGNBIT(a) && arf_cmpabs_mag(b, arb_radref(acb_imagref(z))) <= 0))
    {
        mag_clear(r);
        mag_clear(m);
        return 0;
    }

    ARF_GET_MPN_READONLY(ad, an, a);
    ARF_GET_MPN_READONLY(bd, bn, b);
    asgn = ARF_SGNBIT(a);
    bsgn = ARF_SGNBIT(b);
    two = !(prec <= FLINT_BITS && an == 1 && bn == 1);
    da1 = ad[an - 1];
    da0 = (an == 2) ? ad[0] : 0;
    db1 = bd[bn - 1];
    db0 = (bn == 2) ? bd[0] : 0;
    pi = (f == ACB_ELEM_EXP_PI_I || f == ACB_ELEM_SIN_COS_PI);

    switch (f)
    {
        case ACB_ELEM_EXP:
        case ACB_ELEM_EXP_PI_I:
            st = two ? _acf_small_exp_2(&F1, &F2, da1, da0, ea, asgn, db1, db0, eb, bsgn, pi)
                     : _acf_small_exp_1(&F1, &F2, da1, ea, asgn, db1, eb, bsgn, pi);
            break;
        case ACB_ELEM_SIN_COS:
        case ACB_ELEM_SIN_COS_PI:
            st = two ? _acf_small_sin_cos_2(&F1, &F2, &F3, &F4, da1, da0, ea, asgn, db1, db0, eb, bsgn, pi)
                     : _acf_small_sin_cos_1(&F1, &F2, &F3, &F4, da1, ea, asgn, db1, eb, bsgn, pi);
            break;
        default:
            st = two ? _acf_small_log_2(&F1, &F2, da1, da0, ea, asgn, db1, db0, eb, bsgn)
                     : _acf_small_log_1(&F1, &F2, da1, ea, asgn, db1, eb, bsgn);
            break;
    }

    if (st != 0 || !_frame_nonzero(&F1, two) || !_frame_nonzero(&F2, two) ||
        ((f == ACB_ELEM_SIN_COS || f == ACB_ELEM_SIN_COS_PI) && (!_frame_nonzero(&F3, two) || !_frame_nonzero(&F4, two))))
    {
        mag_clear(r);
        mag_clear(m);
        return 0;
    }

    /* (into temporaries: the outputs may alias z) */
    arb_init(t1);
    arb_init(t2);
    arb_init(t3);
    arb_init(t4);
    _frame_set(t1, &F1, two, prec);
    _frame_set(t2, &F2, two, prec);
    if (f == ACB_ELEM_SIN_COS || f == ACB_ELEM_SIN_COS_PI)
    {
        _frame_set(t3, &F3, two, prec);
        _frame_set(t4, &F4, two, prec);
    }

    /* the propagated radii, from the gradients of the parts (the
       magnitudes of the parts themselves, or of the other pair) */
    if (!mag_is_zero(r))
    {
        mag_t ra, rb, m1, m2, m3, m4, u;
        slong sh;

        mag_init(ra);
        mag_init(rb);
        mag_init(m1);
        mag_init(m2);
        mag_init(m3);
        mag_init(m4);
        mag_init(u);
        mag_set(ra, arb_radref(acb_realref(z)));
        mag_set(rb, arb_radref(acb_imagref(z)));

        if (f == ACB_ELEM_LOG)
        {
            /* Re = log |z|, Im = arg z: the gradients (a, b) / |z|^2 and
               (-b, a) / |z|^2, |z| - |d| at least the larger part less
               |d| (which is below |z| / 4) */
            mag_t lz;
            mag_init(lz);
            if (mag_cmp_2exp_si(r, FLINT_MAX(ea, eb) - 3) >= 0)
            {
                mag_clear(lz);
                mag_clear(ra); mag_clear(rb); mag_clear(m1); mag_clear(m2); mag_clear(m3); mag_clear(m4); mag_clear(u);
                arb_clear(t1); arb_clear(t2); arb_clear(t3); arb_clear(t4);
                mag_clear(r);
                mag_clear(m);
                return 0;
            }
            /* (|z| - |d|)^2 from below */
            arf_get_mag_lower(lz, a);
            mag_mul_lower(lz, lz, lz);
            arf_get_mag_lower(u, b);
            mag_mul_lower(u, u, u);
            mag_add_lower(lz, lz, u);
            mag_sqrt_lower(lz, lz);
            mag_sub_lower(lz, lz, r);
            mag_mul_lower(lz, lz, lz);
            arf_get_mag(m1, a);
            mag_add(m1, m1, ra);            /* |a| + ra */
            arf_get_mag(m2, b);
            mag_add(m2, m2, rb);            /* |b| + rb */
            mag_mul(m3, m1, ra);
            mag_addmul(m3, m2, rb);         /* ra |a| + rb |b| */
            mag_div(m3, m3, lz);
            mag_mul(m4, m2, ra);
            mag_addmul(m4, m1, rb);         /* ra |b| + rb |a| */
            mag_div(m4, m4, lz);
            mag_add(arb_radref(t1), arb_radref(t1), m3);
            mag_add(arb_radref(t2), arb_radref(t2), m4);
            mag_clear(lz);
        }
        else
        {
            /* With P + Qi the first result and P' + Q'i the second (sin
               and cos; for exp both are exp z itself),

                   sin: dP = (P', -Q'), dQ = (Q', P'), dP' = (-P, Q), dQ' = (-Q, -P)
                   exp: dP = (P, -Q), dQ = (Q, P)
                   exp(pi i z): dP = pi (-Q, -P), dQ = pi (P, -Q)

               (the partial derivatives in a and b), so that the radius of
               each part is the radii times the magnitudes of two parts
               (m3, m4 for the first pair, m1, m2 for the second), times
               e^|d| (e^(pi |d|) and pi for the functions of pi z: 2^3 in
               all), plus a second-order term |d|^2 2^sh, 2^sh bounding
               the magnitudes (|exp z| < 2^(k + 1), cosh b < 2^(k + 2)). */
            int pi_f = (f == ACB_ELEM_EXP_PI_I || f == ACB_ELEM_SIN_COS_PI);
            /* (the magnitudes of the parts: the rounded midpoints plus
               their radii) */
            arb_get_mag(m1, t1);
            arb_get_mag(m2, t2);
            if (f == ACB_ELEM_SIN_COS || f == ACB_ELEM_SIN_COS_PI)
            {
                arb_get_mag(m3, t3);
                arb_get_mag(m4, t4);
                sh = F3.k + 3;
            }
            else
            {
                if (f == ACB_ELEM_EXP)
                {
                    mag_set(m3, m1);
                    mag_set(m4, m2);
                }
                else
                {
                    mag_set(m3, m2);
                    mag_set(m4, m1);
                }
                sh = F1.k + 2;
            }
            mag_mul(u, r, r);
            mag_mul_2exp_si(u, u, sh + (pi_f ? 4 : 0));

            /* a part below r 2^sh (near the real axis, say) would lose its
               relative accuracy to the second-order term u: the
               compositions keep it (sharper by the products of the
               parts) */
            mag_mul_2exp_si(m, r, sh + (pi_f ? 4 : 0));
            if (mag_cmp(m1, m) < 0 || mag_cmp(m2, m) < 0 ||
                ((f == ACB_ELEM_SIN_COS || f == ACB_ELEM_SIN_COS_PI) && (mag_cmp(m3, m) < 0 || mag_cmp(m4, m) < 0)))
            {
                mag_clear(ra); mag_clear(rb); mag_clear(m1); mag_clear(m2); mag_clear(m3); mag_clear(m4); mag_clear(u);
                arb_clear(t1); arb_clear(t2); arb_clear(t3); arb_clear(t4);
                mag_clear(r);
                mag_clear(m);
                return 0;
            }

            /* the factor e^|d| < 1 + 2^-7, resp. pi e^(pi |d|) < 3.2, on
               the first-order terms: first pair ra m3 + rb m4 and ra m4
               + rb m3; second pair ra m1 + rb m2 and ra m2 + rb m1 */
            if (pi_f)
            {
                mag_set_d(m, 3.2);
                mag_mul(ra, ra, m);
                mag_mul(rb, rb, m);
            }
            else
            {
                mag_mul_ui(ra, ra, 257);
                mag_mul_2exp_si(ra, ra, -8);
                mag_mul_ui(rb, rb, 257);
                mag_mul_2exp_si(rb, rb, -8);
            }
            mag_mul(m, m3, ra);
            mag_addmul(m, m4, rb);
            mag_add(m, m, u);
            mag_add(arb_radref(t1), arb_radref(t1), m);
            mag_mul(m, m4, ra);
            mag_addmul(m, m3, rb);
            mag_add(m, m, u);
            mag_add(arb_radref(t2), arb_radref(t2), m);
            if (f == ACB_ELEM_SIN_COS || f == ACB_ELEM_SIN_COS_PI)
            {
                mag_mul(m, m1, ra);
                mag_addmul(m, m2, rb);
                mag_add(m, m, u);
                mag_add(arb_radref(t3), arb_radref(t3), m);
                mag_mul(m, m2, ra);
                mag_addmul(m, m1, rb);
                mag_add(m, m, u);
                mag_add(arb_radref(t4), arb_radref(t4), m);
            }
        }

        mag_clear(ra); mag_clear(rb); mag_clear(m1); mag_clear(m2); mag_clear(m3); mag_clear(m4); mag_clear(u);
    }

    if (f == ACB_ELEM_SIN_COS || f == ACB_ELEM_SIN_COS_PI)
    {
        if (r1 != NULL)
        {
            arb_swap(acb_realref(r1), t1);
            arb_swap(acb_imagref(r1), t2);
        }
        if (r2 != NULL)
        {
            arb_swap(acb_realref(r2), t3);
            arb_swap(acb_imagref(r2), t4);
        }
    }
    else
    {
        arb_swap(acb_realref(r1), t1);
        arb_swap(acb_imagref(r1), t2);
    }
    arb_clear(t3);
    arb_clear(t4);

    arb_clear(t1);
    arb_clear(t2);
    mag_clear(r);
    mag_clear(m);
    return 1;
}

#else

int
_acb_elem_small(acb_t r1, acb_t r2, const acb_t z, int f, slong prec)
{
    return 0;
}

#endif
