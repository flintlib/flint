/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb/elem.h"

void
_arb_set_mpn_err(arb_t res, nn_ptr y, slong len, slong e, ulong elo, ulong ehi, int sgnbit, slong prec)
{
    ulong err = FLINT_MAX(elo, ehi);
    slong fix;
    int inexact;

    while (len > 0 && y[len - 1] == 0)
        len--;

    if (len == 0)
    {
        arf_zero(arb_midref(res));
        mag_set_ui_2exp_si(arb_radref(res), err, e);
        return;
    }

    /* the midpoint rounded to nearest: the radius gains half an ulp
       rather than one (the kernels' errors are far below an ulp) */
    inexact = _arf_set_round_mpn(arb_midref(res), &fix, y, len, sgnbit, prec, ARF_RND_NEAR);
    _fmpz_demote(ARF_EXPREF(arb_midref(res)));
    ARF_EXP(arb_midref(res)) = e + len * FLINT_BITS + fix;

    if (err == 0)
        mag_zero(arb_radref(res));
    else
        mag_set_ui_2exp_si(arb_radref(res), err, e);

    if (inexact)
        arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec + 1);
}

void
_arb_const_set_floor(arb_t res, int id, int k, slong prec)
{
    slong w, fix;

    /* at least 16 guard bits over prec, within the static tables */
    w = (prec + 16 + FLINT_BITS - 1) / FLINT_BITS;
    FLINT_ASSERT(w <= MP_REAL_CONST_STATIC_N);

    /* the value lies in [Y, Y + 1) units of 2^(k - FLINT_BITS w): the
       midpoint rounded to nearest */
    _arf_set_round_mpn(arb_midref(res), &fix, _mp_real_const_ptr(id, w), w, 0, prec, ARF_RND_NEAR);
    _fmpz_demote(ARF_EXPREF(arb_midref(res)));
    ARF_EXP(arb_midref(res)) = fix + k;

    /* half an ulp: Y and Y + 1 round alike for every prec served here
       (with 64- and 32-bit limbs, checked by arb/test/t-const_floor), so
       that the midpoint is the correctly rounded value */
    _fmpz_demote(MAG_EXPREF(arb_radref(res)));
    MAG_EXP(arb_radref(res)) = ARF_EXP(arb_midref(res)) - prec;
    MAG_MAN(arb_radref(res)) = MAG_ONE_HALF;
}

void
_arb_const_round_cached(arb_t res, void (*cached)(arb_t, slong), slong prec)
{
    /* the cached value at 32 more bits rounded to nearest: within half an
       ulp plus the cached error (arb_set_round would add a whole ulp) */
    cached(res, prec + 32);
    if (arf_set_round(arb_midref(res), arb_midref(res), prec, ARF_RND_NEAR))
        arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec + 1);
}

/* the short series of the small arguments (arf/elem.c), as a ball */
void
_arb_small_arg_series(arb_t z, nn_srcptr d, slong xn, slong e, int sgnbit, int kind, slong prec)
{
    nn_ptr Y;
    slong n = _arb_elem_limbs(prec) + 2, ue;
    ulong elo, ehi;
    int even;
    TMP_INIT;

    TMP_START;
    Y = TMP_ALLOC((n + 2) * sizeof(ulong));
    even = _arf_small_arg_series_frame(Y, &ue, &elo, &ehi, d, xn, e, kind, n);
    _arb_set_mpn_err(z, Y, n + 1, ue, elo, ehi, even ? 0 : sgnbit, prec);
    TMP_END;
}

#if FLINT_BITS == 64
/* The inverse functions at one and two limbs by the register
   compositions of arf/elem_arc_fast.c: the value at the midpoint, plus
   the radius times a bound D for |f'| over the ball (the derivatives
   being monotone in |x| where they are unbounded: asin, acos, atanh at
   |m| + r, acosh at |m| - r; asinh below 1); 0 when the paths
   do not apply (precision or mantissa beyond two limbs, the ball
   reaching the domain's boundary, exponents beyond the paths') */
int
_arb_arc_fast(arb_t res, const arb_t x, int f, slong prec)
{
    _arf_arc_out_struct out;
    nn_srcptr d;
    slong xn, n, e;
    ulong dd[2];
    mag_t D;
    int ok;

    if (prec > 2 * FLINT_BITS || !arb_is_finite(x) || arf_is_special(arb_midref(x))
        || COEFF_IS_MPZ(ARF_EXP(arb_midref(x))) || ARF_SIZE(arb_midref(x)) > 2)
        return 0;

    e = ARF_EXP(arb_midref(x));
    if (e < -39 || e > 62)
        return 0;

    mag_init(D);
    ok = 1;

    if (mag_is_zero(arb_radref(x)))
    {
        /* the domain (the exact |x| = 1 is handled by the callers) */
        if (f == ARF_ARC_ASIN || f == ARF_ARC_ACOS || f == ARF_ARC_ATANH)
            ok = arf_cmpabs_ui(arb_midref(x), 1) < 0;
        else if (f == ARF_ARC_ACOSH)
            ok = arf_cmp_ui(arb_midref(x), 1) > 0;
    }
    else if (mag_cmp_2exp_si(arb_radref(x), e - 20) >= 0)
    {
        /* (a radius beyond 2^-20 |m|: the compositions in ball
           arithmetic do better) */
        ok = 0;
    }
    else if (f == ARF_ARC_ASINH)
    {
        /* 1 / sqrt(1 + l^2) at l = |m| - r (or 0), decreasing in |x|:
           about 1 / |m| for a large argument */
        mag_t l;
        mag_init(l);
        arf_get_mag_lower(l, arb_midref(x));
        mag_sub_lower(l, l, arb_radref(x));
        mag_mul_lower(l, l, l);
        mag_one(D);
        mag_add_lower(l, l, D);
        mag_rsqrt(D, l);
        mag_clear(l);
    }
    else
    {
        /* the domain and the derivative bound in arf arithmetic (the
           30 bits of a mag cannot resolve 1 - u^2 near 1): asin, acos,
           atanh at u = |m| + r < 1, acosh at l = |m| - r > 1 */
        arf_t u, v;
        mag_t t;
        arf_init(u);
        arf_init(v);
        mag_init(t);
        arf_set_mag(v, arb_radref(x));
        arf_abs(u, arb_midref(x));

        if (f != ARF_ARC_ACOSH)
        {
            arf_add(u, u, v, FLINT_BITS, ARF_RND_UP);
            arf_mul(u, u, u, FLINT_BITS, ARF_RND_UP);
            arf_sub_ui(v, u, 1, FLINT_BITS, ARF_RND_UP);        /* 1 - u^2 from below */
            arf_neg(v, v);
        }
        else
        {
            ok = !ARF_SGNBIT(arb_midref(x));
            arf_sub(u, u, v, FLINT_BITS, ARF_RND_DOWN);
            arf_mul(u, u, u, FLINT_BITS, ARF_RND_DOWN);
            arf_sub_ui(v, u, 1, FLINT_BITS, ARF_RND_DOWN);     /* l^2 - 1 from below */
        }

        if (arf_sgn(v) <= 0)
            ok = 0;
        else
        {
            arf_get_mag_lower(t, v);
            if (f == ARF_ARC_ATANH)
                mag_inv(D, t);          /* 1 / (1 - u^2) */
            else
                mag_rsqrt(D, t);        /* 1 / sqrt(1 - u^2), 1 / sqrt(l^2 - 1) */
        }

        arf_clear(u);
        arf_clear(v);
        mag_clear(t);
    }

    if (ok)
    {
        ARF_GET_MPN_READONLY(d, xn, arb_midref(x));
        n = (xn == 1 && prec <= FLINT_BITS) ? 1 : 2;
        dd[0] = (xn == n) ? d[0] : 0;
        dd[1] = d[xn - 1];

        if (_arf_arc_fast(&out, dd, n, e, ARF_SGNBIT(arb_midref(x)), f, 0) == ARF_ELEM_UNCERTAIN)
            ok = 0;
        else
        {
            mag_t r;
            mag_init_set(r, arb_radref(x));     /* (res may alias x) */
            if (n == 1)
                _arb_small_set_1(res, out.T[2], out.T[1], out.T[0], out.k, out.err, out.sgnbit, prec);
            else
                _arb_small_set_2(res, out.T[3], out.T[2], out.T[1], out.T[0], out.k, out.err, out.sgnbit, prec);
            if (!mag_is_zero(r))
                mag_addmul(arb_radref(res), r, D);
            mag_clear(r);
        }
    }

    mag_clear(D);
    return ok;
}
#else
int
_arb_arc_fast(arb_t res, const arb_t x, int f, slong prec)
{
    return 0;
}
#endif

#if FLINT_BITS == 64
#include "acf/elem.h"

/* log |a + bi| (which = 0) or atan2(b, a) (which = 1) at one and two
   limbs by the register paths shared with acf (acf/elem.h), for nonzero
   midpoints of at most two limbs and radii small against the parts
   (and, for atan2, a ball off the negative real axis); the midpoint
   through the frame, the radii propagated through the gradients
   (a, b) / |z|^2 resp. (-b, a) / |z|^2 over the ball.  Returns 0 when
   the path does not apply. */
int
_arb_log_hypot_atan2_small(arb_t res, const arb_t a, const arb_t b, int which, slong prec)
{
    arf_srcptr am = arb_midref(a), bm = arb_midref(b);
    _acf_frame_struct F;
    nn_srcptr ad, bd;
    slong an, bn, ea, eb;
    ulong da1, da0, db1, db0;
    int two, st;
    mag_t r;

    if (prec > 2 * FLINT_BITS || arf_is_special(am) || arf_is_special(bm)
        || ARF_SIZE(am) > 2 || ARF_SIZE(bm) > 2
        || COEFF_IS_MPZ(ARF_EXP(am)) || COEFF_IS_MPZ(ARF_EXP(bm)))
        return 0;

    ea = ARF_EXP(am);
    eb = ARF_EXP(bm);
    if (FLINT_ABS(ea) >= (WORD(1) << 40) || FLINT_ABS(eb) >= (WORD(1) << 40))
        return 0;

    mag_init(r);
    mag_add(r, arb_radref(a), arb_radref(b));
    if (!mag_is_finite(r) || (!mag_is_zero(r) && mag_cmp_2exp_si(r, FLINT_MAX(ea, eb) - 3) >= 0)
        || (which == 1 && ARF_SGNBIT(am) && arf_cmpabs_mag(bm, arb_radref(b)) <= 0))
    {
        mag_clear(r);
        return 0;
    }

    ARF_GET_MPN_READONLY(ad, an, am);
    ARF_GET_MPN_READONLY(bd, bn, bm);
    two = !(prec <= FLINT_BITS && an == 1 && bn == 1);
    da1 = ad[an - 1];
    da0 = (an == 2) ? ad[0] : 0;
    db1 = bd[bn - 1];
    db0 = (bn == 2) ? bd[0] : 0;

    if (which == 0)
        st = two ? _acf_small_log_abs_2(&F, da1, da0, ea, db1, db0, eb) : _acf_small_log_abs_1(&F, da1, ea, db1, eb);
    else
        st = two ? _acf_small_atan2_2(&F, da1, da0, ea, ARF_SGNBIT(am), db1, db0, eb, ARF_SGNBIT(bm))
                 : _acf_small_atan2_1(&F, da1, ea, ARF_SGNBIT(am), db1, eb, ARF_SGNBIT(bm));
    if (st != 0 || (two ? ((F.y[3] | F.y[2]) == 0) : ((F.y[2] | F.y[1]) == 0)))
    {
        mag_clear(r);
        return 0;
    }

    {
        arb_t t;
        arb_init(t);
        if (two)
            _arb_small_set_2(t, F.y[3], F.y[2], F.y[1], F.y[0], F.k, F.err, F.sgnbit, prec);
        else
            _arb_small_set_1(t, F.y[2], F.y[1], F.y[0], F.k, F.err, F.sgnbit, prec);

        if (!mag_is_zero(r))
        {
            mag_t lz, m1, m2, m;
            mag_init(lz);
            mag_init(m1);
            mag_init(m2);
            mag_init(m);
            arf_get_mag_lower(lz, am);
            mag_mul_lower(lz, lz, lz);
            arf_get_mag_lower(m, bm);
            mag_mul_lower(m, m, m);
            mag_add_lower(lz, lz, m);
            mag_sqrt_lower(lz, lz);
            mag_sub_lower(lz, lz, r);
            mag_mul_lower(lz, lz, lz);
            arf_get_mag(m1, am);
            mag_add(m1, m1, arb_radref(a));
            arf_get_mag(m2, bm);
            mag_add(m2, m2, arb_radref(b));
            if (which == 0)
            {
                mag_mul(m, m1, arb_radref(a));
                mag_addmul(m, m2, arb_radref(b));
            }
            else
            {
                mag_mul(m, m2, arb_radref(a));
                mag_addmul(m, m1, arb_radref(b));
            }
            mag_div(m, m, lz);
            mag_add(arb_radref(t), arb_radref(t), m);
            mag_clear(lz);
            mag_clear(m1);
            mag_clear(m2);
            mag_clear(m);
        }
        arb_swap(res, t);
        arb_clear(t);
    }

    mag_clear(r);
    return 1;
}
#endif
