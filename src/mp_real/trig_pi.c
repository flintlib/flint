/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* sin(pi x), cos(pi x) and tan(pi x) of a ball.

   The midpoint m is dyadic, so its reduction is exact: from the units
   bit and the fraction bits of |m|,

       |m| = a/2 + s t,   a in 0..3, s = +-1, t in [0, 1/4],

   and with v = pi t in [0, pi/4], exactly as for mp_real_sin_cos_bits
   and mp_real_tan_bits after their reduction mod pi/2,

       sin(pi |m|) = s sin v, cos v, -s sin v, -cos v   (a = 0, 1, 2, 3),
       cos(pi |m|) = cos v, -s sin v, -cos v, s sin v,
       tan(pi |m|) = s tan v (a even), -s / tan v (a odd).

   t = 0 gives the exact values 0 and +-1 (and a pole of tan for odd a),
   t = 1/4 the exact tan(pi/4) = 1.  Otherwise t has z leading zero bits,
   and v = (pi/4) 4t goes to the kernel at z more bits of absolute
   precision (the small value then keeps its relative accuracy; the
   kernels' series make small arguments cheap), or for 2z >= prec + 3 to
   the first-order values sin v = tan v = v, cos v = 1 with the Taylor
   remainders in the radius: nothing limits the relative accuracy near
   the zeros of sin and cos or the zeros and poles of tan, at any
   magnitude of m.  t lives in a buffer of the size of m's mantissa and
   everything is written into the outputs directly (temporary balls cost
   more than the evaluation at one or two limbs).

   The radius r of x enters as a radius pi r of v: added to the sine and
   cosine (whose derivatives are bounded by 1), through the pole test
   of mp_real_tan_bits for the tangent. */

/* the exact reduction of |m| (sign ignored): |m| = a/2 + s t, returning
   a + 4 [s < 0], with t = (T, *tl) B^(*te - *tl), T's top limb nonzero
   (*tl = 0 for t = 0), *te <= 0.  T has room for m->size limbs. */
static int
_pi_reduce(nn_ptr T, slong * tl, slong * te, const mp_real_t m)
{
    slong size = m->size, e = m->exp, fl, i;
    ulong half;
    int a, fold;

    /* an integer multiple of B (or zero): even */
    if (size == 0 || e - size >= 1)
    {
        *tl = 0;
        *te = 0;
        return 0;
    }

    if (e < 0)
    {
        /* |m| < B^-1: t = |m| */
        flint_mpn_copyi(T, m->d, size);
        *tl = size;
        *te = e;
        return 0;
    }

    /* the fraction of |m| as the fl limbs below the units limb (the
       frame te = 0), and a = floor(2 |m|) mod 4 */
    fl = size - e;
    a = (e >= 1) ? (int) (2 * (m->d[fl] & 1)) : 0;
    if (fl == 0)
    {
        *tl = 0;
        *te = 0;
        return a;
    }
    flint_mpn_copyi(T, m->d, fl);

    half = T[fl - 1] >> (FLINT_BITS - 1);
    a += (int) half;
    T[fl - 1] &= ~(UWORD(1) << (FLINT_BITS - 1));

    /* fold t > 1/4: t = 1/2 - t, a + 1, s = -1 */
    fold = 0;
    if (T[fl - 1] >> (FLINT_BITS - 2))
    {
        if (T[fl - 1] != (UWORD(1) << (FLINT_BITS - 2)))
            fold = 1;
        else
            for (i = 0; i < fl - 1 && !fold; i++)
                fold = (T[i] != 0);
    }
    if (fold)
    {
        /* 1/2 - t = B^fl / 2 - T */
        mpn_neg(T, T, fl);
        T[fl - 1] += UWORD(1) << (FLINT_BITS - 1);
        a = (a + 1) & 3;
    }

    /* the fraction's leading zero limbs move its exponent */
    *tl = fl;
    while (*tl > 0 && T[*tl - 1] == 0)
        (*tl)--;
    *te = *tl - fl;
    return (a & 3) | (fold << 2);
}

/* 2^et bounds t (2^(et - 1) <= t < 2^et), t nonzero */
FLINT_FORCE_INLINE slong
_pi_t_exp(nn_srcptr T, slong tl, slong te)
{
    return FLINT_BITS * te - (slong) flint_clz(T[tl - 1]);
}

/* (V, N) = floor-ish (pi/4) 4t B^N, below the exact value by less than 3
   ulps: P = floor(pi/4 B^N) times the top K <= N + 1 limbs of T, times 4,
   the product truncated (P, T's truncation and the floor each under one
   ulp since 4t <= 1) */
static void
_pi_v_fixed(nn_ptr V, nn_srcptr T, slong tl, slong te, slong N)
{
    nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, N);
    slong K = FLINT_MIN(tl, N + 1), sh;
    nn_ptr R;
    TMP_INIT;

    /* V = floor(4 P Tk B^(te - K)): the limbs from K - te on of 4 P Tk */
    sh = K - te;
    if (sh >= N + K + 1)
    {
        flint_mpn_zero(V, N);
        return;
    }

    TMP_START;
    R = TMP_ALLOC((N + K + 1) * sizeof(ulong));
    if (N >= K)
        flint_mpn_mul(R, P, N, T + (tl - K), K);
    else
        flint_mpn_mul(R, T + (tl - K), K, P, N);
    R[N + K] = mpn_lshift(R, R, N + K, 2);

    /* limbs [sh, N + K + 1) hold the value's integral part times B^N, of
       which N limbs (4t <= 1 and pi/4 < 1: below B^N) */
    flint_mpn_zero(V, N);
    flint_mpn_copyi(V, R + sh, FLINT_MIN(N, N + K + 1 - sh));
    TMP_END;
}

/* res = (-1)^neg v, v = (pi/4) 4t, as a ball of about w limbs: P (w
   limbs) times the top K <= w + 1 limbs of T, times 4, exact, within
   4 (Tk + 1) units of the product's lowest limb for P's floor, plus
   4 (P + 1) when T was truncated (then K = w + 1): below 16 B^w units
   either way, 16 at the anchor B^(te - w).  Y, ey, err: the same value
   truncated there (for the reciprocal), Y with room for w + 1 limbs. */
static void
_pi_v_ball(mp_real_t res, nn_ptr Y, slong * L, slong * ey, ulong * err,
    nn_srcptr T, slong tl, slong te, slong w, int neg)
{
    nn_srcptr P = _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, w);
    slong K = FLINT_MIN(tl, w + 1);
    nn_ptr R;
    TMP_INIT;

    TMP_START;
    R = TMP_ALLOC((w + K + 1) * sizeof(ulong));
    if (w >= K)
        flint_mpn_mul(R, P, w, T + (tl - K), K);
    else
        flint_mpn_mul(R, T + (tl - K), K, P, w);
    R[w + K] = mpn_lshift(R, R, w + K, 2);

    /* v = R B^(te - K - w) within 16 B^(te - w) */
    if (res != NULL)
    {
        _mp_real_set_mpn_2exp(res, R, w + K + 1, FLINT_BITS * (te - K - w));
        _mp_real_elem_add_rad(res, 16, te - w);
        if (neg)
            res->negative = 1;
    }

    if (Y != NULL)
    {
        /* R truncated below B^K (the anchor): within 17 units */
        *L = w + 1;
        flint_mpn_copyi(Y, R + K, w + 1);
        *ey = te - w;
        *err = 17;
    }
    TMP_END;
}

/* the ball [0 +- pi r] for the radius xerr B^xanc (two limbs) */
static void
_pi_rad(mp_real_t E, ulong xerr, slong xanc)
{
    mp_real_t P;
    mp_real_init(P);
    mp_real_set_ui(E, xerr);
    mp_real_mul_2exp_si(E, E, FLINT_BITS * xanc);
    mp_real_const_pi4(P, 2, 1);
    mp_real_mul(E, E, P, 2);
    mp_real_mul_2exp_si(E, E, 2);
    mp_real_clear(P);
}

/* the relative precision of the small output for an inexact x: at most
   that of pi r (r < 2^ra) against v < 2^-z, a few bits more */
FLINT_FORCE_INLINE slong
_pi_prec(slong prec, ulong xerr, slong xanc, slong z)
{
    if (xerr == 0)
        return prec;
    return FLINT_MIN(prec, FLINT_MAX(-(FLINT_BITS * xanc
        + (slong) FLINT_BIT_COUNT(xerr) + 2) - z + 8, 10));
}

void
mp_real_sin_cos_pi_bits(mp_real_t rs, mp_real_t rc, const mp_real_t x, slong prec)
{
    ulong xerr = x->err;
    slong xanc = x->exp - x->size, tl, te, z, et;
    int code, a, s, xneg = x->negative, negs, negc;
    mp_real_ptr ov, oc;
    nn_ptr T;
    TMP_INIT;

    prec = FLINT_MAX(prec, 2);

    /* x = 0 exactly */
    if (x->size == 0 && x->err == 0)
    {
        if (rs != NULL)
            mp_real_zero(rs);
        if (rc != NULL)
            mp_real_set_ui(rc, 1);
        return;
    }

    /* a radius pi r >= 1/4 (the absolute accuracy of both outputs below
       two bits): [0 +- 1] */
    if (xerr != 0 && FLINT_BITS * xanc + (slong) FLINT_BIT_COUNT(xerr) + 2 > -2)
    {
        if (rs != NULL)
            _mp_real_elem_set_error(rs, 0, 0);
        if (rc != NULL)
            _mp_real_elem_set_error(rc, 0, 0);
        return;
    }

    /* (the outputs may alias x: x is read by the reduction only) */
    TMP_START;
    T = TMP_ALLOC(FLINT_MAX(x->size, 1) * sizeof(ulong));
    code = _pi_reduce(T, &tl, &te, x);
    a = code & 3;
    s = (code & 4) ? -1 : 1;

    /* sin v goes to ov, cos v to oc, with the signs
       sin = s sin v, cos v, -s sin v, -cos v;
       cos = cos v, -s sin v, -cos v, s sin v   (a = 0, 1, 2, 3) */
    ov = (a & 1) ? rc : rs;
    oc = (a & 1) ? rs : rc;
    if (a & 1)
    {
        negs = (a == 3);
        negc = (s < 0) ^ (a == 1);
    }
    else
    {
        negs = (s < 0) ^ (a == 2);
        negc = (a == 2);
    }
    negs ^= xneg;

    if (tl == 0)
    {
        /* sin v = 0, cos v = 1 exactly */
        if (ov != NULL)
            mp_real_zero(ov);
        if (oc != NULL)
        {
            mp_real_set_ui(oc, 1);
            if ((a & 1) ? negs : negc)
                mp_real_neg(oc, oc);
        }
    }
    else
    {
        /* v < 2^(et + 2), z >= 0 */
        et = _pi_t_exp(T, tl, te);
        z = FLINT_MAX(-(et + 2), 0);
        prec = _pi_prec(prec, xerr, xanc, z);

        if (2 * z >= prec + 3)
        {
            /* sin v = v within |v|^3/6 < 2^(3e - 2), 0 <= 1 - cos v <=
               v^2/2 < 2^(2e - 1), e = -z */
            slong e = -z;
            if (ov != NULL)
            {
                _pi_v_ball(ov, NULL, NULL, NULL, NULL, T, tl, te,
                    mp_real_prec_bits(prec), (a & 1) ? negc : negs);
                mp_real_add_error_2exp_si(ov, _mp_real_err_exp_clamp(3 * e - 2, e, prec));
            }
            if (oc != NULL)
            {
                mp_real_set_ui(oc, 1);
                mp_real_add_error_2exp_si(oc, _mp_real_err_exp_clamp(2 * e - 1, 1, prec));
                if ((a & 1) ? negs : negc)
                    mp_real_neg(oc, oc);
            }
        }
        else
        {
            slong N = _mp_real_sin_cos_limbs(prec + z);
            nn_ptr V = TMP_ALLOC(N * sizeof(ulong));
            _pi_v_fixed(V, T, tl, te, N);
            _mp_real_sin_cos_eval(rs, rc, V, N, 3, a, s, xneg);
        }
    }

    TMP_END;

    /* the radius of x, pi r in both outputs */
    if (xerr != 0)
    {
        mp_real_t E;
        mp_real_init(E);
        _pi_rad(E, xerr, xanc);
        if (rs != NULL)
            _mp_real_elem_add_mag(rs, E);
        if (rc != NULL)
            _mp_real_elem_add_mag(rc, E);
        mp_real_clear(E);
    }
}

int
mp_real_tan_pi_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    ulong xerr = x->err;
    slong xanc = x->exp - x->size, tl, te, z, et, wn;
    int code, a, s, neg, ok = 1;
    nn_ptr T;
    TMP_INIT;

    prec = FLINT_MAX(prec, 2);

    /* x = 0 exactly */
    if (x->size == 0 && x->err == 0)
    {
        mp_real_zero(res);
        return 1;
    }

    TMP_START;
    T = TMP_ALLOC(FLINT_MAX(x->size, 1) * sizeof(ulong));
    code = _pi_reduce(T, &tl, &te, x);
    a = code & 3;
    s = (code & 4) ? -1 : 1;

    /* tan(pi m) = s tan v (a even), -s / tan v (a odd), negated for
       m < 0 */
    neg = ((a & 1) ? (s > 0) : (s < 0)) ^ x->negative;

    if (tl == 0)
    {
        /* tan v = 0: a zero of tan, or a pole */
        if (a & 1)
            ok = 0;
        else
            mp_real_zero(res);
    }
    else if (te == 0 && T[tl - 1] == (UWORD(1) << (FLINT_BITS - 2))
        && flint_mpn_zero_p(T, tl - 1))
    {
        /* t = 1/4: tan v = 1/tan v = 1 */
        mp_real_set_ui(res, 1);
        if (neg)
            mp_real_neg(res, res);
    }
    else
    {
        et = _pi_t_exp(T, tl, te);
        z = FLINT_MAX(-(et + 2), 0);
        prec = _pi_prec(prec, xerr, xanc, z);
        wn = mp_real_prec_bits(prec);

        if (2 * z >= prec + 3)
        {
            /* tan v = v (1 + d), |d| < 0.372 v^2 < 2^(2e - 1), e = -z */
            slong e = -z;
            if (!(a & 1))
            {
                _pi_v_ball(res, NULL, NULL, NULL, NULL, T, tl, te, wn, neg);
                mp_real_add_error_2exp_si(res, _mp_real_err_exp_clamp(3 * e - 1, e, prec));
            }
            else
            {
                slong L, ey;
                ulong err;
                nn_ptr Y = TMP_ALLOC((wn + 2) * sizeof(ulong));
                _pi_v_ball(NULL, Y, &L, &ey, &err, T, tl, te, wn, 0);
                ok = _mp_real_inv_mpn(res, Y, L, ey, err, 2 * e - 1, wn, neg);
            }
        }
        else
        {
            /* the kernel within its bound, the 3 ulps of v times
               sec^2 v <= 2 */
            slong N = _mp_real_tan_limbs(prec + z);
            nn_ptr V, y;
            ulong err;

            V = TMP_ALLOC(N * sizeof(ulong));
            y = TMP_ALLOC((N + 1) * sizeof(ulong));
            _pi_v_fixed(V, T, tl, te, N);
            _mp_real_tan_kernel(y, &err, V, N);
            err += 7;

            if (!(a & 1))
            {
                mp_real_fit_length(res, N + 1);
                flint_mpn_copyi(res->d, y, N + 1);
                _mp_real_elem_finish(res, N, err, neg);
            }
            else
                ok = _mp_real_inv_mpn(res, y, N + 1, -N, err, WORD_MIN, wn, neg);
        }
    }

    TMP_END;

    /* the radius of x, a radius pi r of the argument of tan */
    if (ok && xerr != 0)
    {
        mp_real_t E;
        mp_real_init(E);
        _pi_rad(E, xerr, xanc);
        ok = _mp_real_tan_add_arg_rad(res, E);
        mp_real_clear(E);
    }

    if (!ok)
        mp_real_zero(res);

    return ok;
}

/* the wrappers for the other files (see impl.h) */

int
_mp_real_trig_pi_reduce(nn_ptr T, slong * tl, slong * te, const mp_real_t x)
{
    return _pi_reduce(T, tl, te, x);
}

void
_mp_real_trig_pi_v_fixed(nn_ptr V, nn_srcptr T, slong tl, slong te, slong N)
{
    _pi_v_fixed(V, T, tl, te, N);
}
