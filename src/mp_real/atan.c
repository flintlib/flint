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

/* atan of a ball to a relative accuracy of about 2^-prec.

   The midpoint m is evaluated exactly as given; the radius rho enters
   at the end, as rho itself (|atan'| <= 1) or, when the ball stays
   beyond 1 in absolute value, as rho / (|m| - rho)^2 (|atan'(y)| =
   1/(1 + y^2)), and prec is first lowered to what that allows.

   - |m| < 1: m, m - m^3/3 with the remainder in the radius for small
     m; from 32 leading zero bits (by size) the series of
     _mp_real_atan_rs; otherwise the kernel (bitwise up to 600 limbs,
     Newton over the forward functions beyond) at prec + z bits, so
     that the result keeps its relative accuracy near zero;
   - |m| = 1: pi/4;
   - 1 < |m| < 2: pi/4 + atan((|m| - 1)/(|m| + 1)), the argument in
     (0, 1/3);
   - |m| >= 2: pi/2 - atan(1/|m|), 1/|m| <= 1/2.
   The last two take one ball division and pi/4 from the cache; the
   result exceeds pi/4, so absolute accuracy suffices and no
   cancellation arises. */

#define AT_BITWISE_MAX 600

/* the series from z leading zero bits on: the reduced series of
   series_tapered.c / series_rs.c below 32 zero bits, _mp_real_atan_rs
   from 32 (measured end to end, kernel / series time ratios: 1.1-2 at
   1-8 limbs from z = 10, 1.1-1.5 at 10-24 limbs from z = 14-16,
   1.1-1.3 at 32-64 limbs from z = 20-32, 1.05-1.35 at 96-512 limbs from
   z = 40-128) */
static slong
_at_series_min_z(slong n)
{
    if (n <= 8)
        return 10;
    if (n <= 14)
        return 14;
    if (n <= 24)
        return 16;
    if (n <= 40)
        return 20;
    if (n <= 48)
        return 24;
    if (n <= 64)
        return 32;
    if (n <= 96)
        return 40;
    if (n <= 256)
        return 96;
    if (n <= 512)
        return 128;
    return 192;
}

/* kernel limbs for an absolute accuracy of about 2^-p: the bitwise
   bound 4 r + 64, the division and pi/2 (a few ulps) and two bits */
static slong
_at_limbs(slong p)
{
    slong n = (p + 12 + FLINT_BITS - 1) / FLINT_BITS, g;

    if (n <= AT_BITWISE_MAX)
        g = FLINT_BIT_COUNT(4 * (ulong) _mp_real_atan_default_r_inline(n) + 72) + 3;
    else
        g = 8;
    n = (p + g + FLINT_BITS - 1) / FLINT_BITS;
#if FLINT_BITS == 32
    n = FLINT_MAX(n, 2);
#endif
    return n;
}

/* (y, n) = atan(v), v in [0, 1) at n fraction limbs */
static void
_at_kernel(nn_ptr y, ulong * err, nn_srcptr v, slong n)
{
    slong z = _mp_real_elem_lzb(v, n);

    if (z == WORD_MAX)
    {
        flint_mpn_zero(y, n);
        *err = 0;
    }
    else if (z >= _at_series_min_z(n))
    {
        if (z >= 32)
            _mp_real_atan_rs(y, err, v, n);
        else
        {
            _mp_real_series_atan(y, v, n, (flint_bitcnt_t) z, MP_REAL_SERIES_ATAN);
            *err = 6;
        }
    }
    else if (n <= AT_BITWISE_MAX)
    {
        if (!_mp_real_atan_opt(y, err, v, n))
            _mp_real_atan_bitwise_rs(y, err, v, n, 0);
    }
    else
        _mp_real_atan_newton(y, err, v, n);
}

/* (v, n) = 1/|m| for exact |m| >= 2 (m->exp >= 1), within two ulps:
   B^(n + k - e) / M by one approximate reciprocal (flint_mpn_invapprox), M the top k <= n + 1 limbs of the
   mantissa (a relative truncation below B^-n) and e = m->exp */
static void
_at_recip(nn_ptr v, slong n, const mp_real_t m)
{
    slong k = FLINT_MIN(m->size, n + 1), a = n + k - m->exp;
    nn_ptr r;
    TMP_INIT;

    flint_mpn_zero(v, n);
    if (a < k)
        return;                 /* 1/|m| B^n = B^a / M <= 1: v = 0 is
                                   within one ulp (and invapprox wants
                                   a >= k) */

    TMP_START;
    r = TMP_ALLOC((a - k + 3) * sizeof(ulong));
    /* floor(B^a / M) or one more, a - k + 2 <= n + 1 limbs (the top one
       zero for e = 1, since 1/|m| <= 1/2) */
    flint_mpn_invapprox(r, m->d + m->size - k, k, a);
    flint_mpn_copyi(v, r, FLINT_MIN(n, a - k + 2));
    TMP_END;
}

/* atan(m), m exact and nonzero, to about prec bits */
static void
_at_mid(mp_real_t res, const mp_real_t m, slong prec)
{
    slong emid, n;
    int neg = m->negative;

    emid = FLINT_BITS * (m->exp - 1) + FLINT_BIT_COUNT(m->d[m->size - 1]);

    if (m->exp <= 0)
    {
        /* |m| < 2^emid <= 2^-z */
        slong z = -emid;
        slong wp = mp_real_prec_bits(prec + 4);

        if (2 * z >= prec + 3)
        {
            /* |atan m - m| <= |m|^3/3 */
            _mp_real_elem_set_trunc(res, m, wp);
            mp_real_add_error_2exp_si(res, 3 * emid - 1);
            return;
        }

        if (4 * z >= prec + 3)
        {
            /* m - m^3/3 in fixed point at prec + z + 8 bits, remainder
               |m|^5/5 in the radius; the truncations below 4 ulps */
            nn_ptr v, q;
            int trunc;
            TMP_INIT;

            n = (prec + z + 8 + FLINT_BITS - 1) / FLINT_BITS;
            TMP_START;
            v = TMP_ALLOC(3 * n * sizeof(ulong));
            q = v + n;
            trunc = _mp_real_elem_copy(v, n, n, m);
            flint_mpn_sqrhigh(q + n, v, n);
            flint_mpn_mulhigh_n(q, q + n, v, n);
            mpn_divrem_1(q, 0, q, n, 3);
            mp_real_fit_length(res, n + 1);
            mpn_sub_n(res->d, v, q, n);
            res->d[n] = 0;
            _mp_real_elem_finish(res, n, 3 + trunc, neg);
            mp_real_add_error_2exp_si(res, 5 * emid - 2);
            TMP_END;
            return;
        }

        {
            nn_ptr v;
            ulong err;
            int trunc;
            TMP_INIT;

            nn_srcptr vv;

            n = _at_limbs(prec + z);
            TMP_START;
            v = TMP_ALLOC(n * sizeof(ulong));
            vv = _mp_real_elem_frame(v, n, m, res->d, NULL, &trunc);
            mp_real_fit_length(res, n + 1);
            _at_kernel(res->d, &err, vv, n);
            res->d[n] = 0;
            _mp_real_elem_finish(res, n, err + trunc, neg);
            TMP_END;
        }
        return;
    }

    if (m->exp == 1 && m->d[m->size - 1] == 1
        && flint_mpn_zero_p(m->d, m->size - 1))
    {
        n = mp_real_prec_bits(prec + 2);
        mp_real_const_pi4(res, n, 1);
        if (neg)
            mp_real_neg(res, res);
        return;
    }

    if (m->exp == 1 && m->d[m->size - 1] == 1)
    {
        /* 1 < |m| < 2: pi/4 + atan(w), w = (|m| - 1)/(|m| + 1) in
           (0, 1/3), by one approximate division */
        nn_ptr v, y;
        ulong ev, ea;
        TMP_INIT;

        n = _at_limbs(prec);
        TMP_START;
        v = TMP_ALLOC(3 * (n + 1) * sizeof(ulong));
        y = v + n + 1;
        {
            /* s = |m| - 1 (truncated, one ulp), w = s/(2 + s) within 2 */
            nn_ptr t = y + n + 1;
            _mp_real_elem_copy(t, n + 1, n, m);
            _mp_real_elem_half_ratio(v, t, n, 0);
            ev = 2;
        }

        _at_kernel(y, &ea, v, n);

        /* pi/4 within one ulp below */
        mp_real_fit_length(res, n + 1);
        res->d[n] = mpn_add_n(res->d, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, n), y, n);
        _mp_real_elem_finish(res, n, ea + ev + 1, neg);

        TMP_END;
        return;
    }

    {
        /* |m| >= 2: pi/2 - atan(1/|m|), 1/|m| <= 1/2 */
        nn_ptr v, y;
        ulong ev, ea;
        TMP_INIT;

        n = _at_limbs(prec);
        TMP_START;
        v = TMP_ALLOC(2 * (n + 1) * sizeof(ulong));
        y = v + n + 1;
        _at_recip(v, n, m);
        ev = 2;

        _at_kernel(y, &ea, v, n);

        /* P = 2 floor(pi/4 B^n) within 2 ulps below pi/2 */
        mp_real_fit_length(res, n + 1);
        res->d[n] = mpn_lshift(res->d, _mp_real_const_ptr(MP_REAL_CONST_ID_PI4, n), n, 1);
        mpn_sub(res->d, res->d, n + 1, y, n);
        _mp_real_elem_finish(res, n, ea + ev + 2, neg);

        TMP_END;
    }
}

void
mp_real_atan_bits(mp_real_t res, const mp_real_t x, slong prec)
{
    mp_real_struct mid;
    ulong xerr = x->err;
    slong xanc = x->exp - x->size, e;

    prec = FLINT_MAX(prec, 2);

    if (x->size == 0)
    {
        /* |atan y| <= min(|y|, pi/2) */
        if (xerr == 0)
            mp_real_zero(res);
        else
        {
            e = mp_real_abs_bound_lt_2exp_si(x);
            _mp_real_elem_set_error(res, 0, FLINT_MIN(e, 1));
        }
        return;
    }

    mid = *x;
    mid.err = 0;

    if (xerr == 0)
    {
        _at_mid(res, &mid, prec);
        return;
    }
    else
    {
        slong rel = mp_real_rel_radius_lt_2exp_si(x), acc;
        slong emid = FLINT_BITS * (x->exp - 1) + FLINT_BIT_COUNT(x->d[x->size - 1]);
        double rd = (double) xerr, v = 0.0;
        slong a = 0;

        rel = FLINT_MAX(rel, -WORD_MAX / 4);

        /* rho < 2^rel |m|: relative to |atan m| ~ |m| for |m| < 1, and
           rho / (|m| - rho)^2 < 2^(rel - emid + 3) relative to
           atan |m| > pi/4 when the ball stays beyond 1 */
        if (emid <= 0)
            acc = -rel;
        else
            acc = -(rel - emid + 3);
        prec = FLINT_MIN(prec, FLINT_MAX(acc + 8, 10));

        /* the radius, read before res is written (it may alias x): when
           |m| - rho >= 1, rho / (|m| - rho)^2 <= (rd / (mlo (1 - r))^2)
           B^(xanc - 2 (exp - 1)) with |m| >= mlo B^(exp - 1) and
           r >= rho / (mlo B^(exp - 1)); else rho itself.  Double
           arithmetic, rounded up. */
        if (x->exp >= 1)
        {
            double mlo = _mp_real_mag_lo(x), r, lo;
            slong b = xanc - (x->exp - 1);

            /* as for log: r >= 2^64 from b = 2 on; lo >= 1 means
               mlo (1 - r) (1 - 2^-50) B^(exp - 1) >= 1, which for
               exp >= 2 holds whenever r < 1 (then 1 - r >= 2^-53) */
            r = (b >= 2) ? 0x1p64
                : _mp_real_elem_scale_up(rd / mlo, b) * (1.0 + 0x1p-50);
            lo = mlo * (1.0 - r) * (1.0 - 0x1p-50);
            if (r < 1.0 && (x->exp >= 2 || lo >= 1.0))
            {
                double d = mlo * (1.0 - r) * (1.0 - 0x1p-50);
                v = rd / (d * d) * (1.0 + 0x1p-48);
                a = xanc - 2 * (x->exp - 1);
            }
        }

        _at_mid(res, &mid, prec);

        if (v != 0.0)
            _mp_real_elem_add_rad_d(res, v, a);
        else
            _mp_real_elem_add_rad(res, xerr, xanc);
    }
}
