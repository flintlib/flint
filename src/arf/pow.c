/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "arf/elem.h"

/*
    arf_pow: x^y = exp(y log |x|) by the kernel _mp_real_pow_fp
    (mp_real.h) with the sign (-1)^y for x < 0 and an integer y (nan
    otherwise), after the special values and the cases with an exact
    result (arf/elem.h explains the rounding).

    A dyadic x = m 2^s (m odd) to a dyadic power y = n 2^-f (n odd, f >= 1)
    is rational exactly when m is a 2^f-th power N^(2^f) and 2^f divides
    s (then x^y = N^n 2^(s n / 2^f)); for an integer y = n the same with
    N = m.  Such a result is dyadic for n > 0, with N^n odd: it is
    computed exactly (and rounded once) when N^n has at most about prec
    bits, which covers the exact results and the ties the retries of the
    general code could never decide; longer N^n (whose lowest bit is set)
    and the non-dyadic results for n < 0 are decided by the general code.
    Powers of two to integer powers are exact for any y that fits a word.
*/

/* the exact cases: 1 with the inexact flag or 0 when x^y was computed,
   ARF_ELEM_UNCERTAIN when the result is not known to be rational or is
   too long to compute exactly */
static int
_arf_pow_exact(arf_t res, const arf_t x, const arf_t y, int rsgn, slong prec, arf_rnd_t rnd)
{
    fmpz_t m, s, n, sy, p;
    slong f, mbits, i;
    int status = ARF_ELEM_UNCERTAIN;

    fmpz_init(m);
    fmpz_init(s);
    fmpz_init(n);
    fmpz_init(sy);
    fmpz_init(p);

    arf_get_fmpz_2exp(m, s, x);
    arf_get_fmpz_2exp(n, sy, y);
    fmpz_abs(m, m);

    /* y = n 2^sy with n odd: an integer for sy >= 0 */
    if (fmpz_sgn(sy) >= 0)
    {
        /* (y beyond a word: only a power of two can be exact, when its
           exponent s y is not absurdly long) */
        if (!fmpz_fits_si(sy) || fmpz_get_si(sy) > (WORD(1) << 24) ||
            (fmpz_bits(n) + fmpz_get_si(sy) > FLINT_BITS - 2 && !fmpz_is_one(m)))
            goto cleanup;
        fmpz_mul_2exp(n, n, fmpz_get_si(sy));
        f = 0;
    }
    else
    {
        if (!fmpz_fits_si(sy) || (!fmpz_is_one(m) && -fmpz_get_si(sy) > fmpz_bits(m) + 1))
            goto cleanup;
        f = -fmpz_get_si(sy);

        /* 2^f | s (s != 0 for a power of two, |x| != 1) */
        if (fmpz_val2(s) < (ulong) f && !fmpz_is_zero(s))
            goto cleanup;
        fmpz_fdiv_q_2exp(s, s, f);

        /* m = N^(2^f): f exact square roots (m = 1 trivially) */
        for (i = 0; i < f && !fmpz_is_one(m); i++)
        {
            fmpz_sqrtrem(m, p, m);
            if (!fmpz_is_zero(p))
                goto cleanup;
        }
    }

    /* x^y = (-1)^rsgn N^n 2^(s n), N = m */
    mbits = fmpz_bits(m);

    if (mbits == 1)
    {
        /* a power of two: 2^(s n) */
        fmpz_mul(s, s, n);
        arf_one(res);
        arf_mul_2exp_fmpz(res, res, s);
        status = 0;
    }
    else if (fmpz_bits(n) <= FLINT_BITS - 2 && FLINT_ABS(fmpz_get_si(n)) <= (prec + FLINT_BITS) / (mbits - 1))
    {
        slong ni = fmpz_get_si(n);

        /* (the sign goes in before the rounding, for the directed modes) */
        fmpz_pow_ui(p, m, FLINT_ABS(ni));
        if (rsgn)
            fmpz_neg(p, p);
        fmpz_mul_si(s, s, ni);

        if (ni > 0)
        {
            status = arf_set_round_fmpz(res, p, prec, arf_rnd_relaxed_to_strict(rnd));
            arf_mul_2exp_fmpz(res, res, s);
        }
        else
        {
            /* 1 / (N^|n| 2^(-s n)) */
            arf_t t, u;
            arf_init(t);
            arf_init(u);
            arf_set_fmpz(t, p);
            fmpz_neg(s, s);
            arf_mul_2exp_fmpz(t, t, s);
            arf_one(u);
            status = arf_div(res, u, t, prec, rnd);
            arf_clear(t);
            arf_clear(u);
        }
        rsgn = 0;
    }
    else
        goto cleanup;

    if (rsgn)
        arf_neg(res, res);

cleanup:
    fmpz_clear(m);
    fmpz_clear(s);
    fmpz_clear(n);
    fmpz_clear(sy);
    fmpz_clear(p);
    return status;
}

#if FLINT_BITS == 64
/* the in-register paths at one and two limbs (arf/elem.h), shorter
   mantissas padded: ARF_ELEM_UNCERTAIN when they do not apply or do not
   decide the rounding */
static int
_arf_pow_small(arf_t res, const arf_t x, const arf_t y, int rsgn, slong prec, arf_rnd_t rnd)
{
    nn_srcptr xd, yd;
    slong xn, yn, k;
    ulong y3 = 0, y2 = 0, y1 = 0, y0 = 0, err = 0;

    if (prec > 2 * FLINT_BITS || COEFF_IS_MPZ(ARF_EXP(x)) || COEFF_IS_MPZ(ARF_EXP(y)))
        return ARF_ELEM_UNCERTAIN;

    ARF_GET_MPN_READONLY(xd, xn, x);
    ARF_GET_MPN_READONLY(yd, yn, y);
    if (xn > 2 || yn > 2)
        return ARF_ELEM_UNCERTAIN;

    if (xn == 1 && yn == 1 && prec <= FLINT_BITS)
    {
        if (_arf_small_pow_1(&y2, &y1, &y0, &k, &err, xd[0], ARF_EXP(x), yd[0], ARF_EXP(y), ARF_SGNBIT(y)) == ARF_ELEM_UNCERTAIN)
            return ARF_ELEM_UNCERTAIN;
        return _arf_small_set_1(res, y2, y1, y0, k, err, rsgn, prec, rnd);
    }
    else
    {
        if (_arf_small_pow_2(&y3, &y2, &y1, &y0, &k, &err, xd[xn - 1], (xn == 2) ? xd[0] : 0, ARF_EXP(x),
                yd[yn - 1], (yn == 2) ? yd[0] : 0, ARF_EXP(y), ARF_SGNBIT(y)) == ARF_ELEM_UNCERTAIN)
            return ARF_ELEM_UNCERTAIN;
        return _arf_small_set_2(res, y3, y2, y1, y0, k, err, rsgn, prec, rnd);
    }
}
#else
static int
_arf_pow_small(arf_t res, const arf_t x, const arf_t y, int rsgn, slong prec, arf_rnd_t rnd)
{
    return ARF_ELEM_UNCERTAIN;
}
#endif

/* x^y for finite nonzero x and y, not aliased with res */
static int
_arf_pow(arf_t res, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
{
    nn_srcptr xd, yd;
    slong xn, yn, xe, ye, nk, extra;
    int isint, odd, rsgn, status;

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_pow: exact precision not supported\n");

    isint = arf_is_int(y);
    if (ARF_SGNBIT(x) && !isint)
    {
        arf_nan(res);
        return 0;
    }
    odd = isint && !arf_is_int_2exp_si(y, 1);
    rsgn = ARF_SGNBIT(x) && odd;

    /* |x| = 1 */
    if (arf_cmpabs_ui(x, 1) == 0)
    {
        arf_set_si(res, rsgn ? -1 : 1);
        return 0;
    }

    /* y = +-1, 2, +-1/2 */
    if (arf_is_int_2exp_si(y, -1) && arf_cmpabs_2exp_si(y, 1) <= 0)
    {
        if (arf_equal_si(y, 1))
            return arf_set_round(res, x, prec, arf_rnd_relaxed_to_strict(rnd));
        if (arf_equal_si(y, -1))
            return arf_ui_div(res, 1, x, prec, rnd);
        if (arf_equal_si(y, 2))
            return arf_mul(res, x, x, prec, rnd);
        if (arf_equal_si(y, -2))
        {
            arf_t t;
            int inexact;
            arf_init(t);
            arf_mul(t, x, x, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_ui_div(res, 1, t, prec, rnd);
            arf_clear(t);
            return inexact;
        }
        if (arf_cmpabs_2exp_si(y, -1) == 0)
            return ARF_SGNBIT(y) ? arf_rsqrt(res, x, prec, rnd) : arf_sqrt(res, x, prec, rnd);
    }

    /* the relaxed modes need no exact cases: the fast paths first */
    if ((rnd >= ARF_RND_FAST))
    {
        status = _arf_pow_small(res, x, y, rsgn, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    status = _arf_pow_exact(res, x, y, rsgn, prec, rnd);
    if (status != ARF_ELEM_UNCERTAIN)
        return status;

    ARF_GET_MPN_READONLY(xd, xn, x);
    ARF_GET_MPN_READONLY(yd, yn, y);
    xe = ARF_EXP(x);
    ye = ARF_EXP(y);

    /* (strict modes: the exact cases were settled above) */
    if (!(rnd >= ARF_RND_FAST))
    {
        status = _arf_pow_small(res, x, y, rsgn, prec, rnd);
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    /* bounds on t = y log |x| for x = M 2^E (M in [1/2, 1), |x| != 1):
       |t| < 2^(ye + bits(|E| + 1)) since |log |x|| < |E| + 1, and for
       |E| >= 2, |t| > 2^(ye + bits(|E| - 1) - 3) since |log |x|| >
       (|E| - 1) log 2 and |y| >= 2^(ye - 1), for |E| <= 1 from the
       leading bits of x (_mp_real_pow_t_exceeds).  A tiny t gives 1 + delta,
       |delta| < 2|t|, with the sign of t (log |x| < 0 for E <= 0); a
       huge one is beyond the limit of exp (nan), as are the arguments
       the kernel does not take (|E| >= 2^(FLINT_BITS - 3), a bignum
       exponent of y) */
    {
        fmpz_t E, thi, tlo;
        int tiny, huge, big_e;

        fmpz_init(E);
        fmpz_init(thi);
        fmpz_init(tlo);
        fmpz_abs(E, ARF_EXPREF(x));
        fmpz_add_ui(thi, E, 1);
        fmpz_set_ui(thi, fmpz_bits(thi));
        fmpz_add(thi, thi, ARF_EXPREF(y));
        tiny = (fmpz_cmp_si(thi, -prec - 2) <= 0);
        big_e = (fmpz_cmp_ui(E, 2) >= 0);
        fmpz_sub_ui(tlo, E, 1);
        fmpz_set_ui(tlo, big_e ? fmpz_bits(tlo) : 0);
        fmpz_add(tlo, tlo, ARF_EXPREF(y));
        fmpz_sub_ui(tlo, tlo, 3);
        huge = big_e ? (fmpz_cmp_si(tlo, ARF_EXP_MAGLIM(prec)) >= 0)
            : (!COEFF_IS_MPZ(ye) && _mp_real_pow_t_exceeds(xd, xn, xe, ye, ARF_EXP_MAGLIM(prec)));
        fmpz_clear(E);
        fmpz_clear(thi);
        fmpz_clear(tlo);

        if (tiny)
            return _arf_one_plus_tiny(res, rsgn, ARF_SGNBIT(y) ^ (fmpz_sgn(ARF_EXPREF(x)) <= 0), prec, rnd);

        if (huge || COEFF_IS_MPZ(ye))
        {
            arf_nan(res);
            return 0;
        }

        if (COEFF_IS_MPZ(xe) || xe >= (WORD(1) << (FLINT_BITS - 3)) || xe <= -(WORD(1) << (FLINT_BITS - 3)))
        {
            /* a huge exponent of x with a moderate t (y small): log x by
               arf_log (any exponent) within an ulp at tb + 4 more bits,
               so that t = y L (exact) is within 2^-(wp + 3) of the true
               value and the exponential within 2^-(wp + 2) relatively
               (one unit of its mantissa at wp bits), plus its own ulp */
            arf_t L, T, R;
            slong tb = fmpz_get_si(ARF_EXPREF(y)) + fmpz_bits(ARF_EXPREF(x)) + 2;    /* |t| < 2^tb */

            arf_init(L);
            arf_init(T);
            arf_init(R);

            for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
            {
                nn_srcptr rd;
                nn_ptr out;
                slong rn, m, wp;
                TMP_INIT;

                m = _arf_elem_limbs(prec, rnd, extra);
                wp = FLINT_BITS * m;
                arf_abs(L, x);
                arf_log(L, L, wp + FLINT_MAX(tb, 0) + 4, ARF_RND_FAST);
                arf_mul(T, L, y, ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_exp(R, T, wp, ARF_RND_FAST);

                /* (the test above bounds |t| from below only: a t beyond
                   the limit of exp gives nan, as documented) */
                if (arf_is_special(R))
                {
                    arf_nan(res);
                    status = 0;
                    break;
                }

                ARF_GET_MPN_READONLY(rd, rn, R);
                TMP_START;
                out = TMP_ALLOC((m + 1) * sizeof(ulong));
                _flint_mpn_zero_small(out, m - rn);
                flint_mpn_copyi(out + m - rn, rd, rn);
                /* (the exponent of R may be a bignum: rounded with a zero
                   unit exponent, then shifted) */
                {
                    fmpz_t ue;
                    fmpz_init(ue);
                    fmpz_sub_ui(ue, ARF_EXPREF(R), FLINT_BITS * m);
                    status = _arf_set_mpn_err(res, out, m, 0, ARF_ELEM_EB(3), ARF_ELEM_EB(3), rsgn, prec, rnd);
                    if (status != ARF_ELEM_UNCERTAIN && !arf_is_special(res))
                        fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), ue);
                    fmpz_clear(ue);
                }
                TMP_END;

                if (status != ARF_ELEM_UNCERTAIN)
                    break;
            }

            arf_clear(L);
            arf_clear(T);
            arf_clear(R);
            return status;
        }
    }

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr out, q;
        slong qn, maglim = ARF_EXP_MAGLIM(prec);
        ulong err = 0;
        int code, tsgn;
        TMP_INIT;

        nk = _arf_elem_limbs(prec, rnd, extra);
        TMP_START;
        out = TMP_ALLOC((nk + 2 + (maglim + FLINT_BITS - 1) / FLINT_BITS + 2) * sizeof(ulong));
        q = out + nk + 2;

        code = _mp_real_pow_fp(out, &err, q, &qn, &tsgn, xd, xn, xe, yd, yn, ye, ARF_SGNBIT(y), _arf_elem_bits(prec, rnd, extra), maglim);

        if (code == MP_REAL_POW_HUGE)
        {
            arf_nan(res);
            status = 0;
        }
        else if (code == MP_REAL_POW_TINY)
            status = _arf_one_plus_tiny(res, rsgn, tsgn, prec, rnd);
        else if (code == MP_REAL_POW_ONE)
        {
            arf_set_si(res, rsgn ? -1 : 1);
            status = 0;
        }
        else
        {
            status = _arf_set_mpn_err(res, out, nk + 1, -FLINT_BITS * nk, ARF_ELEM_EB(err), ARF_ELEM_EB(err), rsgn, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN && qn != 0)
            {
                fmpz_t qz;
                fmpz_init(qz);
                fmpz_set_ui_array(qz, q, qn);
                if (code)
                    fmpz_neg(qz, qz);
                fmpz_add(ARF_EXPREF(res), ARF_EXPREF(res), qz);
                fmpz_clear(qz);
            }
        }

        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
}

int
arf_pow(arf_ptr res, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x) || arf_is_special(y))
    {
        int xs = arf_sgn(x), ys = arf_sgn(y);
        int odd = arf_is_int(y) && !arf_is_int_2exp_si(y, 1);

        /* x^0 = 1 and 1^y = 1 (also for nan) */
        if (arf_is_zero(y) || arf_equal_si(x, 1))
        {
            arf_one(res);
            return 0;
        }

        if (arf_is_nan(x) || arf_is_nan(y))
            arf_nan(res);
        else if (arf_is_zero(x))
        {
            /* 0^y: 0 for y > 0, +inf for y < 0 */
            if (ys > 0)
                arf_zero(res);
            else
                arf_pos_inf(res);
        }
        else if (arf_is_inf(x))
        {
            /* |x|^y with the sign for an odd integer y */
            if (ys > 0)
                (xs < 0 && odd) ? arf_neg_inf(res) : arf_pos_inf(res);
            else
                arf_zero(res);
        }
        else
        {
            /* y = +-inf: by |x| against 1 (1 for |x| = 1) */
            int c = arf_cmpabs_ui(x, 1);
            if (c == 0)
                arf_one(res);
            else if ((c > 0) == (ys > 0))
                arf_pos_inf(res);
            else
                arf_zero(res);
        }
        return 0;
    }

    if (res == x || res == y)
    {
        arf_t tx, ty;
        int status;
        arf_init(tx);
        arf_init(ty);
        arf_set(tx, x);
        arf_set(ty, y);
        status = _arf_pow(res, tx, ty, prec, rnd);
        arf_clear(tx);
        arf_clear(ty);
        return status;
    }

    return _arf_pow(res, x, y, prec, rnd);
}

int
arf_exp10(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    arf_t ten;
    int status;

    /* (exp10 of a nan is nan, unlike 10^0 = 1 for any x) */
    if (arf_is_nan(x))
    {
        arf_nan(res);
        return 0;
    }

    arf_init(ten);
    arf_set_ui(ten, 10);
    status = arf_pow(res, ten, x, prec, rnd);
    arf_clear(ten);
    return status;
}
