/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "arf/elem.h"

/*
    arf_sinh, arf_cosh, arf_sinh_cosh, arf_tanh and the reciprocals
    arf_coth, arf_sech, arf_csch (arf/elem.h explains the rounding).

    - tiny x: sinh x = x (1 + delta), tanh x = x (1 - delta), csch x =
      (1/x)(1 - delta), coth x = (1/x)(1 + delta) with delta < x^2 / 2,
      cosh x = 1 + delta', sech x = 1 - delta' with delta' < x^2 / 2, on
      the known sides.
    - huge x (beyond 2^MP_REAL_EXP_FP_MAX_EXP, up to 2^ARF_EXP_MAGLIM):
      sinh |x| and cosh x are exp(|x|) / 2 within a relative
      2 exp(-2|x|) < 2^-(2^(FLINT_BITS - 2)), far below a unit of any
      working precision (_arf_exp_huge with one more unit on the known
      side), and sech, csch likewise 2 exp(-|x|); tanh and coth are
      +-1 within such a term for |x| >= prec / 2 + 2 already.
    - the register fast paths of arf/elem.h for mantissas of one and two
      limbs at prec <= FLINT_BITS resp. 2 FLINT_BITS (64-bit).
    - the kernels _mp_real_sinh_cosh_fp and _mp_real_tanh_fp at
      increasing guard bits until the rounding is decided; the
      reciprocals by one approximate division (_arf_set_mpn_err_quot).
*/

#define HYP_SINH 0
#define HYP_COSH 1
#define HYP_TANH 2
#define HYP_COTH 3
#define HYP_SECH 4
#define HYP_CSCH 5

/* sinh (rs), cosh (rc), either NULL, or tanh, coth, sech, csch into rt
   (which); x finite and nonzero, not aliased with the outputs */
static int
_arf_hyp(arf_t rs, arf_t rc, arf_t rt, const arf_t x, int which, slong prec, arf_rnd_t rnd)
{
    nn_srcptr d;
    slong xn, e, nk, extra, w, k, pb;
    ulong err;
    fmpz_t xe;
    int xsgn, status, eshift = 0, recip = (which == HYP_COTH || which == HYP_SECH || which == HYP_CSCH);
    int tsgn = (which == HYP_SECH) ? 0 : ARF_SGNBIT(x);

    if (prec == ARF_PREC_EXACT)
        flint_throw(FLINT_ERROR, "arf_sinh_cosh: exact precision not supported\n");

    xsgn = ARF_SGNBIT(x);
    e = ARF_EXP(x);
    ARF_GET_MPN_READONLY(d, xn, x);

    if (FLINT_UNLIKELY(COEFF_IS_MPZ(e)))
    {
        if (fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            /* tiny: evaluated as if the exponent were e (below every
               threshold), the difference restored on the results that
               carry it (sinh, tanh) */
            e = -prec - FLINT_BITS * (xn + 2);
            eshift = 1;
            fmpz_init_set(xe, ARF_EXPREF(x));
        }
        else
            e = WORD_MAX;    /* huge: nan below */
    }

    if (e > MP_REAL_EXP_FP_MAX_EXP)
    {
        if (e > ARF_EXP_MAGLIM(prec))
        {
            if (rs != NULL) arf_nan(rs);
            if (rc != NULL) arf_nan(rc);
            if (rt != NULL) arf_nan(rt);
            return 0;
        }

        switch (which)
        {
            case HYP_TANH:
            case HYP_COTH:
                return _arf_one_plus_tiny(rt, xsgn, which == HYP_TANH, prec, rnd);
            case HYP_SECH:
            case HYP_CSCH:
                /* 2 exp(-|x|) (1 -+ delta) */
                return _arf_exp_huge(rt, d, xn, e, 1, which == HYP_SECH, which == HYP_CSCH, 1, tsgn, prec, rnd);
            default:
                status = 1;
                if (rs != NULL)
                    status = _arf_exp_huge(rs, d, xn, e, 0, 1, 0, -1, xsgn, prec, rnd);
                if (rc != NULL)
                    status = _arf_exp_huge(rc, d, xn, e, 0, 0, 1, -1, 0, prec, rnd);
                return status;
        }
    }

    /* tiny x: the bounds above, delta, delta' < x^2 / 2 < 2^(2e - 1) */
    if (2 * e <= -prec - 2)
    {
        status = 1;
        switch (which)
        {
            case HYP_TANH:
                status = _arf_set_tiny_rel_err(rt, x, -1, 1 - 2 * e, prec, rnd);
                break;
            case HYP_COTH:
                status = _arf_set_tiny_recip_rel_err(rt, x, 1, 1 - 2 * e, prec, rnd);
                break;
            case HYP_SECH:
                status = _arf_one_plus_tiny(rt, 0, 1, prec, rnd);
                break;
            case HYP_CSCH:
                status = _arf_set_tiny_recip_rel_err(rt, x, -1, 1 - 2 * e, prec, rnd);
                break;
            default:
                if (rs != NULL)
                    status = _arf_set_tiny_rel_err(rs, x, 1, 1 - 2 * e, prec, rnd);
                if (status != ARF_ELEM_UNCERTAIN && rc != NULL)
                    status = _arf_one_plus_tiny(rc, 0, 0, prec, rnd);
                break;
        }

        if (status != ARF_ELEM_UNCERTAIN)
        {
            if (eshift)
                fmpz_clear(xe);
            return status;
        }
    }

#if FLINT_BITS == 64
    if (xn <= 2 && prec <= 2 * FLINT_BITS && e >= ((xn == 1 && prec <= FLINT_BITS) ? -31 : -40))
    {
        /* (a one-limb mantissa at more than FLINT_BITS bits padded to two) */
        ulong t3 = 0, t2 = 0, t1 = 0, t0 = 0, c3 = 0, c2 = 0, c1 = 0, c0 = 0;
        ulong d1 = d[xn - 1], d0 = (xn == 2) ? d[0] : 0;
        int one = (xn == 1 && prec <= FLINT_BITS);

        if (which <= HYP_COSH)
        {
            status = 1;
            if (one)
            {
                _arf_small_sinh_cosh_1((rs != NULL) ? &t2 : NULL, &t1, &t0, (rc != NULL) ? &c2 : NULL, &c1, &c0, &k, d[0], e);
                if (rs != NULL)
                    status = _arf_small_set_1(rs, t2, t1, t0, k, ARF_SMALL_HYP_ERR_2, xsgn, prec, rnd);
                if (status != ARF_ELEM_UNCERTAIN && rc != NULL)
                    status = _arf_small_set_1(rc, c2, c1, c0, k, ARF_SMALL_HYP_ERR_2, 0, prec, rnd);
            }
            else
            {
                _arf_small_sinh_cosh_2((rs != NULL) ? &t3 : NULL, &t2, &t1, &t0, (rc != NULL) ? &c3 : NULL, &c2, &c1, &c0, &k, d1, d0, e);
                if (rs != NULL)
                    status = _arf_small_set_2(rs, t3, t2, t1, t0, k, ARF_SMALL_HYP_ERR_3, xsgn, prec, rnd);
                if (status != ARF_ELEM_UNCERTAIN && rc != NULL)
                    status = _arf_small_set_2(rc, c3, c2, c1, c0, k, ARF_SMALL_HYP_ERR_3, 0, prec, rnd);
            }
        }
        else
        {
            int wh = (which == HYP_TANH) ? ARF_SMALL_HYP_TANH : (which == HYP_COTH) ? ARF_SMALL_HYP_COTH
                   : (which == HYP_SECH) ? ARF_SMALL_HYP_SECH : ARF_SMALL_HYP_CSCH;

            if (one)
            {
                if (_arf_small_tanh_1(&t2, &t1, &t0, &k, &err, d[0], e, wh) == 0)
                    status = _arf_small_set_1(rt, t2, t1, t0, k, err, tsgn, prec, rnd);
                else
                    status = ARF_ELEM_UNCERTAIN;
            }
            else
            {
                if (_arf_small_tanh_2(&t3, &t2, &t1, &t0, &k, &err, d1, d0, e, wh) == 0)
                    status = _arf_small_set_2(rt, t3, t2, t1, t0, k, err, tsgn, prec, rnd);
                else
                    status = ARF_ELEM_UNCERTAIN;
            }
        }

        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }
#endif

    /* small x (4e <= -prec - 2) by the series in fixed point (arf/elem.c):
       sinh x = x + x^3/6 within |x|^5 / 120 < 2^(5e - 6), cosh x = 1 +
       x^2/2 within x^4 / 24 < 2^(4e - 4), tanh x = x - x^3/3 within
       2 |x|^5 / 15 < 2^(5e - 2) */
    if (!eshift && e <= -((prec + 5) / 4) && which <= HYP_TANH)
    {
        status = 1;
        if (which == HYP_TANH)
            status = _arf_small_arg_series(rt, d, xn, e, xsgn, ARF_SMALL_SERIES_ATAN, -2, prec, rnd);
        else
        {
            if (rs != NULL)
                status = _arf_small_arg_series(rs, d, xn, e, xsgn, ARF_SMALL_SERIES_SINH, -6, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN && rc != NULL)
                status = _arf_small_arg_series(rc, d, xn, e, 0, ARF_SMALL_SERIES_COSH, -4, prec, rnd);
        }
        if (status != ARF_ELEM_UNCERTAIN)
            return status;
    }

    /* |tanh x| = 1 - 2 / (E^2 + 1) with 2 / (E^2 + 1) < 2 exp(-2 |x|) =
       2^(1 - 2.885 |x|) <= 2^-(prec + 2) (|x| from the top limb in double
       arithmetic, rounded down; and |coth x| = 1 + 2 / (E^2 - 1)
       likewise) */
    if ((which == HYP_TANH || which == HYP_COTH) && e >= 1 && (e > 62
        || 2.8853900817779268 * (1 - 1e-15) * ldexp((double) d[xn - 1], (int) (e - FLINT_BITS)) - 1 >= (double) (prec + 2)))
        return _arf_one_plus_tiny(rt, xsgn, which == HYP_TANH, prec, rnd);

    for (extra = 0; ; extra = (extra == 0) ? 32 : 2 * extra)
    {
        nn_ptr ys, yc;
        slong ylen;
        TMP_INIT;

        pb = _arf_elem_bits(prec, rnd, extra);
        nk = _mp_real_elem_limbs(pb);
        /* the limbs of the kernels: nk + ceil((1 - e) / FLINT_BITS) for
           |x| < 1 (-e <= (prec + 1) / 2 unless a tiny argument fell
           through, then at most prec + FLINT_BITS (xn + 2)) */
        ylen = nk + 2 + ((e < 0) ? (1 - e + FLINT_BITS - 1) / FLINT_BITS : 1);

        TMP_START;
        ys = TMP_ALLOC(2 * ylen * sizeof(ulong));
        yc = ys + ylen;

        if (which == HYP_TANH || which == HYP_COTH)
        {
            _mp_real_tanh_fp(ys, &w, &err, d, xn, e, pb);
            if (which == HYP_TANH)
                status = _arf_set_mpn_err(rt, ys, w + 1, -FLINT_BITS * w, ARF_ELEM_EB(err), ARF_ELEM_EB(err), xsgn, prec, rnd);
            else
                status = _arf_set_mpn_err_quot(rt, ys, w + 1, -FLINT_BITS * w, err, NULL, 0, 0, 1, xsgn, nk, prec, rnd);
        }
        else if (recip)
        {
            k = _mp_real_sinh_cosh_fp((which == HYP_CSCH) ? ys : NULL, (which == HYP_SECH) ? yc : NULL, &w, &err, d, xn, e, pb);
            status = _arf_set_mpn_err_quot(rt, (which == HYP_CSCH) ? ys : yc, w + 1, k - FLINT_BITS * w, err, NULL, 0, 0, 1, tsgn, nk, prec, rnd);
        }
        else
        {
            k = _mp_real_sinh_cosh_fp((rs != NULL) ? ys : NULL, (rc != NULL) ? yc : NULL, &w, &err, d, xn, e, pb);
            status = 1;
            if (rs != NULL)
                status = _arf_set_mpn_err(rs, ys, w + 1, k - FLINT_BITS * w, ARF_ELEM_EB(err), ARF_ELEM_EB(err), xsgn, prec, rnd);
            if (status != ARF_ELEM_UNCERTAIN && rc != NULL)
                status = _arf_set_mpn_err(rc, yc, w + 1, k - FLINT_BITS * w, ARF_ELEM_EB(err), ARF_ELEM_EB(err), 0, prec, rnd);
        }

        TMP_END;

        if (status != ARF_ELEM_UNCERTAIN)
        {
            if (eshift)
            {
                /* the results scaling like x (sinh, tanh) or 1/x (coth,
                   csch) carry the exponent: the fake e restored */
                arf_ptr r = (which == HYP_SINH) ? rs : (which == HYP_COSH || which == HYP_SECH) ? NULL : rt;
                if (r != NULL)
                {
                    fmpz_sub_si(xe, xe, e);
                    if (which == HYP_COTH || which == HYP_CSCH)
                        fmpz_sub(ARF_EXPREF(r), ARF_EXPREF(r), xe);
                    else
                        fmpz_add(ARF_EXPREF(r), ARF_EXPREF(r), xe);
                }
                fmpz_clear(xe);
            }
            return status;
        }
    }
}

/* the special values: sinh, tanh keep 0 and +-inf (tanh: +-1), cosh
   gives 1 and +inf; coth(+-inf) = +-1, sech(0) = 1, sech(+-inf) = 0,
   csch(+-inf) = 0, poles and nan give nan */
static int
_arf_hyp_special(arf_t res, const arf_t x, int which)
{
    if (arf_is_nan(x) || (arf_is_zero(x) && (which == HYP_COTH || which == HYP_CSCH)))
        arf_nan(res);
    else if (arf_is_zero(x))
        arf_set_si(res, (which == HYP_COSH || which == HYP_SECH) ? 1 : 0);
    else
    {
        /* infinite */
        switch (which)
        {
            case HYP_SINH: arf_set(res, x); break;
            case HYP_COSH: arf_pos_inf(res); break;
            case HYP_TANH:
            case HYP_COTH: arf_set_si(res, arf_sgn(x)); break;
            default: arf_zero(res); break;
        }
    }
    return 0;
}

#define HYP_FUNC(name, which) \
int \
name(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd) \
{ \
    if (arf_is_special(x)) \
        return _arf_hyp_special(res, x, which); \
    if (res == x) \
    { \
        arf_t t; \
        int status; \
        arf_init(t); \
        arf_set(t, x); \
        status = ((which) <= HYP_COSH) \
            ? _arf_hyp(((which) == HYP_SINH) ? res : NULL, ((which) == HYP_COSH) ? res : NULL, NULL, t, which, prec, rnd) \
            : _arf_hyp(NULL, NULL, res, t, which, prec, rnd); \
        arf_clear(t); \
        return status; \
    } \
    return ((which) <= HYP_COSH) \
        ? _arf_hyp(((which) == HYP_SINH) ? res : NULL, ((which) == HYP_COSH) ? res : NULL, NULL, x, which, prec, rnd) \
        : _arf_hyp(NULL, NULL, res, x, which, prec, rnd); \
}

HYP_FUNC(arf_sinh, HYP_SINH)
HYP_FUNC(arf_cosh, HYP_COSH)
HYP_FUNC(arf_tanh, HYP_TANH)
HYP_FUNC(arf_coth, HYP_COTH)
HYP_FUNC(arf_sech, HYP_SECH)
HYP_FUNC(arf_csch, HYP_CSCH)

int
arf_sinh_cosh(arf_ptr res1, arf_ptr res2, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        _arf_hyp_special(res1, x, HYP_SINH);
        return _arf_hyp_special(res2, x, HYP_COSH);
    }

    if (res1 == x || res2 == x)
    {
        arf_t t;
        int status;
        arf_init(t);
        arf_set(t, x);
        status = _arf_hyp(res1, res2, NULL, t, HYP_SINH, prec, rnd);
        arf_clear(t);
        return status;
    }

    return _arf_hyp(res1, res2, NULL, x, HYP_SINH, prec, rnd);
}
