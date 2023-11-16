/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "mpn_extras.h"

int _arb_log_ui_smooth(arb_t res, ulong n, slong prec);

#define TMP_ALLOC_LIMBS(size) TMP_ALLOC((size) * sizeof(mp_limb_t))

#if 0
/* requires x != 1 */
static void
arf_log_via_mpfr(arf_t z, const arf_t x, slong prec, arf_rnd_t rnd)
{
    mpfr_t xf, zf;
    mp_ptr zptr, tmp;
    mp_srcptr xptr;
    mp_size_t xn, zn, val;
    TMP_INIT;
    TMP_START;

    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    tmp = TMP_ALLOC(zn * sizeof(mp_limb_t));

    ARF_GET_MPN_READONLY(xptr, xn, x);

    xf->_mpfr_d = (mp_ptr) xptr;
    xf->_mpfr_prec = xn * FLINT_BITS;
    xf->_mpfr_sign = ARF_SGNBIT(x) ? -1 : 1;
    xf->_mpfr_exp = ARF_EXP(x);

    zf->_mpfr_d = tmp;
    zf->_mpfr_prec = prec;
    zf->_mpfr_sign = 1;
    zf->_mpfr_exp = 0;

    mpfr_set_emin(MPFR_EMIN_MIN);
    mpfr_set_emax(MPFR_EMAX_MAX);

    mpfr_log(zf, xf, arf_rnd_to_mpfr(rnd));

    val = 0;
    while (tmp[val] == 0)
        val++;

    ARF_GET_MPN_WRITE(zptr, zn - val, z);
    flint_mpn_copyi(zptr, tmp + val, zn - val);
    if (zf->_mpfr_sign < 0)
        ARF_NEG(z);

    fmpz_set_si(ARF_EXPREF(z), zf->_mpfr_exp);

    TMP_END;
}
#endif

void
arb_log_arf_huge(arb_t z, const arf_t x, slong prec)
{
    arf_t t;
    arb_t c;
    fmpz_t exp;
    slong wp;

    arf_init(t);
    arb_init(c);
    fmpz_init(exp);

    fmpz_neg(exp, ARF_EXPREF(x));
    arf_mul_2exp_fmpz(t, x, exp);

    wp = prec + 4 - fmpz_bits(exp);
    wp = FLINT_MAX(wp, 4);

    arb_log_arf(z, t, wp);
    arb_const_log2(c, prec + 4);
    arb_submul_fmpz(z, c, exp, prec);

    arf_clear(t);
    arb_clear(c);
    fmpz_clear(exp);
}

void
arb_log_arf(arb_t z, const arf_t x, slong prec)
{
    if (arf_is_special(x))
    {
        if (arf_is_pos_inf(x))
            arb_pos_inf(z);
        else
            arb_indeterminate(z);
    }
    else if (ARF_SGNBIT(x))
    {
        arb_indeterminate(z);
    }
    else if (ARF_IS_POW2(x))
    {
        if (fmpz_is_one(ARF_EXPREF(x)))
        {
            arb_zero(z);
        }
        else
        {
            fmpz_t exp;
            fmpz_init(exp);
            _fmpz_add_fast(exp, ARF_EXPREF(x), -1);
            arb_const_log2(z, prec + 2);
            arb_mul_fmpz(z, z, exp, prec);
            fmpz_clear(exp);
        }
    }
    else if (COEFF_IS_MPZ(ARF_EXP(x)))
    {
        arb_log_arf_huge(z, x, prec);
    }
    else
    {
        slong exp, wp, wn, N, r, closeness_to_one;
        mp_srcptr xp;
        mp_size_t xn, tn;
        mp_ptr tmp, w, t, u;
        mp_limb_t p1, q1bits, p2, q2bits, error, error2, cy;
        int negative, inexact, used_taylor_series;
        TMP_INIT;

        exp = ARF_EXP(x);
        negative = 0;

        ARF_GET_MPN_READONLY(xp, xn, x);

        /* compute a c >= 0 such that |x-1| <= 2^(-c) if c > 0 */
        closeness_to_one = 0;

        if (exp == 0)
        {
            slong i;

            closeness_to_one = FLINT_BITS - FLINT_BIT_COUNT(~xp[xn - 1]);

            if (closeness_to_one == FLINT_BITS)
            {
                for (i = xn - 2; i > 0 && xp[i] == LIMB_ONES; i--)
                    closeness_to_one += FLINT_BITS;

                closeness_to_one += (FLINT_BITS - FLINT_BIT_COUNT(~xp[i]));
            }
        }
        else if (exp == 1)
        {
            closeness_to_one = FLINT_BITS - FLINT_BIT_COUNT(xp[xn - 1] & (~LIMB_TOP));

            if (closeness_to_one == FLINT_BITS)
            {
                slong i;

                for (i = xn - 2; xp[i] == 0; i--)
                    closeness_to_one += FLINT_BITS;

                closeness_to_one += (FLINT_BITS - FLINT_BIT_COUNT(xp[i]));
            }

            closeness_to_one--;
        }

        /* if |t-1| <= 0.5               */
        /* |log(1+t) - t| <= t^2         */
        /* |log(1+t) - (t-t^2/2)| <= t^3 */
        if (closeness_to_one > prec + 1)
        {
            inexact = arf_sub_ui(arb_midref(z), x, 1, prec, ARB_RND);
            mag_set_ui_2exp_si(arb_radref(z), 1, -2 * closeness_to_one);
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
            return;
        }
        else if (2 * closeness_to_one > prec + 1)
        {
            arf_t t, u;
            arf_init(t);
            arf_init(u);
            arf_sub_ui(t, x, 1, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul(u, t, t, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_mul_2exp_si(u, u, -1);
            inexact = arf_sub(arb_midref(z), t, u, prec, ARB_RND);
            mag_set_ui_2exp_si(arb_radref(z), 1, -3 * closeness_to_one);
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
            arf_clear(t);
            arf_clear(u);
            return;
        }

        /* Absolute working precision (NOT rounded to a limb multiple) */
        wp = prec + closeness_to_one + 5;


        /* Too high precision to use table */
        if (wp > ARB_LOG_NEWTON_PREC || wp > ARB_LOG_TAB2_PREC)
        {
            /* Special case: check for smooth integers */
            if (xn == 1 && exp <= FLINT_BITS && exp >= 1)
            {
                ulong n;

                n = xp[0] >> (FLINT_BITS - exp);
                if ((n << (FLINT_BITS - exp)) == xp[0])
                {
                    if (_arb_log_ui_smooth(z, n, prec))
                    {
                        return;
                    }
                }
            }

#if 0
            /* The earlier test for COEFF_IS_MPZ(ARF_EXP(x)) rules out
               too large exponents for MPFR, except on Windows 64 where
               MPFR still uses 32-bit exponents. */
            if (exp < MPFR_EMIN_MIN || exp > MPFR_EMAX_MAX)
            {
                arb_log_arf_huge(z, x, prec);
            }
            else
            {
                arf_log_via_mpfr(arb_midref(z), x, prec, ARB_RND);
                arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
            }
#else
            arb_log_arf_newton(z, x, prec);
#endif
            return;
        }

        /* Working precision in limbs */
        wn = (wp + FLINT_BITS - 1) / FLINT_BITS;

        TMP_START;

        tmp = TMP_ALLOC_LIMBS(4 * wn + 3);
        w = tmp;        /* requires wn+1 limbs */
        t = w + wn + 1; /* requires wn+1 limbs */
        u = t + wn + 1; /* requires 2wn+1 limbs */

        /* read x-1 */
        if (xn <= wn)
        {
            flint_mpn_zero(w, wn - xn);
            mpn_lshift(w + wn - xn, xp, xn, 1);
            error = 0;
        }
        else
        {
            mpn_lshift(w, xp + xn - wn, wn, 1);
            error = 1;
        }

        /* First table-based argument reduction */
        if (wp <= ARB_LOG_TAB1_PREC)
            q1bits = ARB_LOG_TAB11_BITS;
        else
            q1bits = ARB_LOG_TAB21_BITS;

        p1 = w[wn-1] >> (FLINT_BITS - q1bits);

        /* Special case: covers logarithms of small integers */
        if (xn == 1 && (w[wn-1] == (p1 << (FLINT_BITS - q1bits))))
        {
            p2 = 0;
            flint_mpn_zero(t, wn);
            used_taylor_series = 0;
            N = r = 0; /* silence compiler warning */
        }
        else
        {
            /* log(1+w) = log(1+p/q) + log(1 + (qw-p)/(p+q)) */
            w[wn] = mpn_mul_1(w, w, wn, UWORD(1) << q1bits) - p1;
            mpn_divrem_1(w, 0, w, wn + 1, p1 + (UWORD(1) << q1bits));
            error += 1;

            /* Second table-based argument reduction (fused with log->atanh
               conversion) */
            if (wp <= ARB_LOG_TAB1_PREC)
                q2bits = ARB_LOG_TAB11_BITS + ARB_LOG_TAB12_BITS;
            else
                q2bits = ARB_LOG_TAB21_BITS + ARB_LOG_TAB22_BITS;

            p2 = w[wn-1] >> (FLINT_BITS - q2bits);

            u[2 * wn] = mpn_lshift(u + wn, w, wn, q2bits);
            flint_mpn_zero(u, wn);
            flint_mpn_copyi(t, u + wn, wn + 1);
            t[wn] += p2 + (UWORD(1) << (q2bits + 1));
            u[2 * wn] -= p2;
            mpn_tdiv_q(w, u, 2 * wn + 1, t, wn + 1);

            /* propagated error from 1 ulp error: 2 atanh'(1/3) = 2.25 */
            error += 3;

            /* |w| <= 2^-r */
            r = _arb_mpn_leading_zeros(w, wn);

            /* N >= (wp-r)/(2r) */
            N = (wp - r + (2*r-1)) / (2*r);
            N = FLINT_MAX(N, 0);

            /* Evaluate Taylor series */
            _arb_atan_taylor_rs(t, &error2, w, wn, N, 0);
            /* Multiply by 2 */
            mpn_lshift(t, t, wn, 1);
            /* Taylor series evaluation error (multiply by 2) */
            error += error2 * 2;

            used_taylor_series = 1;
        }

        /* Size of output number */
        tn = wn;

        /* First table lookup */
        if (p1 != 0)
        {
            if (wp <= ARB_LOG_TAB1_PREC)
                mpn_add_n(t, t, arb_log_tab11[p1] + ARB_LOG_TAB1_LIMBS - tn, tn);
            else
                mpn_add_n(t, t, arb_log_tab21[p1] + ARB_LOG_TAB2_LIMBS - tn, tn);
            error++;
        }

        /* Second table lookup */
        if (p2 != 0)
        {
            if (wp <= ARB_LOG_TAB1_PREC)
                mpn_add_n(t, t, arb_log_tab12[p2] + ARB_LOG_TAB1_LIMBS - tn, tn);
            else
                mpn_add_n(t, t, arb_log_tab22[p2] + ARB_LOG_TAB2_LIMBS - tn, tn);
            error++;
        }

        /* add exp * log(2) */
        exp--;

        if (exp > 0)
        {
            cy = mpn_addmul_1(t, arb_log_log2_tab + ARB_LOG_TAB2_LIMBS - tn, tn, exp);
            t[tn] = cy;
            tn += (cy != 0);
            error += exp;
        }
        else if (exp < 0)
        {
            t[tn] = 0;
            u[tn] = mpn_mul_1(u, arb_log_log2_tab + ARB_LOG_TAB2_LIMBS - tn, tn, -exp);

            if (mpn_cmp(t, u, tn + 1) >= 0)
            {
                mpn_sub_n(t, t, u, tn + 1);
            }
            else
            {
                mpn_sub_n(t, u, t, tn + 1);
                negative = 1;
            }

            error += (-exp);

            tn += (t[tn] != 0);
        }

        /* The accumulated arithmetic error */
        mag_set_ui_2exp_si(arb_radref(z), error, -wn * FLINT_BITS);

        /* Truncation error from the Taylor series */
        if (used_taylor_series)
            mag_add_ui_2exp_si(arb_radref(z), arb_radref(z), 1, -r*(2*N+1) + 1);

        /* Set the midpoint */
        inexact = _arf_set_mpn_fixed(arb_midref(z), t, tn, wn, negative, prec, ARB_RND);
        if (inexact)
            arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);

        TMP_END;
    }
}

