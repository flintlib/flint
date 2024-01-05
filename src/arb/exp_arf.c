/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "thread_support.h"
#include "mpn_extras.h"

#define TMP_ALLOC_LIMBS(__n) TMP_ALLOC((__n) * sizeof(mp_limb_t))

void
arb_exp_arf_huge(arb_t z, const arf_t x, slong mag, slong prec, int minus_one)
{
    arb_t ln2, t, u;
    fmpz_t q;
    slong wp;

    arb_init(ln2);
    arb_init(t);
    arb_init(u);
    fmpz_init(q);

    wp = prec + mag + 10;

    arb_const_log2(ln2, wp);
    arb_set_arf(t, x);
    arb_div(u, t, ln2, mag + 10);
    arf_get_fmpz(q, arb_midref(u), ARF_RND_DOWN);
    arb_submul_fmpz(t, ln2, q, wp);

    if (minus_one)
    {
        arb_exp(z, t, prec + 10);
        arb_mul_2exp_fmpz(z, z, q);
        arb_sub_ui(z, z, 1, prec);
    }
    else
    {
        arb_exp(z, t, prec);
        arb_mul_2exp_fmpz(z, z, q);
    }

    arb_clear(ln2);
    arb_clear(t);
    arb_clear(u);
    fmpz_clear(q);
}

/* |x| >= 2^expbound */
static void
arb_exp_arf_overflow(arb_t z, slong expbound, int negative, int minus_one, slong prec)
{
    if (!negative)
    {
        arf_zero(arb_midref(z));
        mag_inf(arb_radref(z));
    }
    else
    {
        /* x <= -2^expbound   ==>   0 < exp(x) <= 2^(-2^expbound) */
        fmpz_t t;
        fmpz_init(t);

        fmpz_set_si(t, -1);
        fmpz_mul_2exp(t, t, expbound);

        arf_one(arb_midref(z));
        mag_one(arb_radref(z));
        arb_mul_2exp_fmpz(z, z, t);

        if (minus_one)
            arb_sub_ui(z, z, 1, prec);

        fmpz_clear(t);
    }
}

void
arb_exp_arf_generic(arb_t z, const arf_t x, slong prec, int minus_one)
{
    slong mag;

    mag = arf_abs_bound_lt_2exp_si(x);

    /* reduce by log(2) if needed, but avoid computing log(2) unnecessarily at
       extremely high precision */
    if (mag > 64 || (mag > 8 && prec < 1000000))
    {
        arb_exp_arf_huge(z, x, mag, prec, minus_one);
    }
    else
    {
        int want_rs;
        if (prec < 10000 || mag < -prec / 16)
        {
            want_rs = 1;
        }
        else if (arf_bits(x) < prec / 128)
        {
            want_rs = 0;
        }
        else if (flint_get_num_available_threads() == 1)
        {
            want_rs = (prec < 20000) || (prec < 1000000000 && mag < -prec / 800);
        }
        else
        {
            want_rs = (prec < 10000) || (prec < 1000000000 && mag < -prec / 200);
        }

        if (want_rs)
        {
            arb_exp_arf_rs_generic(z, x, prec, minus_one);
        }
        else
        {
            arb_exp_arf_bb(z, x, prec, minus_one);
        }
    }
}

void
arb_exp_arf(arb_t z, const arf_t x, slong prec, int minus_one, slong maglim)
{
    if (arf_is_special(x))
    {
        if (minus_one)
        {
            if (arf_is_zero(x))
                arb_zero(z);
            else if (arf_is_pos_inf(x))
                arb_pos_inf(z);
            else if (arf_is_neg_inf(x))
                arb_set_si(z, -1);
            else
                arb_indeterminate(z);
        }
        else
        {
            if (arf_is_zero(x))
                arb_one(z);
            else if (arf_is_pos_inf(x))
                arb_pos_inf(z);
            else if (arf_is_neg_inf(x))
                arb_zero(z);
            else
                arb_indeterminate(z);
        }
    }
    else if (COEFF_IS_MPZ(ARF_EXP(x)))
    {
        if (fmpz_sgn(ARF_EXPREF(x)) > 0)
        {
            /* huge input */
            arb_exp_arf_overflow(z, maglim, ARF_SGNBIT(x), minus_one, prec);
        }
        else
        {
            /* |exp(x) - (1 + x)| <= |x^2| */
            fmpz_t t;
            int inexact;
            fmpz_init(t);
            fmpz_mul_2exp(t, ARF_EXPREF(x), 1);
            inexact = arf_add_ui(arb_midref(z), x, minus_one ? 0 : 1, prec, ARB_RND);
            mag_one(arb_radref(z));
            mag_mul_2exp_fmpz(arb_radref(z), arb_radref(z), t);
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
            fmpz_clear(t);
        }
    }
    else
    {
        slong exp, wp, wn, N, r, wprounded, finaln;
        fmpz_t n;
        mp_ptr tmp, w, t, u, finalvalue;
        mp_limb_t p1, q1bits, p2, q2bits, error, error2;
        int negative, inexact;
        TMP_INIT;

        exp = ARF_EXP(x);
        negative = ARF_SGNBIT(x);

        /* handle tiny input */
        /* |exp(x) - 1| <= 2|x| */
        if (!minus_one && exp < -prec - 4)
        {
            arf_one(arb_midref(z));
            mag_set_ui_2exp_si(arb_radref(z), 1, exp + 1);
            return;
        }
        /* |exp(x) - (1 + x)| <= |x^2| */
        else if (exp < (minus_one ? -prec - 4 : -(prec / 2) - 4))
        {
            inexact = arf_add_ui(arb_midref(z), x, minus_one ? 0 : 1, prec, ARB_RND);
            mag_set_ui_2exp_si(arb_radref(z), 1, 2 * exp);
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
            return;
        }

        /* handle huge input */
        if (exp > maglim)
        {
            arb_exp_arf_overflow(z, maglim, negative, minus_one, prec);
            return;
        }

        if (prec >= ARB_EXP_LOG_REDUCTION_PREC && prec <= ARB_LOG_REDUCTION_DEFAULT_MAX_PREC)
        {
            arb_exp_arf_log_reduction(z, x, prec, minus_one);
            return;
        }

        /* Absolute working precision (NOT rounded to a limb multiple) */
        wp = prec + 8;
        if (minus_one && exp <= 0)
            wp += (-exp);
        /* Number of limbs */
        wn = (wp + FLINT_BITS - 1) / FLINT_BITS;
        /* Precision rounded to a number of bits */
        wprounded = FLINT_BITS * wn;
        /* Don't be close to the boundary (to allow adding adding the
           Taylor series truncation error without overflow) */
        wp = FLINT_MAX(wp, wprounded - (FLINT_BITS - 4));

        /* Too high precision to use table -- use generic algorithm */
        if (wp > ARB_EXP_TAB2_PREC)
        {
            arb_exp_arf_generic(z, x, prec, minus_one);
            return;
        }

        TMP_START;

        tmp = TMP_ALLOC_LIMBS(4 * wn + 3);
        w = tmp;        /* requires wn+1 limbs */
        t = w + wn + 1; /* requires wn+1 limbs */
        u = t + wn + 1; /* requires 2wn+1 limbs */

        /* reduce modulo log(2) */
        fmpz_init(n);
        if (_arb_get_mpn_fixed_mod_log2(w, n, &error, x, wn) == 0)
        {
            /* may run out of precision for log(2) */
            arb_exp_arf_generic(z, x, prec, minus_one);
            fmpz_clear(n);
            TMP_END;
            return;
        }

        /* err(w) translates to a propagated error bounded by
           err(w) * exp'(x) < err(w) * exp(1) < err(w) * 3 */
        error *= 3;

        /* Table-based argument reduction (1 or 2 steps) */
        if (wp <= ARB_EXP_TAB1_PREC)
        {
            q1bits = ARB_EXP_TAB1_BITS;
            q2bits = 0;

            p1 = w[wn-1] >> (FLINT_BITS - q1bits);
            w[wn-1] -= (p1 << (FLINT_BITS - q1bits));
            p2 = 0;
        }
        else
        {
            q1bits = ARB_EXP_TAB21_BITS;
            q2bits = ARB_EXP_TAB21_BITS + ARB_EXP_TAB22_BITS;

            p1 = w[wn-1] >> (FLINT_BITS - q1bits);
            w[wn-1] -= (p1 << (FLINT_BITS - q1bits));
            p2 = w[wn-1] >> (FLINT_BITS - q2bits);
            w[wn-1] -= (p2 << (FLINT_BITS - q2bits));
        }

        /* |w| <= 2^-r */
        r = _arb_mpn_leading_zeros(w, wn);

        /* Choose number of terms N such that Taylor series truncation
           error is <= 2^-wp */
        N = _arb_exp_taylor_bound(-r, wp);

        if (N < 60)
        {
            /* Evaluate Taylor series */
            _arb_exp_taylor_rs(t, &error2, w, wn, N);
            /* Taylor series evaluation error */
            error += error2;
            /* Taylor series truncation error */
            error += UWORD(1) << (wprounded-wp);
        }
        else  /* Compute cosh(a) from sinh(a) using a square root. */
        {
            /* the summation for sinh is actually done to (2N-1)! */
            N = (N + 1) / 2;

            /* Evaluate Taylor series for sinh */
            _arb_sin_cos_taylor_rs(t, u, &error2, w, wn, N, 1, 0);
            error += error2;
            error += UWORD(1) << (wprounded-wp);

            /* 1 + sinh^2, with wn + 1 limbs */
            flint_mpn_sqr(u, t, wn);
            u[2 * wn] = 1;

            /* cosh, with wn + 1 limbs */
            mpn_sqrtrem(w, u, u, 2 * wn + 1);

            /* exp = sinh + cosh */
            t[wn] = w[wn] + mpn_add_n(t, t, w, wn);

            /* Error for cosh */
            /* When converting sinh to cosh, the error for cosh must be
               smaller than the error for sinh; but we also get 1 ulp
               extra error from the square root. */
            error2 = error + 1;

            /* Error for sinh + cosh */
            error += error2;
        }

        if (wp <= ARB_EXP_TAB1_PREC)
        {
            if (p1 == 0)
            {
                finalvalue = t;
                finaln = wn + 1;
            }
            else
            {
                /* Divide by 2 to get |t| <= 1 (todo: check this?) */
                mpn_rshift(t, t, wn + 1, 1);
                error = (error >> 1) + 2;

                flint_mpn_mul_n(u, t, arb_exp_tab1[p1] + ARB_EXP_TAB1_LIMBS - wn, wn);

                /* (t + err1 * ulp) * (u + err2 * ulp) + 1ulp = t*u +
                   (err1*u + err2*t + t*u*ulp + 1) * ulp
                   note |u| <= 1, |t| <= 1 */
                error += 4;
                finalvalue = u + wn;

                finaln = wn;

                /* we have effectively divided by 2^2 -- todo use inline function */
                fmpz_add_ui(n, n, 2);
            }
        }
        else
        {
            if (p1 == 0 && p2 == 0)
            {
                finalvalue = t;
                finaln = wn + 1;
            }
            else
            {
                /* Divide by 2 to get |t| <= 1 (todo: check this?) */
                mpn_rshift(t, t, wn + 1, 1);
                error = (error >> 1) + 2;

                flint_mpn_mul_n(u, arb_exp_tab21[p1] + ARB_EXP_TAB2_LIMBS - wn,
                             arb_exp_tab22[p2] + ARB_EXP_TAB2_LIMBS - wn, wn);

                /* error of w <= 4 ulp */
                flint_mpn_copyi(w, u + wn, wn);  /* todo: avoid with better alloc */

                flint_mpn_mul_n(u, t, w, wn);

                /* (t + err1 * ulp) * (w + 4 * ulp) + 1ulp = t*u +
                   (err1*w + 4*t + t*w*ulp + 1) * ulp
                   note |w| <= 1, |t| <= 1 */
                error += 6;

                finalvalue = u + wn;
                finaln = wn;

                /* we have effectively divided by 2^3 -- todo use inline function */
                fmpz_add_ui(n, n, 3);
            }
        }

        /* The accumulated arithmetic error */
        mag_set_ui_2exp_si(arb_radref(z), error, -wprounded);

        /* Set the midpoint */
        if (!minus_one)
        {
            inexact = _arf_set_mpn_fixed(arb_midref(z), finalvalue, finaln, wn, 0, prec, ARB_RND);
            if (inexact)
                arf_mag_add_ulp(arb_radref(z), arb_radref(z), arb_midref(z), prec);
        }
        else
        {
            _arf_set_mpn_fixed(arb_midref(z), finalvalue, finaln, wn, 0, finaln * FLINT_BITS, ARB_RND);
        }

        arb_mul_2exp_fmpz(z, z, n);

        if (minus_one)
            arb_sub_ui(z, z, 1, prec);

        TMP_END;
        fmpz_clear(n);
    }
}
