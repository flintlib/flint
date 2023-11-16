/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

#define MAGLIM(prec) FLINT_MAX(128, 2 * (prec))

/* todo: min prec by MAG_BITS everywhere? */
void
arb_exp_wide(arb_t res, const arb_t x, slong prec, slong maglim)
{
    mag_t t, u;

    mag_init(t);
    mag_init(u);

    if (arf_cmpabs_2exp_si(arb_midref(x), 20) < 0
        && mag_cmp_2exp_si(arb_radref(x), 20) < 0)
    {
        if (arf_is_zero(arb_midref(x)))
        {
            if (mag_cmp_2exp_si(arb_radref(x), -10) < 0)
            {
                mag_expm1(arb_radref(res), arb_radref(x));
                arf_one(arb_midref(res));
            }
            else
            {
                mag_expinv_lower(t, arb_radref(x));
                mag_exp(u, arb_radref(x));
                arb_set_interval_mag(res, t, u, prec);
            }
        }
        else if (arb_contains_zero(x))
        {
            arf_get_mag_lower(t, arb_midref(x));
            mag_sub(t, arb_radref(x), t);

            arf_get_mag(u, arb_midref(x));
            mag_add(u, arb_radref(x), u);

            if (arf_sgn(arb_midref(x)) > 0)
            {
                mag_expinv_lower(t, t);
                mag_exp(u, u);
                arb_set_interval_mag(res, t, u, prec);
            }
            else
            {
                mag_expinv_lower(u, u);
                mag_exp(t, t);
                arb_set_interval_mag(res, u, t, prec);
            }
        }
        else if (arf_sgn(arb_midref(x)) < 0)
        {
            arb_get_mag(t, x);
            arb_get_mag_lower(u, x);
            mag_expinv_lower(t, t);
            mag_expinv(u, u);
            arb_set_interval_mag(res, t, u, prec);
        }
        else
        {
            arb_get_mag_lower(t, x);
            arb_get_mag(u, x);
            mag_exp_lower(t, t);
            mag_exp(u, u);
            arb_set_interval_mag(res, t, u, prec);
        }
    }
    else
    {
        /* use arb_exp_arf for accurate argument reduction */
        arf_t q;
        arf_init(q);
        arf_set_mag(q, arb_radref(x));
        arf_add(q, q, arb_midref(x), MAG_BITS, ARF_RND_CEIL);
        arb_exp_arf(res, q, FLINT_MIN(prec, MAG_BITS), 0, maglim);
        arb_get_mag(arb_radref(res), res);
        arf_zero(arb_midref(res));
        arf_clear(q);
    }

    mag_clear(t);
    mag_clear(u);
}

void arb_exp(arb_t res, const arb_t x, slong prec)
{
    slong maglim = MAGLIM(prec);

    if (mag_is_zero(arb_radref(x)))
    {
        arb_exp_arf(res, arb_midref(x), prec, 0, maglim);
    }
    else if (mag_is_inf(arb_radref(x)))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(res);
        else
            arb_zero_pm_inf(res);
    }
    else if (arf_is_special(arb_midref(x)))
    {
        if (arf_is_zero(arb_midref(x)))
            arb_exp_wide(res, x, prec, maglim);
        else if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(res);
        else  /* infinity +/- finite */
            arb_exp_arf(res, arb_midref(x), prec, 0, 1);
    }
    else  /* both finite, non-special */
    {
        slong acc, mexp, rexp;

        mexp = ARF_EXP(arb_midref(x));
        rexp = MAG_EXP(arb_radref(x));

        if (COEFF_IS_MPZ(rexp))
            rexp = (fmpz_sgn(ARF_EXPREF(arb_radref(x))) < 0) ? COEFF_MIN : COEFF_MAX;
        if (COEFF_IS_MPZ(mexp))
            mexp = (fmpz_sgn(MAG_EXPREF(arb_midref(x))) < 0) ? COEFF_MIN : COEFF_MAX;

        if (mexp < -prec && rexp < -prec)
        {
            arb_get_mag(arb_radref(res), x);
            mag_expm1(arb_radref(res), arb_radref(res));
            arf_one(arb_midref(res));
            return;
        }

        acc = -rexp;
        acc = FLINT_MAX(acc, 0);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc < 20 && (rexp >= 0 || mexp <= 10))
        {
            /* may evaluate at endpoints */
            arb_exp_wide(res, x, prec, maglim);
        }
        else
        {
            /* exp(a+b) - exp(a) = exp(a) * (exp(b)-1) */
            mag_t t, u;

            mag_init_set(t, arb_radref(x));
            mag_init(u);

            arb_exp_arf(res, arb_midref(x), prec, 0, maglim);
            mag_expm1(t, t);
            arb_get_mag(u, res);
            mag_addmul(arb_radref(res), t, u);

            mag_clear(t);
            mag_clear(u);
        }
    }
}

void
arb_expm1(arb_t res, const arb_t x, slong prec)
{
    slong maglim = MAGLIM(prec);

    if (mag_is_zero(arb_radref(x)))
    {
        arb_exp_arf(res, arb_midref(x), prec, 1, maglim);
    }
    else if (mag_is_inf(arb_radref(x)))
    {
        if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(res);
        else
            arb_zero_pm_inf(res);
    }
    else if (arf_is_special(arb_midref(x)))
    {
        if (arf_is_zero(arb_midref(x)))
        {
            if (mag_cmp_2exp_si(arb_radref(x), -10) < 0)
            {
                mag_expm1(arb_radref(res), arb_radref(x));
                arf_zero(arb_midref(res));
            }
            else
            {
                arb_exp_wide(res, x, prec, maglim);
                arb_sub_ui(res, res, 1, prec);
            }
        }
        else if (arf_is_nan(arb_midref(x)))
            arb_indeterminate(res);
        else  /* infinity +/- finite */
            arb_exp_arf(res, arb_midref(x), prec, 1, 1);
    }
    else  /* both finite, non-special */
    {
        if (arf_cmpabs_2exp_si(arb_midref(x), 3) < 0 &&
            mag_cmp_2exp_si(arb_radref(x), -3) < 0)
        {
            mag_t t, u, one;
            slong acc, mexp, rexp;

            mexp = ARF_EXP(arb_midref(x));
            rexp = MAG_EXP(arb_radref(x));

            if (COEFF_IS_MPZ(rexp))
                rexp = (fmpz_sgn(ARF_EXPREF(arb_radref(x))) < 0) ? COEFF_MIN : COEFF_MAX;
            if (COEFF_IS_MPZ(mexp))
                mexp = (fmpz_sgn(MAG_EXPREF(arb_midref(x))) < 0) ? COEFF_MIN : COEFF_MAX;

            acc = FLINT_MIN(mexp, 0) - rexp;
            acc = FLINT_MAX(acc, 0);
            acc = FLINT_MIN(acc, prec);
            prec = FLINT_MIN(prec, acc + MAG_BITS);
            prec = FLINT_MAX(prec, 2);

            /* [exp(a+b) - 1] - [exp(a) - 1] = exp(a) * (exp(b)-1) */
            mag_init_set(t, arb_radref(x));
            mag_init(u);
            mag_init(one);
            mag_one(one);

            if (arf_sgn(arb_midref(x)) >= 0)
            {
                arb_exp_arf(res, arb_midref(x), prec, 1, maglim);
                arb_get_mag(u, res);
                mag_add(u, u, one);
            }
            else
            {
                arb_exp_arf(res, arb_midref(x), prec, 1, maglim);
                arb_get_mag_lower(u, res);
                mag_sub(u, one, u);
            }

            mag_expm1(t, t);
            mag_addmul(arb_radref(res), t, u);

            mag_clear(t);
            mag_clear(u);
            mag_clear(one);
        }
        else
        {
            arb_exp(res, x, prec);
            arb_sub_ui(res, res, 1, prec);
        }
    }
}

void arb_exp_invexp(arb_t res, arb_t res2, const arb_t x, slong prec)
{
    slong maglim = MAGLIM(prec);

    if (arf_is_special(arb_midref(x)) || mag_is_special(arb_radref(x)))
    {
        /* [c +/- 0] */
        if (arf_is_finite(arb_midref(x)) && mag_is_zero(arb_radref(x)))
        {
            arb_exp_arf(res, arb_midref(x), prec, 0, maglim);
            arb_inv(res2, res, prec);
        }  /* [nan +/- ?] */
        else if (arf_is_nan(arb_midref(x)))
        {
            arb_indeterminate(res);
            arb_indeterminate(res2);
        }  /* [c +/- inf] */
        else if (mag_is_inf(arb_radref(x)))
        {
            arb_zero_pm_inf(res);
            arb_zero_pm_inf(res2);
        }  /* [+inf +/- c] */
        else if (arf_is_pos_inf(arb_midref(x)))
        {
            arb_pos_inf(res);
            arb_zero(res2);
        }  /* [-inf +/- c] */
        else if (arf_is_neg_inf(arb_midref(x)))
        {
            arb_zero(res);
            arb_pos_inf(res2);
        }
        else
        {
            arb_t t;
            arb_init(t);
            arb_neg(t, x);
            arb_exp_wide(res, x, prec, maglim);
            arb_exp_wide(res2, t, prec, maglim);
            arb_clear(t);
        }
    }
    else  /* both finite, non-special */
    {
        slong acc, mexp, rexp;

        mexp = ARF_EXP(arb_midref(x));
        rexp = MAG_EXP(arb_radref(x));

        if (COEFF_IS_MPZ(rexp))
            rexp = (fmpz_sgn(ARF_EXPREF(arb_radref(x))) < 0) ? COEFF_MIN : COEFF_MAX;
        if (COEFF_IS_MPZ(mexp))
            mexp = (fmpz_sgn(MAG_EXPREF(arb_midref(x))) < 0) ? COEFF_MIN : COEFF_MAX;

        if (mexp < -prec && rexp < -prec)
        {
            arb_get_mag(arb_radref(res), x);
            mag_expm1(arb_radref(res), arb_radref(res));
            arf_one(arb_midref(res));
            arb_set(res2, res);
            return;
        }

        acc = -rexp;
        acc = FLINT_MAX(acc, 0);
        acc = FLINT_MIN(acc, prec);
        prec = FLINT_MIN(prec, acc + MAG_BITS);
        prec = FLINT_MAX(prec, 2);

        if (acc < 20 && (rexp >= 0 || mexp <= 10))
        {
            /* may evaluate at endpoints */
            arb_t t;
            arb_init(t);
            arb_neg(t, x);
            arb_exp_wide(res, x, prec, maglim);
            arb_exp_wide(res2, t, prec, maglim);
            arb_clear(t);
        }
        else
        {
            /* exp(a+b) - exp(a) = exp(a) * (exp(b)-1) */
            mag_t t, u;

            mag_init_set(t, arb_radref(x));
            mag_init(u);

            arb_exp_arf(res, arb_midref(x), prec, 0, maglim);
            arb_inv(res2, res, prec);

            mag_expm1(t, t);

            arb_get_mag(u, res);
            mag_addmul(arb_radref(res), t, u);
            arb_get_mag(u, res2);
            mag_addmul(arb_radref(res2), t, u);

            mag_clear(t);
            mag_clear(u);
        }
    }
}
