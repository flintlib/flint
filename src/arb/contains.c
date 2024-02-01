/*
    Copyright (C) 2012-2016, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "arb.h"

int
arb_contains(const arb_t x, const arb_t y)
{
    arf_t t;
    arf_t u;
    arf_t xr, yr;
    arf_struct tmp[4];
    int left_ok, right_ok;

    if (arb_is_exact(y))
        return arb_contains_arf(x, arb_midref(y));

    if (arf_is_nan(arb_midref(y)))
        return arf_is_nan(arb_midref(x));

    arf_init(t);
    arf_init(u);

    arf_init_set_mag_shallow(xr, arb_radref(x));
    arf_init_set_mag_shallow(yr, arb_radref(y));

    /* fast check */
    arf_sub(t, arb_midref(x), xr, MAG_BITS, ARF_RND_CEIL);
    arf_sub(u, arb_midref(y), yr, MAG_BITS, ARF_RND_FLOOR);
    left_ok = arf_cmp(t, u) <= 0;

    /* exact check */
    if (!left_ok)
    {
        arf_init_set_shallow(tmp + 0, arb_midref(x));
        arf_init_neg_mag_shallow(tmp + 1, arb_radref(x));
        arf_init_neg_shallow(tmp + 2, arb_midref(y));
        arf_init_set_mag_shallow(tmp + 3, arb_radref(y));

        arf_sum(t, tmp, 4, MAG_BITS, ARF_RND_DOWN);
        left_ok = arf_sgn(t) <= 0;
    }

    /* fast check */
    arf_add(t, arb_midref(x), xr, MAG_BITS, ARF_RND_FLOOR);
    arf_add(u, arb_midref(y), yr, MAG_BITS, ARF_RND_CEIL);
    right_ok = (arf_cmp(t, u) >= 0);

    /* exact check */
    if (!right_ok)
    {
        arf_init_set_shallow(tmp + 0, arb_midref(x));
        arf_init_set_mag_shallow(tmp + 1, arb_radref(x));
        arf_init_neg_shallow(tmp + 2, arb_midref(y));
        arf_init_neg_mag_shallow(tmp + 3, arb_radref(y));

        arf_sum(t, tmp, 4, MAG_BITS, ARF_RND_DOWN);
        right_ok = arf_sgn(t) >= 0;
    }

    arf_clear(t);
    arf_clear(u);

    return left_ok && right_ok;
}

/* decide |t| <= xr; t = xm - y */
int
arb_contains_arf(const arb_t x, const arf_t y)
{
    if (!arb_is_finite(x))
    {
        if (arf_is_nan(arb_midref(x)))
            return 1;

        if (arf_is_nan(y))
            return 0;

        if (arf_is_inf(arb_midref(x)) && mag_is_finite(arb_radref(x)))
        {
            return arf_equal(arb_midref(x), y);
        }

        return 1;
    }
    else if (!arf_is_finite(y))
    {
        return 0;
    }
    else if (arb_is_exact(x))
    {
        return arf_equal(arb_midref(x), y);
    }
    else
    {
        arf_t t;
        arf_struct tmp[3];
        int result, inexact;

        arf_init(t);

        inexact = arf_sub(t, arb_midref(x), y, 2 * MAG_BITS, ARF_RND_DOWN);

        if (!inexact)
        {
            result = arf_cmpabs_mag(t, arb_radref(x)) <= 0;
        }
        else
        {
            mag_t a;
            mag_init(a);

            arf_get_mag_lower(a, t);
            if (mag_cmp(a, arb_radref(x)) > 0)
            {
                result = 0;
            }
            else
            {
                arf_get_mag(a, t);
                if (mag_cmp(a, arb_radref(x)) < 0)
                {
                    result = 1;
                }
                else
                {
                    /* y >= xm - xr  <=>  0 >= xm - xr - y */
                    arf_init_set_shallow(tmp + 0, arb_midref(x));
                    arf_init_neg_mag_shallow(tmp + 1,  arb_radref(x));
                    arf_init_neg_shallow(tmp + 2, y);

                    arf_sum(t, tmp, 3, MAG_BITS, ARF_RND_DOWN);
                    result = (arf_sgn(t) <= 0);

                    if (result)
                    {
                        /* y <= xm + xr  <=>  0 <= xm + xr - y */
                        arf_init_set_mag_shallow(tmp + 1,  arb_radref(x));
                        arf_sum(t, tmp, 3, MAG_BITS, ARF_RND_DOWN);
                        result = (arf_sgn(t) >= 0);
                    }
                }
            }

            mag_clear(a);
        }

        arf_clear(t);

        return result;
    }
}

int
arb_contains_fmpq(const arb_t x, const fmpq_t y)
{
    if (fmpz_is_one(fmpq_denref(y)) || !arb_is_finite(x))
    {
        return arb_contains_fmpz(x, fmpq_numref(y));
    }
    else
    {
        arf_t t, xm, xr, ym;
        arf_struct tmp[3];
        int result;

        arf_init(t);
        arf_init(xm);
        arf_init(xr);
        arf_init(ym);

        /* To compare x with p/q, compare qx with p. */
        arf_mul_fmpz(xm, arb_midref(x), fmpq_denref(y), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_mag(xr, arb_radref(x));
        arf_mul_fmpz(xr, xr, fmpq_denref(y), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_fmpz(ym, fmpq_numref(y));

        /* y >= xm - xr  <=>  0 >= xm - xr - y */
        arf_init_set_shallow(tmp + 0, xm);
        arf_init_neg_shallow(tmp + 1, xr);
        arf_init_neg_shallow(tmp + 2, ym);

        arf_sum(t, tmp, 3, 30, ARF_RND_DOWN);
        result = (arf_sgn(t) <= 0);

        if (result)
        {
            /* y <= xm + xr  <=>  0 <= xm + xr - y */
            arf_init_set_shallow(tmp + 1, xr);
            arf_sum(t, tmp, 3, 30, ARF_RND_DOWN);
            result = (arf_sgn(t) >= 0);
        }

        arf_clear(t);
        arf_clear(xm);
        arf_clear(xr);
        arf_clear(ym);

        return result;
    }
}

int
arb_contains_fmpz(const arb_t x, const fmpz_t y)
{
    int ans;
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, y);
    ans = arb_contains_arf(x, t);
    arf_clear(t);
    return ans;
}

int
arb_contains_int(const arb_t x)
{
    if (arf_is_int(arb_midref(x)))
    {
        return 1;
    }
    else if (!arb_is_finite(x))
    {
        return arb_contains_zero(x);
    }
    else if (arb_is_exact(x))
    {
        return 0;
    }
    else if (mag_cmp_2exp_si(arb_radref(x), -1) >= 0)  /* radius >= 1/2 */
    {
        return 1;
    }
    else
    {
        /* radius is < 1/2, so it's enough to test the two integers
            bracketing the midpoint */
        arf_t t;
        int res;
        arf_init(t);
        arf_floor(t, arb_midref(x));

        res = arb_contains_arf(x, t);
        if (!res)
        {
            arf_ceil(t, arb_midref(x));
            res = arb_contains_arf(x, t);
        }

        arf_clear(t);
        return res;
    }
}

int
arb_contains_interior(const arb_t x, const arb_t y)
{
    arf_t t;
    arf_t u;
    arf_t xr, yr;
    arf_struct tmp[4];
    int left_ok, right_ok;

    if (arf_is_nan(arb_midref(x)) || arb_is_exact(x) || !arb_is_finite(y))
        return 0;

    arf_init(t);
    arf_init(u);

    arf_init_set_mag_shallow(xr, arb_radref(x));
    arf_init_set_mag_shallow(yr, arb_radref(y));

    /* fast check */
    arf_sub(t, arb_midref(x), xr, MAG_BITS, ARF_RND_CEIL);
    arf_sub(u, arb_midref(y), yr, MAG_BITS, ARF_RND_FLOOR);
    left_ok = arf_cmp(t, u) < 0;

    /* exact check */
    if (!left_ok)
    {
        arf_init_set_shallow(tmp + 0, arb_midref(x));
        arf_init_neg_mag_shallow(tmp + 1, arb_radref(x));
        arf_init_neg_shallow(tmp + 2, arb_midref(y));
        arf_init_set_mag_shallow(tmp + 3, arb_radref(y));

        arf_sum(t, tmp, 4, MAG_BITS, ARF_RND_DOWN);
        left_ok = arf_sgn(t) < 0;
    }

    /* fast check */
    arf_add(t, arb_midref(x), xr, MAG_BITS, ARF_RND_FLOOR);
    arf_add(u, arb_midref(y), yr, MAG_BITS, ARF_RND_CEIL);
    right_ok = (arf_cmp(t, u) > 0);

    /* exact check */
    if (!right_ok)
    {
        arf_init_set_shallow(tmp + 0, arb_midref(x));
        arf_init_set_mag_shallow(tmp + 1, arb_radref(x));
        arf_init_neg_shallow(tmp + 2, arb_midref(y));
        arf_init_neg_mag_shallow(tmp + 3, arb_radref(y));

        arf_sum(t, tmp, 4, MAG_BITS, ARF_RND_DOWN);
        right_ok = arf_sgn(t) > 0;
    }

    arf_clear(t);
    arf_clear(u);

    return left_ok && right_ok;
}

int
arb_contains_mpfr(const arb_t x, const mpfr_t y)
{
    int ans;
    arf_t t;
    arf_init(t);
    arf_set_mpfr(t, y);
    ans = arb_contains_arf(x, t);
    arf_clear(t);
    return ans;
}

int
arb_contains_si(const arb_t x, slong y)
{
    int ans;
    arf_t t;
    arf_init(t);
    arf_set_si(t, y);
    ans = arb_contains_arf(x, t);
    arf_clear(t);
    return ans;
}
